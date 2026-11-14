function Tbudget = budget_fixe_7_frontieres(img, dossierBenchmark, budgets, dossierSortie, nomDataset)
% BUDGET_FIXE_7_FRONTIERES - Comparaison IoU a budget fixe pour 7 frontieres.
%
% Adaptation de la PARTIE 3 (budget fixe) du script original
% analyser_tore_128_3_points.m, retrouve dans les travaux IMPROVE 2027.
%
% INVARIANTS DU PROTOCOLE ORIGINAL :
%   - banques de boules deja calculees (AUCUN nouveau Delaunay),
%   - audit TD deja calcule, seuil TD <= 0.10,
%   - aucun filtrage TC / chevauchement,
%   - ordre : rayon decroissant, puis centre (x,y,z) croissant,
%   - accepter B si au moins un NOUVEAU voxel de pore y est contenu,
%   - evaluer IoU de l'union COMPLETE des spheres (y compris hors volume),
%   - ne jamais extrapoler un budget non atteint.
%
% APPEL :
%   T = budget_fixe_7_frontieres(img, ...
%       'BENCHMARK_512_NEUF', [1000 5000 10000 25000 50000], ...
%       'BENCHMARK_512_NEUF/BUDGET_FIXE_TD10', 'BD12');
%
% ENTREE img : volume binaire 3D ; 1 ou 255 = PORE / FOREGROUND,
%              0 = SOLIDE. Pour DRP, passer (bin==0), pas bin brut.
%
% Fichiers requis dans dossierBenchmark :
%   boules_tetra_nettoyes/boules_I6.mat ... boules_F.mat
%   TD_TC_step1_sequentiel_audit_k1/audit_TD_I6.mat ... audit_TD_F.mat
%
% SORTIES : CSV, MAT, figure PNG. MATLAB avec tableaux 3D.
% Aucun calcul LBM, de topologie ou de distances n'est effectue.

    if nargin < 2, error('Fournir img et dossierBenchmark.'); end
    if nargin < 3 || isempty(budgets), budgets = [100 250 500 1000 2500]; end
    if nargin < 4 || isempty(dossierSortie)
        dossierSortie = fullfile(dossierBenchmark, 'BUDGET_FIXE_TD10');
    end
    if nargin < 5 || isempty(nomDataset), nomDataset = 'Dataset'; end

    assert(isnumeric(img) || islogical(img), 'Image numerique/logique requise.');
    assert(ndims(img)==3, 'Image 3D requise.');
    assert(all(img(:)==0 | img(:)==1 | img(:)==255), ...
        'Image attendue : pore=1/255, solide=0.');
    ref = logical(img ~= 0);
    clear img
    sz = double(size(ref));
    nPore = nnz(ref);
    assert(nPore>0 && nPore<numel(ref), 'Image sans pore ou sans solide.');

    budgets = unique(double(budgets(:)'));
    assert(~isempty(budgets) && all(isfinite(budgets)) && ...
        all(budgets>0) && all(budgets==floor(budgets)), ...
        'Budgets : entiers strictement positifs.');

    repBoules = fullfile(dossierBenchmark,'boules_tetra_nettoyes');
    repAudit  = fullfile(dossierBenchmark,'TD_TC_step1_sequentiel_audit_k1');
    assert(isfolder(repBoules) && isfolder(repAudit), ...
        'Dossiers de banque / audits TD introuvables dans %s.',dossierBenchmark);
    if ~isfolder(dossierSortie), mkdir(dossierSortie); end
    methodes = {'I6','I18','I26','E6','E18','E26','F'};
    tdPct = 10;
    lignes = struct([]);

    fprintf('\n=== BUDGET FIXE / %s / TD10, SANS TC ===\n',nomDataset);
    fprintf('Dimensions %s, pores %d (%.5f%%)\n',mat2str(sz),nPore,100*nPore/prod(sz));
    fprintf('Dossier benchmark : %s\n',dossierBenchmark);
    fprintf('Budgets : %s\n',mat2str(budgets));

    for m=1:numel(methodes)
        nom = methodes{m}; tMethode = tic;
        fB = fullfile(repBoules,['boules_' nom '.mat']);
        fAudit = fullfile(repAudit,['audit_TD_' nom '.mat']);
        assert(isfile(fB), 'Banque de spheres introuvable : %s',fB);
        assert(isfile(fAudit), 'Audit TD introuvable : %s',fAudit);

        % Meme association identifiants / audit TD que le code du tore.
        B = load(fB,'centres','rayons','classeCentres');
        D = load(fAudit,'TD','idsBruts');
        assert(all(isfield(B,{'centres','rayons','classeCentres'})), ...
            'Champs manquants dans %s',fB);
        assert(all(isfield(D,{'TD','idsBruts'})), ...
            'Champs manquants dans %s',fAudit);
        ids = find(B.classeCentres(:)==1);
        assert(isequal(double(ids(:)),double(D.idsBruts(:))), ...
            'Incoherence IDs/audit TD pour %s.',nom);
        C  = double(B.centres(ids,:));
        RR = double(B.rayons(ids)); RR=RR(:);
        TD = double(D.TD(:));
        clear B D ids
        assert(size(C,2)==3 && size(C,1)==numel(RR) && ...
            numel(RR)==numel(TD), 'Dimensions incompatibles pour %s.',nom);

        bons = find(isfinite(TD) & TD <= tdPct/100+1e-12);
        clear TD
        if ~isempty(bons)
            [~,o] = sortrows([-RR(bons), C(bons,:)], [1 2 3 4]);
            bons = bons(o);
            clear o
        end

        % Meme convention que le benchmark d'origine : tout voxel hors de
        % l'image est SOLIDE VIRTUEL. Le padding evite de tronquer les boules.
        [padL,padH] = padding_admissibles(C(bons,:),RR(bons),sz);
        szPad = sz+padL+padH;
        if prod(szPad)>2^31
            error('%s : domaine etendu trop grand (%s).',nom,mat2str(szPad));
        end
        fprintf('\n[%s] admissibles TD10 = %d | padding [%s ; %s] | grille %s\n', ...
            nom,numel(bons),mat2str(padL),mat2str(padH),mat2str(szPad));

        refPad = false(szPad);
        refPad((1:sz(1))+padL(1),(1:sz(2))+padL(2), ...
               (1:sz(3))+padL(3)) = ref;
        U = false(szPad);
        TP = 0; FP = 0; nChoisi = 0; k = 1;

        for j=1:numel(bons)
            b = bons(j);
            ind = indices_sphere_padded(C(b,:),RR(b),szPad,padL);
            if isempty(ind), continue; end

            % Equivalent EXACT du test du code original du tore :
            % any(ref(ind) & ~U(ind)). Aucun seuil de chevauchement.
            libres = ~U(ind);
            if ~any(refPad(ind) & libres), continue; end

            nv = ind(libres);
            np = nnz(refPad(nv));
            TP = TP + np;
            FP = FP + numel(nv) - np;
            U(nv) = true;
            nChoisi = nChoisi+1;

            if nChoisi == budgets(k)
                IoU = 100*TP/(nPore+FP);
                couverture = 100*TP/nPore;
                debordement = 100*FP/(TP+FP);
                ligne = struct('Dataset',string(nomDataset), ...
                    'Methode',string(nom),'TD_pct',tdPct, ...
                    'Budget',budgets(k),'NBoules',nChoisi, ...
                    'IoU_pct',IoU,'Couverture_pct',couverture, ...
                    'Debordement_pct',debordement, ...
                    'VoxelsPoreCouverts',TP,'VoxelsHorsPore',FP, ...
                    'BudgetAtteint',true);
                if isempty(lignes), lignes=ligne;
                else, lignes(end+1)=ligne; end %#ok<AGROW>
                fprintf('%-4s N=%6d : IoU=%7.3f%% | couverture=%7.3f%% | dehors=%d\n', ...
                    nom,nChoisi,IoU,couverture,FP);
                k=k+1;
                if k>numel(budgets), break; end
            end
            if mod(j,100000)==0
                fprintf('  %s : candidats examines %d/%d, spheres contributives=%d\n', ...
                    nom,j,numel(bons),nChoisi);
            end
        end

        % Les budgets impossibles restent NaN : jamais de remplissage fictif.
        while k<=numel(budgets)
            ligne=struct('Dataset',string(nomDataset), ...
                'Methode',string(nom),'TD_pct',tdPct, ...
                'Budget',budgets(k),'NBoules',nChoisi, ...
                'IoU_pct',NaN,'Couverture_pct',NaN, ...
                'Debordement_pct',NaN, ...
                'VoxelsPoreCouverts',TP,'VoxelsHorsPore',FP, ...
                'BudgetAtteint',false);
            if isempty(lignes), lignes=ligne;
            else, lignes(end+1)=ligne; end %#ok<AGROW>
            warning('%s : budget %d non atteint (seulement %d spheres contributives).', ...
                nom,budgets(k),nChoisi);
            k=k+1;
        end
        fprintf('%s termine en %.1f s.\n',nom,toc(tMethode));
        clear C RR bons refPad U
    end

    Tbudget = struct2table(lignes);
    cheminCSV = fullfile(dossierSortie,'02_comparaison_budget_identique.csv');
    writetable(Tbudget,cheminCSV);
    save(fullfile(dossierSortie,'02_comparaison_budget_identique.mat'), ...
         'Tbudget','budgets','tdPct','nomDataset','-v7.3');
    tracer_budget(Tbudget,nomDataset,dossierSortie);
    fprintf('\nTermine : %s\n',cheminCSV);
end

function [padL,padH] = padding_admissibles(C,R,sz)
    % Identique au main_benchmark_7_frontieres_step1_sequentiel_audit.m.
    if isempty(R), padL=[0 0 0];padH=[0 0 0];return; end
    mn=min(ceil(C-R-1e-12.*max(1,R)),[],1);
    mx=max(floor(C+R+1e-12.*max(1,R)),[],1);
    padL=double(max(0,1-mn));
    padH=double(max(0,mx-sz));
end

function lin = indices_sphere_padded(c,r,szPad,padL)
    % Recopie de l'enumerateur du benchmark : voxel inclus si
    % ||centre_voxel - centre_sphere||^2 <= r^2 + tolerance.
    tol=1e-12*max(1,r*r);
    lo=max([1 1 1],ceil(c-r+padL-1e-12*max(1,r)));
    hi=min(szPad,floor(c+r+padL+1e-12*max(1,r)));
    if any(lo>hi), lin=zeros(0,1);return;end
    x=(lo(1):hi(1))-padL(1)-c(1);
    y=(lo(2):hi(2))-padL(2)-c(2);
    z=(lo(3):hi(3))-padL(3)-c(3);
    if prod(hi-lo+1)<=350000
        d2=reshape(x.^2,[],1,1)+reshape(y.^2,1,[],1)+ ...
            reshape(z.^2,1,1,[]);
        mask=d2<=r*r+tol;
        [ii,jj,kk]=ind2sub(hi-lo+1,find(mask));
        lin=sub2ind(szPad,lo(1)-1+ii,lo(2)-1+jj,lo(3)-1+kk);
    else
        [xx,yy]=ndgrid(lo(1):hi(1),lo(2):hi(2));
        dxy=(xx-padL(1)-c(1)).^2+(yy-padL(2)-c(2)).^2;
        lin=zeros(0,1);
        for k=lo(3):hi(3)
            mask=dxy+(k-padL(3)-c(3))^2<=r*r+tol;
            v=find(mask);
            if isempty(v),continue;end
            [ii,jj]=ind2sub(size(xx),v);
            it=sub2ind(szPad,lo(1)-1+ii,lo(2)-1+jj, ...
                repmat(k,numel(ii),1));
            lin=[lin;it]; %#ok<AGROW>
        end
    end
end

function tracer_budget(T,nomDataset,dossierSortie)
    f=figure('Visible','off','Color','w','Position',[100 100 1000 600]);
    hold on;
    noms=unique(T.Methode,'stable');
    for i=1:numel(noms)
        z=T(T.Methode==noms(i) & T.BudgetAtteint,:);
        if isempty(z),continue;end
        plot(z.Budget,z.IoU_pct,'-o','LineWidth',1.5, ...
            'MarkerSize',5,'DisplayName',char(noms(i)));
    end
    xlabel('Number of selected spheres'); ylabel('IoU (%)');
    title(sprintf('Fixed-sphere-budget comparison - %s (TD = 10%%, no TC)',nomDataset));
    ylim([0 100]); grid on;
    legend('Location','best');
    try
        exportgraphics(f,fullfile(dossierSortie,'figure_02_iou_budget_fixe.png'), ...
            'Resolution',300);
    catch
        saveas(f,fullfile(dossierSortie,'figure_02_iou_budget_fixe.png'));
    end
    close(f);
end
