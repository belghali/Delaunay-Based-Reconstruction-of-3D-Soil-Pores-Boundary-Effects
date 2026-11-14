function [statsTC, statsTD] = selection_TC_volumique_7_frontieres(img, dossier, opts)
% SELECTION_TC_VOLUMIQUE_7_FRONTIERES
% Benchmark des 7 definitions de frontiere : TD volumique puis TC volumique.
%
% ENTREE :
%   img  : volume binaire 3D ; 255 (ou 1) = pore, 0 = solide.
%   dossier : dossier contenant boules_tetra_nettoyes/boules_<M>.mat
%   opts.SeuilsTD : [0 .05 .10 .20] (fractions)
%   opts.SeuilsTC : [0 .05 .10 .20] (fractions)
%   opts.NMonteCarloTD : 8192 tirages par boule pour TD > 0
%   opts.Graine : 42, tirages identiques pour les 7 methodes
%   opts.SousVoxelsTC : 3 sous-voxels par axe (4 conseille pour validation)
%   opts.ToleranceGeometrique : 1e-10
%   opts.VerboseTousLes : 5000 ; 0 = silencieux
%   opts.SauverSelections : true
%
% IMPORTANT :
%  * AUCUN PRE-FILTRAGE PAR TD DISCRET ! La fonction repart des spheres
%    distinctes brutes ayant leur centre strictement dans le pore.
%  * TD=0 : certification de non-penetration dans les cubes solides et
%    de non-depassement du domaine, a tolerance numerique pres.
%  * TD>0 : FRACTION VOLUMIQUE ESTIMEE PAR MONTE CARLO (non exacte).
%  * TC : fraction volumique de la sphere deja couverte par l'union des
%    boules choisies, approximee sur une grille sous-voxelique reguliere.
%    Denominateur = volume echantillonne de toute la sphere, y compris
%    la portion dans le solide ou hors des limites de l'image.
%  * On trie les candidats par rayon decroissant ; egalites departagees
%    par (i,j,k) croissants pour garantir la reproductibilite.
%  * Le domaine est prolonge par une marge solide suffisante pour inclure
%    toutes les spheres admissibles. Le debordement de l'union inclut
%    donc aussi les portions hors du cube initial.
%
% SORTIES :
%  statsTC : une ligne pour chaque methode x TD x TC (112 lignes).
%  statsTD : une ligne par methode x TD (28 lignes).
%  dossier/selection_TC_volumique/statistiques_TC.csv, statistiques_TD.csv
%  dossier/selection_TC_volumique/selection_<M>_TDxx_TCxx.mat
%  dossier/selection_TC_volumique/banque_TDvolume_<M>.mat
%
% USAGE :
%  opts = struct('NMonteCarloTD',8192,'SousVoxelsTC',3);
%  [statsTC,statsTD] = selection_TC_volumique_7_frontieres( ...
%      img,'resultats_delaunay_64',opts);
%
% Pour publication, realiser une sensibilite a NMonteCarloTD et
% SousVoxelsTC : le TD positif et le TC sont des approximations numeriques.

    if nargin < 2 || isempty(dossier), dossier = 'resultats_delaunay_64'; end
    if nargin < 3 || isempty(opts), opts = struct(); end
    opts = opt_defaut(opts,'SeuilsTD',[0 .05 .10 .20]);
    opts = opt_defaut(opts,'SeuilsTC',[0 .05 .10 .20]);
    opts = opt_defaut(opts,'NMonteCarloTD',8192);
    opts = opt_defaut(opts,'Graine',42);
    opts = opt_defaut(opts,'SousVoxelsTC',3);
    opts = opt_defaut(opts,'ToleranceGeometrique',1e-10);
    opts = opt_defaut(opts,'VerboseTousLes',5000);
    opts = opt_defaut(opts,'SauverSelections',true);
    opts = opt_defaut(opts,'MaxPointsGrille',8e7);

    validateattributes(img,{'numeric','logical'},{'real','nonempty','3d'});
    if any(~isfinite(double(img(:)))) || any(img(:)~=0 & img(:)~=1 & img(:)~=255)
        error('img doit etre binaire : 255 ou 1 = pore, 0 = solide.');
    end
    S = (img ~= 0);
    sz = double(size(S));
    validateattributes(opts.SeuilsTD,{'numeric'},{'vector','nonempty','finite','>=',0,'<=',1});
    validateattributes(opts.SeuilsTC,{'numeric'},{'vector','nonempty','finite','>=',0,'<=',1});
    validateattributes(opts.NMonteCarloTD,{'numeric'},{'scalar','integer','>=',100});
    validateattributes(opts.SousVoxelsTC,{'numeric'},{'scalar','integer','>=',1});
    validateattributes(opts.Graine,{'numeric'},{'scalar','integer','nonnegative'});
    validateattributes(opts.ToleranceGeometrique,{'numeric'}, ...
        {'scalar','finite','nonnegative','<=',1e-5});
    validateattributes(opts.VerboseTousLes,{'numeric'},{'scalar','integer','nonnegative'});
    opts.SeuilsTD = unique(double(opts.SeuilsTD(:).'),'sorted');
    opts.SeuilsTC = unique(double(opts.SeuilsTC(:).'),'sorted');
    k = double(opts.SousVoxelsTC);
    nm = 7;
    methodes = {'I6','I18','I26','E6','E18','E26','F'};
    entree = fullfile(dossier,'boules_tetra_nettoyes');
    sortie = fullfile(dossier,'selection_TC_volumique');
    if ~isfolder(entree), error('Dossier inexistant : %s',entree); end
    if ~isfolder(sortie), mkdir(sortie); end

    % Calculer une marge de solide commune aux 7 banques.
    maxR = 0;
    for m = 1:nm
        f = fullfile(entree,['boules_' methodes{m} '.mat']);
        if ~isfile(f), error('Fichier absent : %s',f); end
        D = load(f,'rayons','classeCentres');
        if ~all(isfield(D,{'rayons','classeCentres'}))
            error('Le fichier %s ne contient pas rayons/classeCentres.',f);
        end
        r = double(D.rayons(D.classeCentres == 1));
        if ~isempty(r), maxR = max(maxR,max(r)); end
    end
    marge = ceil(maxR)+2;
    SzPad = sz+2*marge;
    szHi = k*SzPad;
    if prod(szHi) > opts.MaxPointsGrille
        error(['Grille trop grande : %.0f sous-voxels. ' ...
               'Reduire SousVoxelsTC, employer une ROI ou une version par blocs.'],prod(szHi));
    end
    Sp = false(SzPad);
    Sp((1:sz(1))+marge,(1:sz(2))+marge,(1:sz(3))+marge) = S;
    SpHi = repelem(Sp,k,k,k);
    nPoreHi = nnz(SpHi);
    clear Sp

    % Tirages communs, uniformes en volume dans une boule unite.
    etat = rng;
    restauration = onCleanup(@() rng(etat)); %#ok<NASGU>
    rng(opts.Graine,'twister');
    a = rand(opts.NMonteCarloTD,3);
    z = 2*a(:,1)-1;
    theta = 2*pi*a(:,2);
    rho = a(:,3).^(1/3);
    v = rho.*sqrt(max(0,1-z.^2));
    U = [v.*cos(theta),v.*sin(theta),rho.*z];
    clear a z theta rho v

    nomsTD = {}; stTD = [];
    nomsTC = {}; stTC = [];
    fprintf('\n===== TD VOLUMIQUE + TC VOLUMIQUE : 7 METHODES =====\n');
    fprintf('Pore=%d voxels | TD-MC=%d points/boule | TC=%d sous-voxels/axe\n', ...
        nnz(S),size(U,1),k);
    fprintf('Seuils TD : %s | seuils TC : %s\n', ...
        mat2str(opts.SeuilsTD),mat2str(opts.SeuilsTC));

    for m = 1:nm
        nom = methodes{m};
        f = fullfile(entree,['boules_' nom '.mat']);
        B = load(f,'centres','rayons','classeCentres');
        if ~all(isfield(B,{'centres','rayons','classeCentres'}))
            error('%s : contenu incomplet.',nom);
        end
        idsPore = find(B.classeCentres(:)==1);
        centres = double(B.centres(idsPore,:));
        rayons = double(B.rayons(idsPore));
        N = numel(rayons);
        if size(centres,1)~=N || size(centres,2)~=3 || ...
                any(~isfinite(centres(:))) || any(~isfinite(rayons)) || any(rayons<=0)
            error('%s : coordonnees/rayons invalides.',nom);
        end
        fprintf('\n[%d/7] %s : %d centres strictement poreux\n',m,nom,N);

        % Auditer TOUS les centres poreux de la banque brute, sans aucun
        % filtre TD discret. TD continu estime par MC, + test exact TD=0.
        TDvolume_MC = zeros(N,1);
        EstInscrite = false(N,1);
        ttd = tic;
        for b = 1:N
            c = centres(b,:); r = rayons(b);
            X = bsxfun(@plus, U*r, c);
            V = floor(X+0.5);
            dans = V(:,1)>=1 & V(:,1)<=sz(1) & ...
                   V(:,2)>=1 & V(:,2)<=sz(2) & ...
                   V(:,3)>=1 & V(:,3)<=sz(3);
            nP = 0;
            if any(dans)
                iv = sub2ind(sz,V(dans,1),V(dans,2),V(dans,3));
                nP = nnz(S(iv));
            end
            TDvolume_MC(b) = (size(U,1)-nP)/size(U,1);
            EstInscrite(b) = sphere_inscrite_exacte(S,sz,c,r,opts.ToleranceGeometrique);
            if opts.VerboseTousLes>0 && mod(b,opts.VerboseTousLes)==0
                fprintf(' TD : %d/%d | %.1f s\n',b,N,toc(ttd));
            end
        end
        TempsTD_s = toc(ttd);
        if opts.SauverSelections
            provenance = 'Centres strictement poreux ; aucun prefiltrage TD discret';
            save(fullfile(sortie,['banque_TDvolume_' nom '.mat']), ...
                'centres','rayons','idsPore','TDvolume_MC','EstInscrite', ...
                'provenance','opts','-v7.3');
        end
        fprintf(' TD volumique moyen estime : %.4f %% | inscrites : %d/%d\n', ...
            100*mean(TDvolume_MC),nnz(EstInscrite),N);

        for it = 1:numel(opts.SeuilsTD)
            td = opts.SeuilsTD(it);
            if td == 0
                admissible = EstInscrite;
            else
                admissible = (TDvolume_MC <= td);
            end
            ix = find(admissible);
            nAdmissibles = numel(ix);
            nomsTD{end+1,1} = nom; %#ok<AGROW>
            stTD(end+1,:) = [100*td,N,nAdmissibles, ...
                N-nAdmissibles,TempsTD_s]; %#ok<AGROW>
            fprintf(' TD %5.1f %% : %6d admissibles\n',100*td,nAdmissibles);

            % Tri deterministe : rayon decroissant puis i,j,k croissants.
            if nAdmissibles>0
                [~,ordre] = sortrows([-rayons(ix),centres(ix,:)], [1 2 3 4]);
                ix = ix(ordre);
            end

            for jt = 1:numel(opts.SeuilsTC)
                tc = opts.SeuilsTC(jt);
                ttc = tic;
                couverture = false(szHi);
                sel = zeros(nAdmissibles,1,'uint32');
                chevRetenus = zeros(nAdmissibles,1);
                nbSel = 0;
                nAucun = 0;
                for b = 1:nAdmissibles
                    ii = ix(b);
                    lin = indices_sphere_sousvox(centres(ii,:),rayons(ii), ...
                        szHi,marge,k);
                    if isempty(lin), nAucun=nAucun+1; continue; end
                    % TC = Vol(B intersection union precedentes) / Vol(B).
                    % Tous les points de la sphere sont pris en compte,
                    % y compris dans la marge solide.
                    overlap = nnz(couverture(lin))/numel(lin);
                    if overlap <= tc + 1e-12
                        nbSel = nbSel+1;
                        sel(nbSel) = uint32(ii);
                        chevRetenus(nbSel) = overlap;
                        couverture(lin) = true;
                    end
                    if opts.VerboseTousLes>0 && mod(b,opts.VerboseTousLes)==0
                        fprintf('   TD=%g%% TC=%g%% : %d/%d | retenues=%d\n', ...
                            100*td,100*tc,b,nAdmissibles,nbSel);
                    end
                end
                sel = double(sel(1:nbSel));
                chevRetenus = chevRetenus(1:nbSel);
                TempsTC_s = toc(ttc);
                TP = nnz(couverture & SpHi);
                FP = nnz(couverture & ~SpHi);
                denomU = TP+FP;
                % Tous les volumes sont en sous-voxels de meme volume.
                couverturePore_pct = 100*TP/max(1,nPoreHi);
                debordUnion_pct = 100*FP/max(1,denomU);
                precision_pct = 100*TP/max(1,denomU);
                iou_pct = 100*TP/max(1,nPoreHi+FP);
                nomsTC{end+1,1} = nom; %#ok<AGROW>
                stTC(end+1,:) = [100*td,100*tc,N,nAdmissibles,nbSel, ...
                    couverturePore_pct,debordUnion_pct,precision_pct, ...
                    iou_pct,TP/k^3,FP/k^3,TempsTD_s,TempsTC_s]; %#ok<AGROW>
                fprintf('    TC %5.1f %% : %5d boules | couverture %.2f %% | ' ...
                    'TD union %.2f %% | IoU %.2f %% | %.2f s\n', ...
                    100*tc,nbSel,couverturePore_pct,debordUnion_pct,iou_pct,TempsTC_s);
                if opts.SauverSelections
                    centresSelectionnes = centres(sel,:);
                    rayonsSelectionnes = rayons(sel);
                    indicesBanqueBrute = uint32(idsPore(sel));
                    TDVolumeSelection = TDvolume_MC(sel);
                    TCSelection = chevRetenus;
                    meta = struct('Methode',nom,'SeuilTD',td,'SeuilTC',tc, ...
                        'Origine',f,'NMonteCarloTD',opts.NMonteCarloTD, ...
                        'SousVoxelsTC',k,'Graine',opts.Graine, ...
                        'TDzero','Controle geometrique deterministe', ...
                        'TDpositif','Monte Carlo dans la boule', ...
                        'TC','Union sous-voxelique dans boite prolongee');
                    fsave = fullfile(sortie,sprintf( ...
                        'selection_%s_TD%02d_TC%02d.mat',nom,round(100*td),round(100*tc)));
                    save(fsave,'centresSelectionnes','rayonsSelectionnes', ...
                        'indicesBanqueBrute','TDVolumeSelection', ...
                        'TCSelection','meta','-v7.3');
                end
            end
        end
    end
    statsTD = table(nomsTD,stTD(:,1),stTD(:,2),stTD(:,3), ...
        stTD(:,4),stTD(:,5),'VariableNames', ...
        {'Methode','SeuilTD_pct','CentresPore','AdmissiblesTD', ...
         'RejeteesTD','TempsAuditTD_s'});
    statsTC = table(nomsTC,stTC(:,1),stTC(:,2),stTC(:,3),stTC(:,4), ...
        stTC(:,5),stTC(:,6),stTC(:,7),stTC(:,8),stTC(:,9), ...
        stTC(:,10),stTC(:,11),stTC(:,12),stTC(:,13), ...
        'VariableNames',{'Methode','SeuilTD_pct','SeuilTC_pct', ...
        'CentresPore','AdmissiblesTD','BoulesSelectionnees', ...
        'CouverturePore_pct','DebordementUnion_pct','Precision_pct', ...
        'IoU_pct','VolPoreCouvert_vox3','VolHorsPore_vox3', ...
        'TempsAuditTD_s','TempsTC_s'});
    writetable(statsTD,fullfile(sortie,'statistiques_TD.csv'));
    writetable(statsTC,fullfile(sortie,'statistiques_TC.csv'));
    save(fullfile(sortie,'statistiques_TD_TC.mat'), ...
        'statsTD','statsTC','opts','-v7.3');
    fprintf('\nTermine : %d configurations. Dossier : %s\n',height(statsTC),sortie);
end

function tf = sphere_inscrite_exacte(S,sz,c,r,tol)
% Test deterministe au sens du modele voxel = cube de cote 1.
    d = tol*max(1,r);
    if any(c-r < 0.5-d) || any(c+r > sz+0.5+d)
        tf = false; return;
    end
    lo = max([1 1 1],ceil(c-r-0.5));
    hi = min(sz,floor(c+r+0.5));
    tf = true;
    if any(lo>hi), return; end
    [I,J,K] = ndgrid(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3));
    ids = sub2ind(sz,I(:),J(:),K(:));
    mask = ~S(ids);
    if ~any(mask), return; end
    Q = double([I(mask),J(mask),K(mask)]);
    D = max(abs(bsxfun(@minus,Q,c))-0.5,0);
    dist2 = sum(D.^2,2);
    tf = ~any(dist2 < r*r - tol*max(1,r*r));
end

function idx = indices_sphere_sousvox(c,r,szHi,marge,k)
% Sous-voxel global p : x(p)=0.5-marge+(p-0.5)/k.
    lo = ceil(k*(c-r-0.5+marge)+0.5);
    hi = floor(k*(c+r-0.5+marge)+0.5);
    if any(lo<1) || any(hi>szHi)
        error('Marge insuffisante pour la sphere de rayon %.6f.',r);
    end
    if any(lo>hi), idx=zeros(0,1); return; end
    dx = (lo(1):hi(1))';
    dy = (lo(2):hi(2))';
    dz = (lo(3):hi(3))';
    x = 0.5-marge+(dx-0.5)/k-c(1);
    y = 0.5-marge+(dy-0.5)/k-c(2);
    z = 0.5-marge+(dz-0.5)/k-c(3);
    d2 = reshape(x.^2,[],1,1)+reshape(y.^2,1,[],1)+reshape(z.^2,1,1,[]);
    local = find(d2 <= r*r + 1e-12*max(1,r*r));
    if isempty(local), idx=zeros(0,1); return; end
    [ii,jj,kk] = ind2sub(size(d2),local);
    idx = sub2ind(szHi,lo(1)-1+ii,lo(2)-1+jj,lo(3)-1+kk);
end

function o = opt_defaut(o,n,v)
    if ~isfield(o,n) || isempty(o.(n)), o.(n)=v; end
end
