function [statsAvantTC, statsComparaison] = couverture_union_avant_TC_7_frontieres(img, dossier, opts)
% COUVERTURE_UNION_AVANT_TC_7_FRONTIERES
%
% Pour les 7 methodes et 4 seuils TD volumiques, mesure la couverture
% obtenue avec TOUTES les boules admissibles, SANS filtrage TC.
% Compare ensuite cette couverture aux 112 selections TC precedentes.
%
% Important :
%   - 255 ou 1 = PORE ; 0 = SOLIDE.
%   - TD = 0 : utilise EstInscrite (controle geometrique deterministe).
%   - TD = 5, 10, 20 % : utilise TDvolumeGrille, PAS TD discret.
%   - Toutes les mesures d'union sont faites sur la meme grille reguliere
%     de k^3 sous-voxels par voxel que le programme precedent.
%   - Aucun Monte Carlo, aucun nouveau Delaunay, aucun nouveau TD.
%   - Il s'agit d'une approximation volumique deterministe sur grille.
%
% Fichiers d'entree :
%   <dossier>/selection_TC_volumique_grille_complete/
%       banque_TDvolume_I6.mat ... banque_TDvolume_F.mat
%       statistiques_TD_TC.mat
%
% Fichiers de sortie (sans ecraser les precedents) :
%   <dossier>/selection_TC_volumique_grille_complete/couverture_avant_TC/
%       statistiques_avant_TC.csv
%       comparaison_avant_apres_TC.csv
%       resultats_avant_TC.mat
%       union_<methode>_TDxx.mat (si SauverMasques = true)
%
% Exemple :
%   opts = struct();
%   opts.SousVoxelsParAxe = 4;
%   opts.SauverMasques = false;
%   [statsAvantTC, statsComparaison] = ...
%       couverture_union_avant_TC_7_frontieres( ...
%           img, 'resultats_delaunay_64', opts);
%   disp(statsAvantTC);
%   disp(statsComparaison);

    if nargin < 2 || isempty(dossier)
        dossier = 'resultats_delaunay_64';
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    opts = defaut(opts, 'SousVoxelsParAxe', 4);
    opts = defaut(opts, 'SeuilsTD', [0 0.05 0.10 0.20]);
    opts = defaut(opts, 'VerboseTousLes', 5000);
    opts = defaut(opts, 'SauverMasques', false);
    opts = defaut(opts, 'MaxPointsGrille', 8e7);

    validateattributes(img, {'numeric','logical'}, ...
        {'real','nonempty','3d'});
    if any(~isfinite(double(img(:)))) || ...
            any(img(:)~=0 & img(:)~=1 & img(:)~=255)
        error('Image binaire attendue : 255/1 = pore, 0 = solide.');
    end
    validateattributes(opts.SousVoxelsParAxe, {'numeric'}, ...
        {'scalar','integer','positive'});
    validateattributes(opts.SeuilsTD, {'numeric'}, ...
        {'vector','nonempty','finite','>=',0,'<=',1});
    validateattributes(opts.VerboseTousLes, {'numeric'}, ...
        {'scalar','integer','nonnegative'});

    S = (img ~= 0);
    sz = double(size(S));
    k = double(opts.SousVoxelsParAxe);
    seuilsTD = unique(double(opts.SeuilsTD(:).'), 'sorted');
    methodes = {'I6','I18','I26','E6','E18','E26','F'};

    dossierTC = fullfile(dossier, 'selection_TC_volumique_grille_complete');
    sortie = fullfile(dossierTC, 'couverture_avant_TC');

    if ~isfolder(dossierTC)
        error('Dossier TD/TC inexistant : %s', dossierTC);
    end
    if ~isfolder(sortie)
        mkdir(sortie);
    end

    %% 1) Lire les banques TD existantes pour retrouver la meme marge.
    maxR = 0;
    for im = 1:numel(methodes)
        f = fullfile(dossierTC, ...
            ['banque_TDvolume_' methodes{im} '.mat']);
        if ~isfile(f)
            error('Banque TD manquante : %s', f);
        end
        D = load(f, 'rayons', 'opts');
        if ~isfield(D,'rayons') || isempty(D.rayons)
            error('Rayons absents ou vides : %s', f);
        end
        if ~isfield(D,'opts') || ...
                ~isfield(D.opts,'SousVoxelsParAxe')
            error('Resolution d''origine non enregistree : %s', f);
        end
        if D.opts.SousVoxelsParAxe ~= k
            error(['Grille incompatible pour %s : TD/TC calcule avec ' ...
                   'k=%d, mais k=%d a ete demande.'], ...
                   methodes{im}, D.opts.SousVoxelsParAxe, k);
        end
        maxR = max(maxR, max(double(D.rayons(:))));
    end

    % Reconstituer exactement le repere de la grille precedente.
    marge = ceil(maxR) + 2;
    taillePad = sz + 2*marge;
    tailleHi = taillePad * k;
    nbPtsHi = prod(tailleHi);

    if nbPtsHi > opts.MaxPointsGrille
        error(['Grille trop grande : %.0f sous-voxels. ' ...
               'Cette version est prevue pour un petit volume/une ROI.'], ...
               nbPtsHi);
    end

    Sp = false(taillePad);
    Sp((1:sz(1))+marge, ...
       (1:sz(2))+marge, ...
       (1:sz(3))+marge) = S;
    SpHi = repelem(Sp, k, k, k);
    clear Sp

    nPoreHi = nnz(SpHi);
    if nPoreHi == 0
        error('Le volume ne contient aucun voxel poreux.');
    end

    %% 2) Calculer 28 unions, sans selection TC.
    nLignes = numel(methodes)*numel(seuilsTD);
    noms = cell(nLignes,1);
    mesures = nan(nLignes,13);
    ligne = 0;

    fprintf('\n===== COUVERTURE AVANT TC : UNION DE TOUTES LES BOULES =====\n');
    fprintf('Pore : %d voxels | k = %d | seuils TD : %s\n', ...
        nnz(S), k, mat2str(seuilsTD));
    fprintf('Aucun Monte Carlo et aucun recalcul du TD.\n');

    for im = 1:numel(methodes)
        nom = methodes{im};
        f = fullfile(dossierTC, ['banque_TDvolume_' nom '.mat']);
        B = load(f, 'centres', 'rayons', 'TDvolumeGrille', 'EstInscrite');

        champs = {'centres','rayons','TDvolumeGrille','EstInscrite'};
        if ~all(isfield(B,champs))
            error('%s : banque TD incomplete.', nom);
        end

        centres = double(B.centres);
        rayons = double(B.rayons(:));
        tauxTD = double(B.TDvolumeGrille(:));
        inscrites = logical(B.EstInscrite(:));
        N = numel(rayons);

        if size(centres,1)~=N || size(centres,2)~=3 || ...
                numel(tauxTD)~=N || numel(inscrites)~=N || ...
                any(~isfinite(centres(:))) || any(~isfinite(rayons)) || ...
                any(rayons<=0) || any(~isfinite(tauxTD)) || ...
                any(tauxTD<0 | tauxTD>1)
            error('%s : donnees incompatibles ou invalides.', nom);
        end

        fprintf('\n[%d/7] %s : %d spheres candidates\n', im, nom, N);

        % Les seuils sont croissants. On ajoute uniquement les nouvelles
        % boules. Si les ensembles ne sont pas imbriques (cas limite lie
        % au controle geometrique TD=0), on reconstruit l'union.
        unionBoules = false(tailleHi);
        dejaAjoutees = false(N,1);
        tMethode = tic;

        for it = 1:numel(seuilsTD)
            td = seuilsTD(it);
            if td == 0
                admises = inscrites;
            else
                admises = (tauxTD <= td + 1e-12);
            end

            reconstruire = any(dejaAjoutees & ~admises);
            if reconstruire
                unionBoules(:) = false;
                dejaAjoutees(:) = false;
            end

            aAjouter = find(admises & ~dejaAjoutees);
            tAjout = tic;

            for j = 1:numel(aAjouter)
                b = aAjouter(j);
                lin = indices_sphere_sousvox( ...
                    centres(b,:), rayons(b), tailleHi, marge, k);
                if isempty(lin)
                    error('%s : sphere %d sans point sur la grille.', nom, b);
                end
                unionBoules(lin) = true;

                if opts.VerboseTousLes > 0 && ...
                        mod(j,opts.VerboseTousLes)==0
                    fprintf('   TD %g %% : %d/%d nouvelles boules | %.1f s\n', ...
                        td*100, j, numel(aAjouter), toc(tAjout));
                end
            end
            tempsAjout = toc(tAjout);
            dejaAjoutees = admises;

            % Mesures sur les memes sous-voxels que le calcul TC.
            TP = nnz(unionBoules & SpHi);
            FP = nnz(unionBoules & ~SpHi);
            totalUnion = TP+FP;
            FN = nPoreHi-TP;

            couverture_pct = 100*TP/nPoreHi;
            if totalUnion > 0
                debordUnion_pct = 100*FP/totalUnion;
                precision_pct = 100*TP/totalUnion;
            else
                debordUnion_pct = 0;  % Convention : union vide.
                precision_pct = NaN;  % Precision non definie si union vide.
            end
            iou_pct = 100*TP/(nPoreHi+FP);

            ligne = ligne+1;
            noms{ligne} = nom;
            mesures(ligne,:) = [ ...
                100*td, N, nnz(admises), numel(aAjouter), ...
                couverture_pct, debordUnion_pct, precision_pct, ...
                iou_pct, TP/k^3, FP/k^3, FN/k^3, ...
                tempsAjout, toc(tMethode)];

            fprintf(['   TD %5.1f %% | %6d boules | ' ...
                     'couverture %6.2f %% | debordement union %5.2f %% | ' ...
                     'IoU %6.2f %%\n'], ...
                     100*td, nnz(admises), couverture_pct, ...
                     debordUnion_pct, iou_pct);

            if opts.SauverMasques
                masqueUnion = unionBoules; %#ok<NASGU>
                meta = struct('Methode',nom,'SeuilTD',td, ...
                    'SousVoxelsParAxe',k,'Marge',marge, ...
                    'TailleImage',sz,'TailleGrille',tailleHi);
                fmasque = fullfile(sortie, sprintf( ...
                    'union_%s_TD%02d.mat', nom, round(100*td)));
                save(fmasque, 'masqueUnion','meta','-v7.3');
            end
        end
    end

    statsAvantTC = table( ...
        noms, mesures(:,1), mesures(:,2), mesures(:,3), ...
        mesures(:,4), mesures(:,5), mesures(:,6), ...
        mesures(:,7), mesures(:,8), mesures(:,9), ...
        mesures(:,10), mesures(:,11), mesures(:,12), ...
        mesures(:,13), ...
        'VariableNames', { ...
        'Methode','SeuilTD_pct','BoulesCandidates', ...
        'BoulesAdmissibles','NouvellesBoulesAjoutees', ...
        'CouverturePore_pct','DebordementUnion_pct', ...
        'Precision_pct','IoU_pct', ...
        'VolPoreCouvert_vox3','VolHorsPore_vox3', ...
        'VolPoreNonCouvert_vox3','TempsAjout_s', ...
        'TempsCumuleMethode_s'});

    %% 3) Comparaison avec les 112 selections TC deja calculees.
    cheminStats = fullfile(dossierTC, 'statistiques_TD_TC.mat');
    statsComparaison = table();

    if ~isfile(cheminStats)
        warning(['Statistiques TC absentes : comparaison non produite. ' ...
                 'Les 28 unions avant TC sont bien calculees.']);
    else
        X = load(cheminStats, 'statsTC');
        if ~isfield(X,'statsTC')
            error('La variable statsTC manque dans %s.', cheminStats);
        end
        statsComparaison = X.statsTC;
        n = height(statsComparaison);

        statsComparaison.CouvertureAvantTC_pct = nan(n,1);
        statsComparaison.CouverturePerdue_points = nan(n,1);
        statsComparaison.CouvertureConservee_pct = nan(n,1);
        statsComparaison.ReductionBoules_pct = nan(n,1);
        statsComparaison.DebordementAvantTC_pct = nan(n,1);
        statsComparaison.IoUAvantTC_pct = nan(n,1);

        for j = 1:n
            nom = statsComparaison.Methode{j};
            td = statsComparaison.SeuilTD_pct(j);
            ix = find(strcmp(statsAvantTC.Methode,nom) & ...
                abs(statsAvantTC.SeuilTD_pct-td)<1e-8, 1);

            if isempty(ix)
                error('Ligne avant TC absente : %s, TD=%g %%.', nom, td);
            end
            nAdmis = statsAvantTC.BoulesAdmissibles(ix);
            if nAdmis ~= statsComparaison.AdmissiblesTD(j)
                error(['Nombre de boules different entre banques TD ' ...
                       'et statistiques TC : %s, TD=%g %%.'], nom, td);
            end

            cAvant = statsAvantTC.CouverturePore_pct(ix);
            cApres = statsComparaison.CouverturePore_pct(j);
            statsComparaison.CouvertureAvantTC_pct(j) = cAvant;
            statsComparaison.CouverturePerdue_points(j) = ...
                cAvant-cApres;
            if cAvant > 0
                statsComparaison.CouvertureConservee_pct(j) = ...
                    100*cApres/cAvant;
            end
            if nAdmis > 0
                statsComparaison.ReductionBoules_pct(j) = ...
                    100*(1-statsComparaison.BoulesSelectionnees(j)/nAdmis);
            end
            statsComparaison.DebordementAvantTC_pct(j) = ...
                statsAvantTC.DebordementUnion_pct(ix);
            statsComparaison.IoUAvantTC_pct(j) = ...
                statsAvantTC.IoU_pct(ix);

            if cApres > cAvant+1e-6
                warning(['Couverture apres TC superieure a celle de toutes ' ...
                         'les boules : %s TD=%g TC=%g. Verifier les grilles.'], ...
                    nom, td, statsComparaison.SeuilTC_pct(j));
            end
        end
    end

    %% 4) Sauvegardes CSV et MAT.
    writetable(statsAvantTC, ...
        fullfile(sortie,'statistiques_avant_TC.csv'));
    if ~isempty(statsComparaison)
        writetable(statsComparaison, ...
            fullfile(sortie,'comparaison_avant_apres_TC.csv'));
    end
    save(fullfile(sortie,'resultats_avant_TC.mat'), ...
        'statsAvantTC','statsComparaison','opts','-v7.3');

    fprintf('\n===== TERMINE =====\n');
    fprintf('%d unions avant TC calculees.\n',height(statsAvantTC));
    if ~isempty(statsComparaison)
        fprintf('%d comparaisons avec TC calculees.\n', ...
            height(statsComparaison));
    end
    fprintf('Resultats : %s\n',sortie);
end

function idx = indices_sphere_sousvox(c,r,szHi,marge,k)
% Exactement la meme conversion de coordonnees que dans la routine TD/TC.
% Centre de sous-voxel d'indice p : 0.5-marge+(p-0.5)/k.
    lo = ceil(k*(c-r-0.5+marge)+0.5);
    hi = floor(k*(c+r-0.5+marge)+0.5);
    if any(lo<1) || any(hi>szHi)
        error('Marge insuffisante pour une sphere de rayon %.6f.',r);
    end
    if any(lo>hi)
        idx = zeros(0,1);
        return;
    end

    dx = (lo(1):hi(1))';
    dy = (lo(2):hi(2))';
    dz = (lo(3):hi(3))';
    x = 0.5-marge+(dx-0.5)/k-c(1);
    y = 0.5-marge+(dy-0.5)/k-c(2);
    z = 0.5-marge+(dz-0.5)/k-c(3);

    d2 = reshape(x.^2,[],1,1) + ...
         reshape(y.^2,1,[],1) + ...
         reshape(z.^2,1,1,[]);
    local = find(d2 <= r*r+1e-12*max(1,r*r));
    if isempty(local)
        idx = zeros(0,1);
        return;
    end
    [ii,jj,kk] = ind2sub(size(d2),local);
    idx = sub2ind(szHi,lo(1)-1+ii,lo(2)-1+jj,lo(3)-1+kk);
end

function o = defaut(o,nom,valeur)
    if ~isfield(o,nom) || isempty(o.(nom))
        o.(nom) = valeur;
    end
end
