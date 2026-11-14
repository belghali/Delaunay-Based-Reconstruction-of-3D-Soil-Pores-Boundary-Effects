function statsTD = filtrer_boules_TD_7_frontieres(img, dossier, opts)
% FILTRER_BOULES_TD_7_FRONTIERES
% Etape 5 : filtrage par centre et taux de debordement TD.
%
% Image : 255 (ou 1/logical) = pore ; 0 = solide.
% Les voxels hors du volume sont consideres solides.
% Entree (7 fichiers de l'etape 4) :
%   dossier/boules_tetra_nettoyes/boules_I6.mat, etc.
% Sortie :
%   dossier/filtrage_TD/boules_TD_I6.mat, etc.
%   dossier/filtrage_TD/statistiques_TD.csv et .mat
%
% TD = (# centres de voxels dans la sphere et hors du pore) /
%      (# centres de voxels dans la sphere), y compris hors volume.
% Il s'agit du TD DISCRET sur les centres de voxels, pas du volume
% continu reel d'intersection sphere/solide.
%
% opts.SeuilsTD        = [0 0.05 0.10 0.20] (fractions, pas %)
% opts.CentresAdmis    = 'pore' par defaut ; 'pore_interface' optionnel
% opts.TolRasterRel    = 1e-10 ; tolerance de discretisation d2 <= r^2 + tol
% opts.VerboseTousLes = 5000 ; progression tous les N candidats
%
% Appel :
%   statsTD = filtrer_boules_TD_7_frontieres(img, ...
%                'resultats_delaunay_64');

    if nargin < 2 || isempty(dossier)
        dossier = 'resultats_delaunay_64';
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = defaut(opts, 'SeuilsTD', [0 0.05 0.10 0.20]);
    opts = defaut(opts, 'CentresAdmis', 'pore');
    opts = defaut(opts, 'TolRasterRel', 1e-10);
    opts = defaut(opts, 'VerboseTousLes', 5000);

    validateattributes(img, {'logical','numeric'}, {'3d','real','nonempty'});
    if any(img(:) ~= 0 & img(:) ~= 1 & img(:) ~= 255)
        error('Image binaire attendue : 255/1 = pore ; 0 = solide.');
    end
    if any(~isfinite(img(:)))
        error('L''image contient des valeurs non finies.');
    end
    S = img ~= 0;
    sz = size(S);
    validateattributes(opts.SeuilsTD, {'numeric'}, ...
        {'vector','nonempty','finite','>=',0,'<=',1});
    opts.SeuilsTD = unique(double(opts.SeuilsTD(:).'), 'sorted');
    if ~ismember(opts.CentresAdmis, {'pore','pore_interface'})
        error('opts.CentresAdmis = pore ou pore_interface.');
    end
    validateattributes(opts.TolRasterRel, {'numeric'}, ...
        {'scalar','real','finite','>=',0,'<=',1e-4});
    validateattributes(opts.VerboseTousLes, {'numeric'}, ...
        {'scalar','integer','nonnegative'});

    entree = fullfile(dossier, 'boules_tetra_nettoyes');
    sortie = fullfile(dossier, 'filtrage_TD');
    if ~exist(sortie,'dir'), mkdir(sortie); end
    methodes = {'I6','I18','I26','E6','E18','E26','F'};
    nm = numel(methodes);
    ns = numel(opts.SeuilsTD);

    % Table longue : une ligne par methode et seuil TD.
    Methode = cell(nm*ns,1);
    SeuilTD_pct = zeros(nm*ns,1);
    N_Distinctes = zeros(nm*ns,1);
    N_CentresPore = zeros(nm*ns,1);
    N_CentresInterface = zeros(nm*ns,1);
    N_RejetCentre = zeros(nm*ns,1);
    N_EvalueesTD = zeros(nm*ns,1);
    N_RetenuesTD = zeros(nm*ns,1);
    N_RejetTD = zeros(nm*ns,1);
    RetenuesSurDistinctes_pct = zeros(nm*ns,1);
    RetenuesSurEvaluees_pct = zeros(nm*ns,1);
    TempsTD_s = zeros(nm*ns,1);

    fprintf('\n===== ETAPE 5 : CENTRES STRICTS ET FILTRAGE TD =====\n');
    fprintf('Dimensions : %s | pore = %d voxels\n',mat2str(sz),nnz(S));
    fprintf('Centres admis : %s\n',opts.CentresAdmis);
    fprintf('Seuils TD : %s (fractions)\n',mat2str(opts.SeuilsTD));

    for m = 1:nm
        nom = methodes{m};
        fichier = fullfile(entree, ['boules_' nom '.mat']);
        if ~isfile(fichier), error('Fichier introuvable : %s',fichier); end
        B = load(fichier,'centres','rayons','classeCentres', ...
            'idTetNettoye','idTetOriginal');
        N = numel(B.rayons);
        if size(B.centres,1) ~= N || numel(B.classeCentres) ~= N
            error('%s : donnees incoherentes.',nom);
        end
        nPore = nnz(B.classeCentres == 1);
        nInterface = nnz(B.classeCentres == 2);
        if strcmp(opts.CentresAdmis,'pore')
            garderCentre = (B.classeCentres == 1);
        else
            garderCentre = (B.classeCentres == 1 | B.classeCentres == 2);
        end
        indicesOriginaux = uint32(find(garderCentre));
        centres = double(B.centres(garderCentre,:));
        rayons = double(B.rayons(garderCentre));
        classeCentres = B.classeCentres(garderCentre);
        idTetNettoye = B.idTetNettoye(garderCentre);
        idTetOriginal = B.idTetOriginal(garderCentre);
        M = numel(rayons);
        TD = nan(M,1);
        nVoxSphere = zeros(M,1,'uint32');
        nVoxPore = zeros(M,1,'uint32');
        temps = tic;
        fprintf('\n[%d/7] %s : %d spheres distinctes, %d centres admis\n', ...
            m,nom,N,M);

        for b = 1:M
            c = centres(b,:);
            r = rayons(b);
            if ~(isfinite(r) && r > 0 && all(isfinite(c)))
                error('Sphere non valide %s, indice %d.',nom,b);
            end
            % Marge evitant de manquer un centre de voxel sur la sphere
            % en cas d'arrondi flottant. Le test final reste sur d2.
            marge = 1e-8 * max(1,r);
            lo = ceil(c-r-marge);
            hi = floor(c+r+marge);
            [X,Y,Z] = ndgrid(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3));
            d2 = (X-c(1)).^2 + (Y-c(2)).^2 + (Z-c(3)).^2;
            masqueBoule = d2 <= r*r + opts.TolRasterRel*max(1,r*r);
            ii = X(masqueBoule);
            jj = Y(masqueBoule);
            kk = Z(masqueBoule);
            nTotal = numel(ii);
            if nTotal == 0
                error('Sphere %s/%d sans aucun centre de voxel couvert.',nom,b);
            end
            dansVolume = ii>=1 & ii<=sz(1) & ...
                         jj>=1 & jj<=sz(2) & ...
                         kk>=1 & kk<=sz(3);
            nP = 0;
            if any(dansVolume)
                ids = sub2ind(sz,ii(dansVolume),jj(dansVolume),kk(dansVolume));
                nP = nnz(S(ids));
            end
            TD(b) = (nTotal - nP) / nTotal;
            nVoxSphere(b) = uint32(nTotal);
            nVoxPore(b) = uint32(nP);
            if opts.VerboseTousLes > 0 && mod(b,opts.VerboseTousLes) == 0
                fprintf('   %d / %d | %.1f s\n',b,M,toc(temps));
            end
        end
        t = toc(temps);
        indicesParSeuil = cell(ns,1);
        for s = 1:ns
            ids = find(TD <= opts.SeuilsTD(s) + 1e-12);
            indicesParSeuil{s} = uint32(ids);
            row = (m-1)*ns + s;
            Methode{row} = nom;
            SeuilTD_pct(row) = 100*opts.SeuilsTD(s);
            N_Distinctes(row) = N;
            N_CentresPore(row) = nPore;
            N_CentresInterface(row) = nInterface;
            N_RejetCentre(row) = N-M;
            N_EvalueesTD(row) = M;
            N_RetenuesTD(row) = numel(ids);
            N_RejetTD(row) = M-numel(ids);
            RetenuesSurDistinctes_pct(row) = 100*numel(ids)/max(1,N);
            RetenuesSurEvaluees_pct(row) = 100*numel(ids)/max(1,M);
            TempsTD_s(row) = t;
            fprintf('  TD <= %5.1f%% : %6d / %6d (%.2f%% centres admis)\n', ...
                SeuilTD_pct(row),numel(ids),M,RetenuesSurEvaluees_pct(row));
        end

        metadata = struct('methode',nom,'source',fichier, ...
            'ConventionCentre',opts.CentresAdmis, ...
            'SeuilsTD',opts.SeuilsTD, ...
            'TolRasterRel',opts.TolRasterRel, ...
            'ConventionBords','Tous les centres de voxels hors image sont solides.', ...
            'InterpretationTD','Fraction discrete de centres de voxels hors pore.');
        save(fullfile(sortie,['boules_TD_' nom '.mat']), ...
            'centres','rayons','TD','nVoxSphere','nVoxPore', ...
            'classeCentres','indicesOriginaux','idTetNettoye', ...
            'idTetOriginal','indicesParSeuil','metadata','-v7.3');
        fprintf('  Temps TD total : %.3f s\n',t);
    end

    statsTD = table(Methode,SeuilTD_pct,N_Distinctes, ...
        N_CentresPore,N_CentresInterface,N_RejetCentre, ...
        N_EvalueesTD,N_RetenuesTD,N_RejetTD, ...
        RetenuesSurDistinctes_pct,RetenuesSurEvaluees_pct,TempsTD_s);
    disp(statsTD);
    writetable(statsTD,fullfile(sortie,'statistiques_TD.csv'));
    save(fullfile(sortie,'statistiques_TD.mat'),'statsTD','opts','-v7.3');
    fprintf('\nResultats sauvegardes dans : %s\n',sortie);
end

function opts = defaut(opts,champ,valeur)
    if ~isfield(opts,champ) || isempty(opts.(champ))
        opts.(champ) = valeur;
    end
end
