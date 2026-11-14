function statsN = nettoyer_tetra_7_frontieres(img, dossier, opts)
% NETTOYER_TETRA_7_FRONTIERES  Nettoyage des tetraedres pour 7 frontieres.
%
%   statsN = nettoyer_tetra_7_frontieres(img, dossier)
%   statsN = nettoyer_tetra_7_frontieres(img, dossier, opts)
%
% ENTREES
%   img     volume binaire 3D : 255 (ou 1/true) = PORE, 0 = solide.
%   dossier dossier de l'etape 2 contenant delaunay_I6.mat ...
%   opts.SeuilQualite      : 1e-3 par defaut (qualite mean-ratio).
%   opts.SeuilDetRelatif   : 1e-12 par defaut (degenerescence numerique).
%   opts.ToleranceFace    : 1e-10 voxel par defaut.
%   opts.GarderInterface  : true par defaut (centroide sur interface).
%   opts.TailleBloc       : 100000 tetraedres par defaut.
%
% METHODE
%   1) Supprimer tetraedres de volume nul / numeriquement degenere.
%   2) Supprimer les tetraedres q < SeuilQualite.
%   3) Examiner le centroide de chaque tetraedre restant.
%      - PORE : centroide dans l'union de voxels poreux, pas sur interface.
%      - INTERFACE : centroide au contact exact de pore et solide.
%      - SOLIDE : centroide uniquement dans le solide.
%      - HORS : centroide en dehors du domaine voxelise.
%   Le centre d'un voxel (i,j,k) est (i,j,k); son cube occupe
%   [i-0.5,i+0.5] x [j-0.5,j+0.5] x [k-0.5,k+0.5].
%   L'exterieur du volume est traite comme du solide virtuel.
%
% IMPORTANT
%   Le critere de centroide ne garantit PAS que tout le tetraedre soit
%   inclus dans le pore. C'est un filtre d'admissibilite, pas un test
%   exact d'intersection tetraedre / solide.
%   Ne pas filtrer les sommets : les methodes E* et F ont volontairement
%   des sommets en dehors du pore ou sur son interface.
%
% SORTIES
%   dossier/nettoyes/delaunay_nettoye_<methode>.mat
%      points, tetra, qualite, tetra_ids_origine, classe_centroides, metadata
%   dossier/nettoyes/statistiques_nettoyage_tetra.csv / .mat
%
% EXEMPLE
%   statsN = nettoyer_tetra_7_frontieres(img,'resultats_delaunay_64');

    if nargin < 2 || isempty(dossier)
        dossier = 'resultats_delaunay_64';
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = defaut(opts,'SeuilQualite',1e-3);
    opts = defaut(opts,'SeuilDetRelatif',1e-12);
    opts = defaut(opts,'ToleranceFace',1e-10);
    opts = defaut(opts,'GarderInterface',true);
    opts = defaut(opts,'TailleBloc',100000);

    validateattributes(img,{'numeric','logical'}, ...
        {'3d','nonempty','real'},mfilename,'img');
    if any(img(:) ~= 0 & img(:) ~= 1 & img(:) ~= 255)
        error('Image attendue : 255/1 = pore ; 0 = solide.');
    end
    validateattributes(opts.SeuilQualite,{'numeric'}, ...
        {'scalar','real','finite','>=',0,'<=',1});
    validateattributes(opts.SeuilDetRelatif,{'numeric'}, ...
        {'scalar','real','finite','>=',0});
    validateattributes(opts.ToleranceFace,{'numeric'}, ...
        {'scalar','real','finite','>=',0,'<',1e-4});
    validateattributes(opts.TailleBloc,{'numeric'}, ...
        {'scalar','integer','positive'});
    if ~isscalar(opts.GarderInterface) || ~islogical(logical(opts.GarderInterface))
        error('opts.GarderInterface doit etre un booleen scalaire.');
    end
    opts.GarderInterface = logical(opts.GarderInterface);

    S = (img ~= 0);
    sz = size(S);
    methodes = {'I6','I18','I26','E6','E18','E26','F'};
    nm = numel(methodes);

    sortie = fullfile(dossier,'nettoyes');
    if ~exist(sortie,'dir'), mkdir(sortie); end

    NombrePoints = zeros(nm,1);
    TetraInitiaux = zeros(nm,1);
    Degeneres = zeros(nm,1);
    PresquePlats = zeros(nm,1);
    CentroideDansPore = zeros(nm,1);
    CentroideSurInterface = zeros(nm,1);
    CentroideDansSolide = zeros(nm,1);
    CentroideHorsDomaine = zeros(nm,1);
    TetraConserves = zeros(nm,1);
    TauxConservation_pct = zeros(nm,1);
    QualiteMoyenneFinale = nan(nm,1);
    QualiteMedianeFinale = nan(nm,1);
    TempsNettoyage_s = zeros(nm,1);

    fprintf('\n===== ETAPE 3 : NETTOYAGE DES TETRAEDRES =====\n');
    fprintf('Pore : %d voxels, dimensions %s\n',nnz(S),mat2str(sz));
    fprintf('Seuil qualite q < %.3g, determinant relatif <= %.3g\n', ...
        opts.SeuilQualite,opts.SeuilDetRelatif);
    fprintf('Centroide sur interface conserve : %d\n',opts.GarderInterface);

    for m = 1:nm
        nom = methodes{m};
        fichier = fullfile(dossier,['delaunay_' nom '.mat']);
        if ~isfile(fichier)
            error('Impossible de trouver : %s',fichier);
        end
        data = load(fichier,'points','tetra','qualite');
        if ~isfield(data,'qualite')
            error('Champ qualite absent dans %s.',fichier);
        end
        points = double(data.points);
        tetra0 = data.tetra;
        q0 = double(data.qualite(:));
        clear data

        N = size(tetra0,1);
        if size(points,2) ~= 3 || size(tetra0,2) ~= 4 || numel(q0) ~= N
            error('%s : dimensions incoherentes.',nom);
        end
        NombrePoints(m) = size(points,1);
        TetraInitiaux(m) = N;
        fprintf('\n[%d/%d] %s : %d tetraedres initiaux\n',m,nm,nom,N);
        horloge = tic;

        % Stockage des indices originaux seulement, afin de ne pas dupliquer
        % inutilement l'ensemble des tetraedres pendant les boucles.
        idsGardes = zeros(N,1,'uint32');
        classesGardees = zeros(N,1,'uint8');
        cumulGardes = 0;

        for debut = 1:opts.TailleBloc:N
            finBloc = min(N,debut+opts.TailleBloc-1);
            ids = (debut:finBloc).';
            T = double(tetra0(ids,:));
            A = points(T(:,1),:);
            B = points(T(:,2),:);
            C = points(T(:,3),:);
            D = points(T(:,4),:);

            AB = B-A; AC = C-A; AD = D-A;
            BC = C-B; BD = D-B; CD = D-C;
            v6 = abs(sum(AB .* cross(AC,AD,2),2));
            L2max = max([sum(AB.^2,2),sum(AC.^2,2), ...
                         sum(AD.^2,2),sum(BC.^2,2), ...
                         sum(BD.^2,2),sum(CD.^2,2)],[],2);
            detRel = v6 ./ max(L2max.^(3/2),realmin);
            degen = ~isfinite(detRel) | ...
                    detRel <= opts.SeuilDetRelatif;
            Degeneres(m) = Degeneres(m) + nnz(degen);

            qb = q0(ids);
            plat = ~degen & (~isfinite(qb) | qb < opts.SeuilQualite);
            PresquePlats(m) = PresquePlats(m) + nnz(plat);

            okGeometrie = ~degen & ~plat;
            if ~any(okGeometrie), continue; end
            G = (A(okGeometrie,:) + B(okGeometrie,:) + ...
                 C(okGeometrie,:) + D(okGeometrie,:)) / 4;
            classe = classer_centroides_pore(S,G,opts.ToleranceFace);
            CentroideDansPore(m) = CentroideDansPore(m) + nnz(classe == 1);
            CentroideSurInterface(m) = CentroideSurInterface(m) + nnz(classe == 2);
            CentroideDansSolide(m) = CentroideDansSolide(m) + nnz(classe == 0);
            CentroideHorsDomaine(m) = CentroideHorsDomaine(m) + nnz(classe == 3);

            idsGeometriques = ids(okGeometrie);
            garder = (classe == 1) | ...
                     (opts.GarderInterface & classe == 2);
            nk = nnz(garder);
            if nk > 0
                range = (cumulGardes+1):(cumulGardes+nk);
                idsGardes(range) = uint32(idsGeometriques(garder));
                classesGardees(range) = classe(garder);
                cumulGardes = cumulGardes + nk;
            end
        end

        idsGardes = idsGardes(1:cumulGardes);
        classesGardees = classesGardees(1:cumulGardes);
        tetra = tetra0(double(idsGardes),:);
        qualite = q0(double(idsGardes));
        tetra_ids_origine = idsGardes;
        classe_centroides = classesGardees;
        clear tetra0 q0 idsGardes classesGardees

        TetraConserves(m) = size(tetra,1);
        if N > 0
            TauxConservation_pct(m) = 100*TetraConserves(m)/N;
        end
        if ~isempty(qualite)
            QualiteMoyenneFinale(m) = mean(qualite);
            QualiteMedianeFinale(m) = median(qualite);
        end
        TempsNettoyage_s(m) = toc(horloge);

        % Controle comptable des categories, exclusives avant selection.
        total = Degeneres(m) + PresquePlats(m) + ...
                CentroideDansPore(m) + CentroideSurInterface(m) + ...
                CentroideDansSolide(m) + CentroideHorsDomaine(m);
        assert(total == N, 'Erreur de comptage : %s',nom);

        metadata = struct( ...
            'methode',nom, ...
            'fichierSource',fichier, ...
            'seuilQualite',opts.SeuilQualite, ...
            'seuilDetRelatif',opts.SeuilDetRelatif, ...
            'garderInterface',opts.GarderInterface, ...
            'toleranceFace',opts.ToleranceFace, ...
            'note',['Filtrage par centroide dans union des cubes du pore ; ' ...
                    'ne garantit pas la contenance integrale du tetraedre.']);
        fichierSortie = fullfile(sortie,['delaunay_nettoye_' nom '.mat']);
        save(fichierSortie,'points','tetra','qualite', ...
            'tetra_ids_origine','classe_centroides','metadata','-v7.3');

        fprintf('  Degeneres           : %d\n',Degeneres(m));
        fprintf('  Presque plats       : %d\n',PresquePlats(m));
        fprintf('  Centroides pore     : %d\n',CentroideDansPore(m));
        fprintf('  Centroides interface: %d\n',CentroideSurInterface(m));
        fprintf('  Centroides solide   : %d\n',CentroideDansSolide(m));
        fprintf('  Hors domaine        : %d\n',CentroideHorsDomaine(m));
        fprintf('  CONSERVES           : %d / %d (%.2f%%)\n', ...
            TetraConserves(m),N,TauxConservation_pct(m));
        fprintf('  Temps               : %.3f s\n',TempsNettoyage_s(m));
    end

    Methode = methodes(:);
    statsN = table(Methode,NombrePoints,TetraInitiaux,Degeneres, ...
        PresquePlats,CentroideDansPore,CentroideSurInterface, ...
        CentroideDansSolide,CentroideHorsDomaine,TetraConserves, ...
        TauxConservation_pct,QualiteMoyenneFinale, ...
        QualiteMedianeFinale,TempsNettoyage_s);

    fprintf('\n===== RECAPITULATIF DES 7 METHODES =====\n');
    disp(statsN);
    writetable(statsN,fullfile(sortie,'statistiques_nettoyage_tetra.csv'));
    save(fullfile(sortie,'statistiques_nettoyage_tetra.mat'), ...
        'statsN','opts','-v7.3');
    fprintf('\nFichiers nettoyes enregistres dans : %s\n',sortie);
end

function opts = defaut(opts,champ,valeur)
    if ~isfield(opts,champ) || isempty(opts.(champ))
        opts.(champ) = valeur;
    end
end

function classe = classer_centroides_pore(S,G,tol)
% Classification du centroide dans l'union des CUBES voxelises.
% Code : 0 solide, 1 pore strict (au sens de la grille voxelisee),
%        2 interface pore-solide, 3 en dehors du volume numerique.
% Un centre sur une face/arete/sommet est classe en examinant les
% 2, 4 ou 8 voxels qui partagent ce point.
    n = size(G,1);
    classe = zeros(n,1,'uint8');
    if n == 0, return; end
    sz = double(size(S));

    hors = any(bsxfun(@lt,G,0.5-tol) | ...
               bsxfun(@gt,G,sz+0.5+tol),2);

    frac = G-floor(G);
    surDemi = abs(frac-0.5) <= tol;
    ia = floor(G+0.5);
    ib = ia;
    ia(surDemi) = floor(G(surDemi));
    ib(surDemi) = ia(surDemi)+1;

    aPore = false(n,1);
    aSolide = false(n,1);
    for b1 = 0:1
        for b2 = 0:1
            for b3 = 0:1
                I = ia(:,1) + b1*(ib(:,1)-ia(:,1));
                J = ia(:,2) + b2*(ib(:,2)-ia(:,2));
                K = ia(:,3) + b3*(ib(:,3)-ia(:,3));
                inside = I>=1 & I<=sz(1) & J>=1 & J<=sz(2) & ...
                         K>=1 & K<=sz(3);
                isPore = false(n,1);
                if any(inside)
                    idx = sub2ind(size(S),I(inside),J(inside),K(inside));
                    isPore(inside) = S(idx);
                end
                aPore = aPore | isPore;
                aSolide = aSolide | ~isPore;
            end
        end
    end
    classe(aPore & ~aSolide) = uint8(1);
    classe(aPore & aSolide) = uint8(2);
    classe(hors) = uint8(3);
end
