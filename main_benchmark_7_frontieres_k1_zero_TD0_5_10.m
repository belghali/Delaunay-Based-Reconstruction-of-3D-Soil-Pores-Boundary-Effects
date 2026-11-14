function resultats = main_benchmark_7_frontieres_k1_zero_TD0_5_10(img,dossier)
% MAIN_BENCHMARK_7_FRONTIERES_K1_ZERO_TD0_5_10
% Pipeline COMPLET repartant de ZERO, n'importe quelle image 3D.
%
% 1) Sept frontieres I6/I18/I26/E6/E18/E26/F (extérieur=solide)
% 2) Delaunay sur la totalite de chaque frontiere
% 3) Nettoyage des tetraedres (qualite + centroide)
% 4) Spheres circonscrites, dedoublonnage, centres dans le pore
% 5) TD sur les centres des voxels, TD=0 geometrique
% 6) Selection gloutonne TC pour 84 combinaisons
% 7) Unions de toutes les spheres avant TC, 21 combinaisons
% 8) Comparaison avant/apres TC et verification E6/E18/E26
%
% AUCUNE SUBDIVISION, k=1. AUCUN RESULTAT ANTERIEUR REUTILISE.
% AUCUNE REPRISE, MEME APRES INTERRUPTION. Nouveau dossier obligatoire.
% NOTE : TD>0, TC et couverture sont des mesures par centres de voxels;
% ce ne sont pas des integrales exactes sur des cubes partiels.
%
% Appel :
% resultats=main_benchmark_7_frontieres_k1_zero_TD0_5_10(img, ...
%                  'resultats_sol_512_k1_nouveau');
%
% Image : 255 ou 1 = pore ; 0 = solide. Conversion si 0 = pore :
%   img = uint8(imgInitial==0)*255;

    if nargin<1, error('Fournir une image 3D.'); end
    validateattributes(img,{'logical','numeric'}, ...
        {'real','nonempty','3d'},mfilename,'img');
    if any(img(:)~=0 & img(:)~=1 & img(:)~=255)
        error('Image binaire attendue : 0=solide, 1 ou 255=pore.');
    end
    S=(img~=0);
    if ~any(S(:)) || all(S(:))
        error('Le volume doit contenir du pore et du solide.');
    end
    if nargin<2 || isempty(dossier)
        dims=size(S);
        dossier=sprintf('benchmark_%dx%dx%d_k1_NEUF',dims(1),dims(2),dims(3));
    end
    % Rien ne doit etre ecrase ou reutilise. Un nouveau dossier par essai.
    if isfolder(dossier)
        fichiers=dir(dossier);
        autres=~ismember({fichiers.name},{'.','..'});
        if any(autres)
            error(['Le dossier %s n''est pas vide. ' ...
                'Choisir un NOUVEAU dossier pour repartir de zero.'],dossier);
        end
    else
        mkdir(dossier);
    end
    addpath(fileparts(mfilename('fullpath')));

    optsNet=struct('SeuilQualite',1e-3, ...
        'SeuilDetRelatif',1e-12,'ToleranceFace',1e-10, ...
        'GarderInterface',true,'TailleBloc',100000);
    optsBoules=struct('TailleBloc',50000, ...
        'ToleranceDoublon',1e-7,'ToleranceInterface',1e-9, ...
        'ToleranceResiduelle',1e-8);
    optsTC=struct('SeuilsTD',[0 .05 .10], ...
        'SeuilsTC',[0 .05 .10 .20], ...
        'ToleranceGeometrique',1e-10, ...
        'VerboseTousLes',5000,'SauverSelections',true);

    fprintf('\n==========================================================\n');
    fprintf(' BENCHMARK COMPLET - K=1 - AUCUNE REPRISE\n');
    fprintf(' Volume : %s | pore : %d voxels\n',mat2str(size(S)),nnz(S));
    fprintf(' Dossier de sortie neuf : %s\n',dossier);
    fprintf(' 7 methodes / TD 0,5,10 / TC 0,5,10,20\n');
    fprintf('==========================================================\n');
    t0=tic;

    fprintf('\n===== 1/5 : EXTRACTION DES SEPT FRONTIERES =====\n');
    [P,statF]=extraire_7_frontieres(img);
    save(fullfile(dossier,'frontieres_7_methodes.mat'), ...
        'P','statF','-v7.3');
    writetable(statF,fullfile(dossier,'statistiques_frontieres.csv'));

    fprintf('\n===== 2/5 : DELAUNAY GLOBAL =====\n');
    statD=analyser_delaunay_7_frontieres(P,dossier);
    clear P

    fprintf('\n===== 3/5 : NETTOYAGE DES TETRAEDRES =====\n');
    statN=nettoyer_tetra_7_frontieres(img,dossier,optsNet);

    fprintf('\n===== 4/5 : CALCUL DES SPHERES =====\n');
    statB=calculer_boules_7_frontieres(img,dossier,optsBoules);

    fprintf('\n===== 5/5 : TD, TC ET COUVERTURE AVANT TC =====\n');
    [statTC,statTD,statAvant,statComp]= ...
        selection_TD_TC_voxels_7_frontieres(img,dossier,optsTC);

    %% Controle de l'equivalence des banques E, a 1e-7 voxel.
    nomsE={'E6','E18','E26'};
    cles=cell(1,3);
    for i=1:3
        f=fullfile(dossier,'boules_tetra_nettoyes', ...
            ['boules_' nomsE{i} '.mat']);
        B=load(f,'centres','rayons','classeCentres');
        ids=B.classeCentres(:)==1;
        cles{i}=sortrows(round([double(B.centres(ids,:)), ...
            double(B.rayons(ids))]*1e7));
    end
    equivalenceE=struct( ...
        'E6_E18',isequal(cles{1},cles{2}), ...
        'E6_E26',isequal(cles{1},cles{3}), ...
        'E18_E26',isequal(cles{2},cles{3}), ...
        'Tolerance_voxel',1e-7);
    fprintf('\nE6=E18 : %d | E6=E26 : %d | E18=E26 : %d\n', ...
        equivalenceE.E6_E18,equivalenceE.E6_E26, ...
        equivalenceE.E18_E26);
    save(fullfile(dossier,'equivalence_E.mat'),'equivalenceE');

    resultats=struct();
    resultats.Dimensions=size(S);
    resultats.Porosite_pct=100*nnz(S)/numel(S);
    resultats.Resolution='k=1 sans subdivision';
    resultats.Frontieres=statF;
    resultats.Delaunay=statD;
    resultats.Nettoyage=statN;
    resultats.Boules=statB;
    resultats.TD=statTD;
    resultats.TC=statTC;
    resultats.AvantTC=statAvant;
    resultats.Comparaison=statComp;
    resultats.EquivalenceE=equivalenceE;
    resultats.TempsTotal_s=toc(t0);
    resultats.Dossier=dossier;

    save(fullfile(dossier,'bilan_complet_k1.mat'), ...
        'resultats','optsNet','optsBoules','optsTC','-v7.3');
    fprintf('\nBENCHMARK TERMINE : %d TD, %d TC, %d avant TC\n', ...
        height(statTD),height(statTC),height(statAvant));
    fprintf('Temps total : %.2f s | Dossier : %s\n', ...
        resultats.TempsTotal_s,dossier);
end


% =========================================================================
% PARTIE : extraire_7_frontieres.m
% =========================================================================
function [P, stats] = extraire_7_frontieres(img)
%EXTRAIRE_7_FRONTIERES Sept representations de la frontiere d'un volume 3D.
%
% [P, stats] = extraire_7_frontieres(img)
%
% Entree :
%   img : image binaire 3D, avec 255 (ou 1/true) = PORE/OBJET
%         et 0 = SOLIDE/FOND.
%
% Sorties :
%   P.I6, P.I18, P.I26 : frontieres INTERIEURES (centres de voxels poreux)
%   P.E6, P.E18, P.E26 : frontieres EXTERIEURES (centres de voxels solides)
%   P.F              : centres des faces pore-solide, y compris les faces
%                      de decoupe du volume.
%   Chaque champ est une matrice N x 3 de coordonnees [i,j,k] (double),
%   dans le repere MATLAB de l'image originale (1-based).
%
% Convention de bord : en dehors du tableau, la phase est supposée solide.
% Pour l'implementer, une couche de solide virtuelle entoure temporairement
% l'image. Les points E peuvent donc avoir 0 ou dim+1 comme coordonnee;
% les centres de faces F peuvent avoir 0.5 ou dim+0.5.
%
% stats.PointsLiesAuBord : pour I, points ajoutes par la convention de bord;
% pour E, points situees dans la couche virtuelle; pour F, faces situees
% sur les six limites du volume.
%
% Requiert Image Processing Toolbox (imdilate).

    if ~(islogical(img) || isnumeric(img)) || ~isreal(img) || ...
            isempty(img) || ndims(img) ~= 3
        error('img doit etre une image binaire 3D non vide.');
    end
    if any(img(:) ~= 0 & img(:) ~= 1 & img(:) ~= 255)
        error('Les valeurs permises sont 0 et 255 (ou 0 et 1).');
    end

    S = (img ~= 0);              % true = pore, false = solide
    sz = size(S);

    % Une couche de voxels solides virtuels tout autour de l'image.
    Sp = false(sz + 2);
    Sp(2:end-1, 2:end-1, 2:end-1) = S;

    % Voisinages 6, 18, 26 (centre inclus, sans consequence sur le resultat).
    [di,dj,dk] = ndgrid(-1:1, -1:1, -1:1);
    d1 = abs(di) + abs(dj) + abs(dk);
    V = {d1 <= 1, d1 <= 2, true(3,3,3)};
    suffixes = {'6','18','26'};

    P = struct();
    bordI = zeros(3,1);
    bordE = zeros(3,1);

    for t = 1:3
        voisinage = V{t};

        % Tout point de pore adjacent a un solide (y compris hors-volume).
        MI = Sp & imdilate(~Sp, voisinage);
        % Tout point de solide adjacent a un pore (y compris couche virtuelle).
        ME = ~Sp & imdilate(Sp, voisinage);

        nomI = ['I' suffixes{t}];
        nomE = ['E' suffixes{t}];
        P.(nomI) = masque_vers_points(MI) - 1;
        P.(nomE) = masque_vers_points(ME) - 1;

        % Vrais nouveaux points interieurs dus aux limites du volume.
        MI_sans_bord = S & imdilate(~S, voisinage);
        MI_dans_img = MI(2:end-1, 2:end-1, 2:end-1);
        bordI(t) = nnz(MI_dans_img & ~MI_sans_bord);

        % Points exterieurs situes dans la couche virtuelle.
        Q = P.(nomE);
        hors = false(size(Q,1),1);
        for a = 1:3
            hors = hors | Q(:,a) < 1 | Q(:,a) > sz(a);
        end
        bordE(t) = nnz(hors);
    end

    % Points au centre de chaque face separant deux voxels de phases
    % differentes, y compris les faces entre pore et solide virtuel.
    faces = cell(3,1);
    for a = 1:3
        low = {':',':',':'};
        high = {':',':',':'};
        low{a} = 1:size(Sp,a)-1;
        high{a} = 2:size(Sp,a);
        A = Sp(low{:});
        B = Sp(high{:});
        [i,j,k] = ind2sub(size(A), find(A ~= B));
        Q = double([i,j,k]);
        Q(:,a) = Q(:,a) + 0.5;
        faces{a} = Q - 1;  % coordonnees dans l'image non rembourree
    end
    P.F = vertcat(faces{:});

    % Faces placees sur les limites physiques du domaine numerique.
    Q = P.F;
    sur_limite = false(size(Q,1),1);
    for a = 1:3
        sur_limite = sur_limite | Q(:,a) == 0.5 | Q(:,a) == sz(a)+0.5;
    end
    bordF = nnz(sur_limite);

    methodes = {'I6';'I18';'I26';'E6';'E18';'E26';'F'};
    NombrePoints = zeros(7,1);
    for t = 1:7
        NombrePoints(t) = size(P.(methodes{t}),1);
    end
    PointsLiesAuBord = [bordI; bordE; bordF];
    stats = table(methodes, NombrePoints, PointsLiesAuBord, ...
        'VariableNames', {'Methode','NombrePoints','PointsLiesAuBord'});

    fprintf('\nDimensions du volume : %s\n', mat2str(sz));
    fprintf('Voxels poreux        : %d\n', nnz(S));
    fprintf('Convention de bord   : exterieur = solide\n');
    disp(stats);
end

function Q = masque_vers_points(M)
    [i,j,k] = ind2sub(size(M), find(M));
    Q = double([i,j,k]);
end


% =========================================================================
% PARTIE : analyser_delaunay_7_frontieres.m
% =========================================================================
function statsD = analyser_delaunay_7_frontieres(P, dossier)
% ANALYSER_DELAUNAY_7_FRONTIERES
% Triangule de facon identique les 7 jeux de points de frontiere et
% compare la qualite des tetraedres. Ne calcule PAS encore les boules.
%
% ENTREE
%   P       structure creee par extraire_7_frontieres(img), avec les
%           champs I6, I18, I26, E6, E18, E26 et F.
%   dossier dossier de sortie (defaut : 'resultats_delaunay_64')
%
% SORTIE
%   statsD  tableau des statistiques. Un fichier .mat par methode est
%           ecrit avec points, tetraedres, qualite et parametres.
%
% CONVENTIONS
%   - Coordinates [i,j,k] MATLAB, espace isotrope (voxels).
%   - Points virtuels exterieurs (0 ou dim+1) conserves tels quels.
%   - Les 7 methodes utilisent delaunayTriangulation sans changement
%     de points, de perturbation aleatoire ou de filtrage initial.
%   - La triangulation remplit l'enveloppe convexe des points; elle
%     n'est PAS une triangulation du seul espace poral.
%   - La qualite tetraedrique utilise la metrique mean-ratio :
%       q = 12*(3*V)^(2/3) / somme_des_6_aretes_au_carre
%     q=1 pour un tetraedre regulier, q approche 0 s'il est plat.
%
% EXEMPLE
%   [P,~] = extraire_7_frontieres(img);
%   statsD = analyser_delaunay_7_frontieres(P);

    if nargin < 2 || isempty(dossier)
        dossier = 'resultats_delaunay_64';
    end
    if ~exist(dossier,'dir')
        mkdir(dossier);
    end

    methodes = {'I6','I18','I26','E6','E18','E26','F'};
    seuilQualite = 1e-3;
    tailleBloc = 200000;
    n = numel(methodes);

    NombrePoints = zeros(n,1);
    NombreTetraedres = zeros(n,1);
    QualiteMin = nan(n,1);
    QualiteMediane = nan(n,1);
    QualiteMoyenne = nan(n,1);
    TetraPresquePlats = zeros(n,1);
    PourcentagePresquePlats = nan(n,1);
    TempsDelaunay_s = zeros(n,1);
    TempsQualite_s = zeros(n,1);

    fprintf('\n===== ETAPE 2 : DELAUNAY DES 7 FRONTIERES =====\n');
    fprintf('Seuil de qualite q < %.4g\n',seuilQualite);

    for m = 1:n
        nom = methodes{m};
        points = double(P.(nom));

        if size(points,2) ~= 3 || size(points,1) < 4 || ...
                any(~isfinite(points(:)))
            error('Ensemble %s : il faut au moins 4 points 3D finis.',nom);
        end

        if size(unique(points,'rows'),1) ~= size(points,1)
            error('Ensemble %s : doublons trouves ; corriger avant Delaunay.',nom);
        end

        NombrePoints(m) = size(points,1);
        fprintf('\n[%d/%d] %s : %d points\n',m,n,nom,NombrePoints(m));

        ticD = tic;
        DT = delaunayTriangulation(points);
        TempsDelaunay_s(m) = toc(ticD);

        % Utiliser les points tels que references par la triangulation.
        points = DT.Points;
        tetra = DT.ConnectivityList;
        NombreTetraedres(m) = size(tetra,1);
        clear DT

        ticQ = tic;
        qualite = zeros(NombreTetraedres(m),1);

        for debut = 1:tailleBloc:NombreTetraedres(m)
            finBloc = min(NombreTetraedres(m),debut+tailleBloc-1);
            ids = debut:finBloc;
            idsTet = tetra(ids,:);

            A = points(idsTet(:,1),:);
            B = points(idsTet(:,2),:);
            C = points(idsTet(:,3),:);
            D = points(idsTet(:,4),:);

            AB = B-A; AC = C-A; AD = D-A;
            BC = C-B; BD = D-B; CD = D-C;

            triple = sum(AB .* cross(AC,AD,2),2);
            vol = abs(triple)/6;

            sommeL2 = sum(AB.^2 + AC.^2 + AD.^2 + ...
                         BC.^2 + BD.^2 + CD.^2,2);

            q = 12 * (3*vol).^(2/3) ./ max(sommeL2,realmin);
            qualite(ids) = min(1,max(0,q));
        end
        TempsQualite_s(m) = toc(ticQ);

        if ~isempty(qualite)
            QualiteMin(m) = min(qualite);
            QualiteMediane(m) = median(qualite);
            QualiteMoyenne(m) = mean(qualite);
            TetraPresquePlats(m) = nnz(qualite < seuilQualite);
            PourcentagePresquePlats(m) = ...
                100*TetraPresquePlats(m)/NombreTetraedres(m);
        end

        fprintf('  Tetraedres          : %d\n',NombreTetraedres(m));
        fprintf('  Qualite moyenne     : %.5f\n',QualiteMoyenne(m));
        fprintf('  Qualite mediane     : %.5f\n',QualiteMediane(m));
        fprintf('  Quasi-plats q<%.4g : %d (%.4f%%)\n', ...
                seuilQualite,TetraPresquePlats(m),PourcentagePresquePlats(m));
        fprintf('  Temps Delaunay      : %.3f s\n',TempsDelaunay_s(m));
        fprintf('  Temps qualite       : %.3f s\n',TempsQualite_s(m));

        % Sauvegarde sans perte (double pour les points et qualites).
        metadata = struct('methode',nom, ...
                          'seuilQualite',seuilQualite, ...
                          'versionMATLAB',version, ...
                          'note', ['Delaunay global de la frontiere; ' ...
                                   'tetraedres non filtres pore/solide']);
        fichier = fullfile(dossier,['delaunay_' nom '.mat']);
        save(fichier,'points','tetra','qualite','metadata','-v7.3');
    end

    Methode = methodes(:);
    statsD = table(Methode,NombrePoints,NombreTetraedres, ...
        QualiteMin,QualiteMediane,QualiteMoyenne, ...
        TetraPresquePlats,PourcentagePresquePlats, ...
        TempsDelaunay_s,TempsQualite_s);

    disp(statsD);
    writetable(statsD,fullfile(dossier,'comparaison_delaunay.csv'));
    save(fullfile(dossier,'comparaison_delaunay.mat'),'statsD');
    fprintf('\nResultats sauvegardes dans : %s\n',dossier);
end


% =========================================================================
% PARTIE : nettoyer_tetra_7_frontieres.m
% =========================================================================
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


% =========================================================================
% PARTIE : calculer_boules_7_frontieres.m
% =========================================================================
function statsB = calculer_boules_7_frontieres(img, dossier, opts)
% CALCULER_BOULES_7_FRONTIERES
% Calcule les spheres circonscrites des tetraedres DEJA NETTOYES pour
% les sept representations I6, I18, I26, E6, E18, E26 et F.
%
% UTILISATION :
%   statsB = calculer_boules_7_frontieres(img,'resultats_delaunay_64');
%
% CONVENTIONS :
%   img = 255 (ou 1/true) pour le pore, 0 pour le solide.
%   Les coordonnees [i,j,k] sont celles des centres des voxels MATLAB.
%   Un voxel de centre entier (i,j,k) occupe un cube de cote 1.
%   L'exterieur du volume est considere comme solide virtuel.
%
% FICHIERS LUS :
%   dossier/nettoyes/delaunay_nettoye_<methode>.mat
%   Variables : points, tetra, tetra_ids_origine, [qualite].
%
% FICHIERS CREES :
%   dossier/boules_tetra_nettoyes/boules_<methode>.mat
%   dossier/boules_tetra_nettoyes/statistiques_boules.csv
%   dossier/boules_tetra_nettoyes/statistiques_boules.mat
%
% VARIABLES SAUVEGARDEES POUR CHAQUE METHODE :
%   centres        [K x 3] double : centres des spheres distinctes
%   rayons         [K x 1] double : rayons des spheres distinctes
%   idTetNettoye   [K x 1] uint32: ligne dans le fichier nettoye
%   idTetOriginal  [K x 1] uint32: ligne dans le Delaunay original
%   classeCentres  [K x 1] uint8 : 0 solide, 1 pore strict,
%                                  2 interface, 3 hors domaine
%   metadata       structure des parametres et conventions.
%
% IMPORTANT : PAS de filtrage TD, PAS de selection TC et PAS de controle
% de contenance integrale des spheres dans le pore.
% 'unique' regroupe les quadruplets [centre,rayon] quantifies a une
% resolution fixee ; ce n'est pas un clustering par distance euclidienne.

    if nargin < 2 || isempty(dossier)
        dossier = 'resultats_delaunay_64';
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = option_defaut(opts,'TailleBloc',50000);
    opts = option_defaut(opts,'ToleranceDoublon',1e-7);
    opts = option_defaut(opts,'ToleranceInterface',1e-9);
    opts = option_defaut(opts,'ToleranceResiduelle',1e-8);

    validateattributes(img,{'logical','numeric'}, ...
        {'3d','real','nonempty'},mfilename,'img');
    if any(img(:) ~= 0 & img(:) ~= 1 & img(:) ~= 255)
        error('Image binaire attendue : 255/1 = pore, 0 = solide.');
    end
    validateattributes(opts.TailleBloc,{'numeric'}, ...
        {'scalar','integer','positive'});
    validateattributes(opts.ToleranceDoublon,{'numeric'}, ...
        {'scalar','real','finite','positive'});
    validateattributes(opts.ToleranceInterface,{'numeric'}, ...
        {'scalar','real','finite','>=',0,'<',1e-4});
    validateattributes(opts.ToleranceResiduelle,{'numeric'}, ...
        {'scalar','real','finite','>=',0});

    S = (img ~= 0);
    methodes = {'I6','I18','I26','E6','E18','E26','F'};
    nm = numel(methodes);
    dossierEntree = fullfile(dossier,'nettoyes');
    dossierSortie = fullfile(dossier,'boules_tetra_nettoyes');
    if ~exist(dossierSortie,'dir'), mkdir(dossierSortie); end

    N_Tetra = zeros(nm,1);
    N_SpheresFinies = zeros(nm,1);
    N_SpheresInvalides = zeros(nm,1);
    N_SpheresDistinctes = zeros(nm,1);
    N_Doublons = zeros(nm,1);
    Doublons_pct = zeros(nm,1);
    CentresPore = zeros(nm,1);
    CentresInterface = zeros(nm,1);
    CentresSolide = zeros(nm,1);
    CentresHorsDomaine = zeros(nm,1);
    CentresPore_pct = nan(nm,1);
    RayonMin = nan(nm,1);
    RayonMedian = nan(nm,1);
    RayonMoyen = nan(nm,1);
    RayonMax = nan(nm,1);
    ErreurCircumMax = nan(nm,1);
    SpheresErreurResiduelle = zeros(nm,1);
    TempsCalcul_s = zeros(nm,1);
    TempsTotal_s = zeros(nm,1);

    fprintf('\n===== ETAPE 4 : SPHERES DES TETRAEDRES NETTOYES =====\n');
    fprintf('Pore : %d voxels, dimensions %s\n',nnz(S),mat2str(size(S)));
    fprintf('Tolerance quantification doublons : %.2g voxel\n', ...
        opts.ToleranceDoublon);

    for m = 1:nm
        nom = methodes{m};
        fichier = fullfile(dossierEntree, ...
            ['delaunay_nettoye_' nom '.mat']);
        if ~isfile(fichier)
            error('Fichier nettoye introuvable : %s',fichier);
        end
        tTotal = tic;
        data = load(fichier,'points','tetra','tetra_ids_origine');
        points = double(data.points);
        tetra = double(data.tetra);
        N = size(tetra,1);
        N_Tetra(m) = N;
        if size(points,2) ~= 3 || size(tetra,2) ~= 4
            error('%s : format de points/tetra invalide.',nom);
        end
        if isfield(data,'tetra_ids_origine') && ...
                numel(data.tetra_ids_origine) == N
            idsOrigine = uint32(data.tetra_ids_origine(:));
        else
            error('%s : tetra_ids_origine manquant ou invalide.',nom);
        end
        clear data
        fprintf('\n[%d/7] %s : %d tetraedres nettoyes\n',m,nom,N);

        C = zeros(N,3);
        R = zeros(N,1);
        tCalc = tic;
        if N > 0
            % triangulation() ne recalcule pas Delaunay.
            TR = triangulation(tetra,points);
            for premier = 1:opts.TailleBloc:N
                ids = (premier:min(N,premier+opts.TailleBloc-1)).';
                [Cb,Rb] = circumcenter(TR,ids);
                C(ids,:) = Cb;
                R(ids) = Rb;
            end
            clear TR
        end
        TempsCalcul_s(m) = toc(tCalc);

        valide = all(isfinite(C),2) & isfinite(R) & R > 0;
        N_SpheresFinies(m) = nnz(valide);
        N_SpheresInvalides(m) = N - N_SpheresFinies(m);
        idsValides = find(valide);
        C = C(valide,:);
        R = R(valide);

        if isempty(R)
            centres = zeros(0,3);
            rayons = zeros(0,1);
            idTetNettoye = zeros(0,1,'uint32');
            idTetOriginal = zeros(0,1,'uint32');
            classeCentres = zeros(0,1,'uint8');
        else
            % Dedoublonnage geometrique sur les 4 valeurs : i, j, k, r.
            valeurMax = max(abs([C,R]),[],'all');
            if valeurMax/opts.ToleranceDoublon >= flintmax
                error('%s : precision insuffisante pour les cles de doublons.',nom);
            end
            cles = round([C,R] / opts.ToleranceDoublon);
            [~,ia] = unique(cles,'rows','stable');
            clear cles
            centres = C(ia,:);
            rayons = R(ia);
            idTetNettoye = uint32(idsValides(ia));
            idTetOriginal = idsOrigine(double(idTetNettoye));
            classeCentres = classer_centres_cubes(S,centres, ...
                opts.ToleranceInterface);

            RayonMin(m) = min(rayons);
            RayonMedian(m) = median(rayons);
            RayonMoyen(m) = mean(rayons);
            RayonMax(m) = max(rayons);

            % Controle numerique : chaque sphere doit passer par les 4
            % sommets de son tetraedre de provenance.
            Ti = tetra(double(idTetNettoye),:);
            resid = zeros(numel(rayons),1);
            for v = 1:4
                V = points(Ti(:,v),:);
                d = sqrt(sum((V-centres).^2,2));
                resid = max(resid,abs(d-rayons));
            end
            ErreurCircumMax(m) = max(resid);
            SpheresErreurResiduelle(m) = nnz( ...
                resid > opts.ToleranceResiduelle * max(1,rayons));
            clear Ti resid
        end
        clear C R idsValides

        N_SpheresDistinctes(m) = numel(rayons);
        N_Doublons(m) = N_SpheresFinies(m)-N_SpheresDistinctes(m);
        if N_SpheresFinies(m)>0
            Doublons_pct(m) = 100*N_Doublons(m)/N_SpheresFinies(m);
        end
        CentresPore(m) = nnz(classeCentres == 1);
        CentresInterface(m) = nnz(classeCentres == 2);
        CentresSolide(m) = nnz(classeCentres == 0);
        CentresHorsDomaine(m) = nnz(classeCentres == 3);
        if N_SpheresDistinctes(m)>0
            CentresPore_pct(m) = 100*CentresPore(m)/N_SpheresDistinctes(m);
        end
        assert(CentresPore(m)+CentresInterface(m)+ ...
            CentresSolide(m)+CentresHorsDomaine(m) == ...
            N_SpheresDistinctes(m),'Erreur de classification des centres.');

        metadata = struct('methode',nom, ...
            'source',fichier, ...
            'unite','voxel', ...
            'ToleranceDoublon',opts.ToleranceDoublon, ...
            'ToleranceInterface',opts.ToleranceInterface, ...
            'ToleranceResiduelle',opts.ToleranceResiduelle, ...
            'convention',['Centre de voxel (i,j,k), cube de cote 1; ' ...
                          'exterieur du volume = solide virtuel.'], ...
            'note',['Spheres circonscrites aux tetraedres filtres ' ...
                    'par centroide. Ni TD ni TC appliques.']);
        save(fullfile(dossierSortie,['boules_' nom '.mat']), ...
            'centres','rayons','idTetNettoye', ...
            'idTetOriginal','classeCentres','metadata','-v7.3');

        TempsTotal_s(m) = toc(tTotal);
        fprintf('  Spheres finies / invalides : %d / %d\n', ...
            N_SpheresFinies(m),N_SpheresInvalides(m));
        fprintf('  Distinctes / doublons      : %d / %d (%.2f%%)\n', ...
            N_SpheresDistinctes(m),N_Doublons(m),Doublons_pct(m));
        fprintf('  Centres pore/interface     : %d / %d\n', ...
            CentresPore(m),CentresInterface(m));
        fprintf('  Centres solide/hors        : %d / %d\n', ...
            CentresSolide(m),CentresHorsDomaine(m));
        fprintf('  Rayons min/med/moy/max     : %.4f / %.4f / %.4f / %.4f\n', ...
            RayonMin(m),RayonMedian(m),RayonMoyen(m),RayonMax(m));
        fprintf('  Erreur circonscription max : %.3g\n',ErreurCircumMax(m));
        fprintf('  Temps calcul / total       : %.3f / %.3f s\n', ...
            TempsCalcul_s(m),TempsTotal_s(m));
    end

    Methode = methodes(:);
    statsB = table(Methode,N_Tetra,N_SpheresFinies, ...
        N_SpheresInvalides,N_SpheresDistinctes,N_Doublons, ...
        Doublons_pct,CentresPore,CentresInterface, ...
        CentresSolide,CentresHorsDomaine,CentresPore_pct, ...
        RayonMin,RayonMedian,RayonMoyen,RayonMax, ...
        ErreurCircumMax,SpheresErreurResiduelle, ...
        TempsCalcul_s,TempsTotal_s);
    disp(statsB);
    writetable(statsB,fullfile(dossierSortie,'statistiques_boules.csv'));
    save(fullfile(dossierSortie,'statistiques_boules.mat'),'statsB','opts');
    fprintf('\nResultats enregistres dans : %s\n',dossierSortie);
end

function opts = option_defaut(opts,nom,valeur)
    if ~isfield(opts,nom) || isempty(opts.(nom))
        opts.(nom) = valeur;
    end
end

function classe = classer_centres_cubes(S,C,tol)
% Classe les centres par rapport a l'union des cubes voxelises.
% 0 = solide, 1 = pore strict, 2 = interface, 3 = hors domaine.
% En cas de coordonnee demi-entiere, considere les deux voxels voisins.
% L'exterieur du domaine est du solide virtuel.
    n = size(C,1);
    classe = zeros(n,1,'uint8');
    if n == 0, return; end
    sz = double(size(S));
    hors = any(bsxfun(@lt,C,0.5-tol) | ...
               bsxfun(@gt,C,sz+0.5+tol),2);
    fraction = C-floor(C);
    surDemi = abs(fraction-0.5) <= tol;
    bas = floor(C+0.5);
    haut = bas;
    bas(surDemi) = floor(C(surDemi));
    haut(surDemi) = bas(surDemi)+1;

    dansPore = false(n,1);
    dansSolide = false(n,1);
    for a = 0:1
        for b = 0:1
            for c = 0:1
                ii = bas(:,1)+a*(haut(:,1)-bas(:,1));
                jj = bas(:,2)+b*(haut(:,2)-bas(:,2));
                kk = bas(:,3)+c*(haut(:,3)-bas(:,3));
                valide = ii>=1 & ii<=sz(1) & jj>=1 & jj<=sz(2) & ...
                         kk>=1 & kk<=sz(3);
                valeurs = false(n,1);
                if any(valide)
                    ind = sub2ind(size(S),ii(valide),jj(valide),kk(valide));
                    valeurs(valide) = S(ind);
                end
                dansPore = dansPore | valeurs;
                dansSolide = dansSolide | ~valeurs;
            end
        end
    end
    classe(dansPore & ~dansSolide) = uint8(1);
    classe(dansPore & dansSolide) = uint8(2);
    classe(hors) = uint8(3);
end


% =========================================================================
% PARTIE : selection_TD_TC_voxels_7_frontieres.m
% =========================================================================
function [statsTC,statsTD,statsAvantTC,statsComparaison] = selection_TD_TC_voxels_7_frontieres(img,dossier,opts)
% SELECTION_TD_TC_VOXELS_7_FRONTIERES
% AUCUNE SUBDIVISION : k = 1, un centre de voxel par voxel.
% 7 methodes x 3 TD x 4 TC = 84 configurations et 21 unions avant TC.
%
% TD=0 : test geometrique deterministe de penetration dans les cubes solides.
% TD>0 : part des CENTRES de voxels non poreux dans chaque sphere.
% TC : part des centres de voxels deja couverts par l'union precedente.
%    L'ordre et le critere de la selection gloutonne sont conserves.
%
% IMPORTANT : k=1 n'est PAS une integration exacte des volumes continus.
% Les voxels hors image sont traites comme du solide virtuel.
% Les 4 seuils TC sont evalues simultanement pour chaque seuil TD.
% Aucun fichier preexistant n'est lu pour TD/TC; banques de spheres du
% lancement courant seulement. Pas de checkpoint, pas de reprise.
%
% FICHIERS d'entree :
%   dossier/boules_tetra_nettoyes/boules_I6.mat, etc.
%
% FICHIERS de sortie :
%   dossier/TD_TC_voxels_k1/statistiques_TD.csv
%   dossier/TD_TC_voxels_k1/statistiques_TC.csv
%   dossier/TD_TC_voxels_k1/statistiques_avant_TC.csv
%   dossier/TD_TC_voxels_k1/comparaison_avant_apres_TC.csv
%   dossier/TD_TC_voxels_k1/selection_...mat

    if nargin<3 || isempty(opts), opts=struct(); end
    opts = defaut_selection(opts,'SeuilsTD',[0 .05 .10]);
    opts = defaut_selection(opts,'SeuilsTC',[0 .05 .10 .20]);
    opts = defaut_selection(opts,'VerboseTousLes',5000);
    opts = defaut_selection(opts,'ToleranceGeometrique',1e-10);
    opts = defaut_selection(opts,'SauverSelections',true);

    S = (img~=0);
    sz = double(size(S));
    methodes = {'I6','I18','I26','E6','E18','E26','F'};
    seuilsTD = unique(double(opts.SeuilsTD(:).'),'sorted');
    seuilsTC = unique(double(opts.SeuilsTC(:).'),'sorted');
    if ~isequal(seuilsTD,[0 .05 .10]) || ~isequal(seuilsTC,[0 .05 .10 .20])
        error('Ce benchmark attend TD : [0 0.05 0.10] et TC : [0 0.05 0.10 0.20].');
    end
    rep = fullfile(dossier,'TD_TC_voxels_k1');
    if ~isfolder(rep), mkdir(rep); end

    % Structures de sortie (lignes converties en tables a la fin).
    TDr = cell(21,9);
    TCr = cell(84,13);
    Ar = cell(21,10);
    Cr = cell(84,8);
    nTDrow=0; nTCrow=0; nArow=0; nCrow=0;
    nPore = nnz(S);

    fprintf('\n=====================================================\n');
    fprintf(' TD / TC SANS SUBDIVISION (k=1) - DEPART DE ZERO\n');
    fprintf(' Image : %s, pore : %d voxels\n',mat2str(sz),nPore);
    fprintf(' 7 methodes, 3 TD, 4 TC => 84 configurations\n');
    fprintf('=====================================================\n');

    for m=1:7
        nom=methodes{m};
        tMethode=tic;
        fichier=fullfile(dossier,'boules_tetra_nettoyes',['boules_' nom '.mat']);
        if ~isfile(fichier), error('Fichier introuvable : %s',fichier); end
        D=load(fichier,'centres','rayons','classeCentres');
        idsOrigine=find(D.classeCentres(:)==1);
        C=double(D.centres(idsOrigine,:));
        R=double(D.rayons(idsOrigine));
        N=numel(R);
        if isempty(R)
            warning('%s : aucune sphere a centre pore.',nom);
        end
        if size(C,1)~=N || size(C,2)~=3 || any(~isfinite(C(:))) || ...
                any(~isfinite(R)) || any(R<=0)
            error('%s : banque de spheres invalide.',nom);
        end

        % On entoure la phase pore d'une marge solide. Ainsi, chaque centre
        % de voxel exterieur couvert par une sphere entre dans TD et TC.
        % La marge est propre a la methode pour ne pas allouer trop de RAM.
        if isempty(R), marge=2; else, marge=ceil(max(R))+2; end
        taillePad=sz+2*marge;
        nbPoints=prod(taillePad);
        if ~isfinite(nbPoints) || nbPoints > flintmax
            error('%s : grille depassant la precision d''indexation MATLAB.',nom);
        end
        Sp=false(taillePad);
        Sp((1:sz(1))+marge,(1:sz(2))+marge,(1:sz(3))+marge)=S;
        fprintf('\n[%d/7] %s : %d boules pore | grille etendue %s (%.2f Gio/masque)\n', ...
            m,nom,N,mat2str(taillePad),nbPoints/2^30);

        %% Audit TD : une fois par sphere, sur TOUS les centres entiers
        TD=nan(N,1);
        NombrePointsSphere=zeros(N,1,'uint32');
        EstInscrite=false(N,1);
        tTD=tic;
        for i=1:N
            ind=points_de_voxels_dans_sphere(C(i,:),R(i),taillePad,marge);
            nb=numel(ind);
            if nb==0
                % Sphere trop petite pour contenir un centre de voxel.
                continue;
            end
            NombrePointsSphere(i)=uint32(nb);
            nbSolide=nb-nnz(Sp(ind));
            TD(i)=nbSolide/nb;

            % Optimisation SANS CHANGER le test TD=0 : si une sphere
            % contient deja un centre de voxel solide, elle ne peut pas
            % etre geometriquement inscrite dans le domaine poreux.
            if nbSolide==0
                EstInscrite(i)=sphere_inscrite_dans_cubes(S,sz,C(i,:), ...
                    R(i),opts.ToleranceGeometrique);
            end
            if opts.VerboseTousLes>0 && mod(i,opts.VerboseTousLes)==0
                fprintf('  TD %s : %d/%d | %.1f s\n',nom,i,N,toc(tTD));
            end
        end
        tempsTD=toc(tTD);
        fprintf('  TD termine : %.1f s | inscrites: %d | sans point: %d\n', ...
            tempsTD,nnz(EstInscrite),nnz(NombrePointsSphere==0));

        save(fullfile(rep,['TD_banque_' nom '.mat']), ...
            'C','R','idsOrigine','TD','EstInscrite', ...
            'NombrePointsSphere','marge','-v7.3');

        %% Boucle des 3 seuils TD; 4 seuils TC EN PARALLELE LOGIQUE
        for it=1:numel(seuilsTD)
            td=seuilsTD(it);
            if td==0
                admissible=EstInscrite & NombrePointsSphere>0;
            else
                admissible=isfinite(TD) & (TD<=td+1e-12);
            end
            ids=find(admissible);
            nAdm=numel(ids);
            if nAdm>0
                % Identique au tri glouton historique : rayon decroissant,
                % coordonnees croissantes pour departager les egalites.
                [~,ordre]=sortrows([-R(ids),C(ids,:)],[1 2 3 4]);
                ids=ids(ordre);
            end

            nTDrow=nTDrow+1;
            TDr(nTDrow,:)={nom,100*td,N,nAdm,N-nAdm, ...
                mean(TD,'omitnan')*100,nnz(EstInscrite), ...
                nnz(NombrePointsSphere==0),tempsTD};
            fprintf('  TD=%g%% : %d boules admissibles\n',100*td,nAdm);

            % 4 unites de couverture TC + union integrale sans TC.
            % Les masques TC vivent SEULEMENT pendant ce seuil TD.
            U=cell(1,4);
            for j=1:4, U{j}=false(taillePad); end
            Uavant=false(taillePad);
            sel=cell(1,4); tcAcceptes=cell(1,4);
            nSel=zeros(1,4);
            TP=zeros(1,4); FP=zeros(1,4);
            for j=1:4
                sel{j}=zeros(nAdm,1,'uint32');
                tcAcceptes{j}=zeros(nAdm,1);
            end
            TPavant=0; FPavant=0;
            tGroupe=tic;

            for ii=1:nAdm
                b=ids(ii);
                lin=points_de_voxels_dans_sphere(C(b,:),R(b),taillePad,marge);
                if isempty(lin), continue; end
                nSphere=numel(lin);

                % Union avant TC : toutes les spheres admissibles.
                nouveaux=lin(~Uavant(lin));
                if ~isempty(nouveaux)
                    nNvPore=nnz(Sp(nouveaux));
                    TPavant=TPavant+nNvPore;
                    FPavant=FPavant+numel(nouveaux)-nNvPore;
                    Uavant(nouveaux)=true;
                end

                % Quatre selections gloutonnes simultanees :
                % on ne recalcule pas les points de la sphere 4 fois.
                for j=1:4
                    deja=U{j}(lin);
                    recouvrement=nnz(deja)/nSphere;
                    if recouvrement<=seuilsTC(j)+1e-12
                        nSel(j)=nSel(j)+1;
                        sel{j}(nSel(j))=uint32(b);
                        tcAcceptes{j}(nSel(j))=recouvrement;
                        nouveaux=lin(~deja);
                        if ~isempty(nouveaux)
                            gainPore=nnz(Sp(nouveaux));
                            TP(j)=TP(j)+gainPore;
                            FP(j)=FP(j)+numel(nouveaux)-gainPore;
                            U{j}(nouveaux)=true;
                        end
                    end
                end

                if opts.VerboseTousLes>0 && mod(ii,opts.VerboseTousLes)==0
                    fprintf('    %s TD=%g%% : %d/%d | TC0=%d TC5=%d TC10=%d TC20=%d | %.1f s\n', ...
                        nom,100*td,ii,nAdm,nSel,toc(tGroupe));
                end
            end
            dureeGroupe=toc(tGroupe);

            cAvant=100*TPavant/max(1,nPore);
            dAvant=100*FPavant/max(1,TPavant+FPavant);
            iAvant=100*TPavant/max(1,nPore+FPavant);
            nArow=nArow+1;
            Ar(nArow,:)={nom,100*td,nAdm,cAvant,dAvant,iAvant, ...
                TPavant,FPavant,dureeGroupe,marge};

            fprintf('    AVANT TC : couverture %.2f%%, debordement %.2f%%\n', ...
                cAvant,dAvant);

            for j=1:4
                tc=seuilsTC(j);
                cApres=100*TP(j)/max(1,nPore);
                dApres=100*FP(j)/max(1,TP(j)+FP(j));
                ioU=100*TP(j)/max(1,nPore+FP(j));
                precision=100*TP(j)/max(1,TP(j)+FP(j));
                nTCrow=nTCrow+1;
                TCr(nTCrow,:)={nom,100*td,100*tc,N,nAdm,nSel(j), ...
                    cApres,dApres,precision,ioU,TP(j),FP(j),dureeGroupe};
                nCrow=nCrow+1;
                if cAvant>0, conserve=100*cApres/cAvant; else, conserve=NaN; end
                if nAdm>0, reduction=100*(1-nSel(j)/nAdm); else, reduction=NaN; end
                Cr(nCrow,:)={nom,100*td,100*tc,cAvant,cApres, ...
                    cAvant-cApres,conserve,reduction};
                fprintf('      TC=%2g%% : %6d boules | couv=%.2f%% | TDunion=%.2f%% | IoU=%.2f%%\n', ...
                    100*tc,nSel(j),cApres,dApres,ioU);

                if opts.SauverSelections
                    localIds=double(sel{j}(1:nSel(j)));
                    centresSelectionnes=C(localIds,:);
                    rayonsSelectionnes=R(localIds);
                    indicesBanqueBrute=idsOrigine(localIds);
                    TDSelection=TD(localIds);
                    TCSelection=tcAcceptes{j}(1:nSel(j));
                    meta=struct('Methode',nom,'TD',td,'TC',tc, ...
                        'k',1,'TypeTDPositif','centres de voxels', ...
                        'TypeTDZero','geometrie exacte des cubes', ...
                        'TypeTC','centres de voxels, glouton rayon decroissant');
                    sortie=fullfile(rep,sprintf('selection_%s_TD%02d_TC%02d.mat', ...
                        nom,round(100*td),round(100*tc)));
                    save(sortie,'centresSelectionnes','rayonsSelectionnes', ...
                        'indicesBanqueBrute','TDSelection','TCSelection', ...
                        'meta','-v7.3');
                end
            end
            clear U Uavant sel tcAcceptes
        end
        clear Sp C R D
        fprintf('  %s termine : %.1f s\n',nom,toc(tMethode));
    end

    %% Construire les tables / sauvegarde
    colonnesTD={'SeuilTD_pct','CentresPore','AdmissiblesTD', ...
        'RejeteesTD','TDMoyenVoxel_pct','N_InscritesExactes', ...
        'N_SansCentreVoxel','TempsAuditTD_s'};
    colonnesTC={'SeuilTD_pct','SeuilTC_pct','CentresPore', ...
        'AdmissiblesTD','BoulesSelectionnees','CouverturePore_pct', ...
        'DebordementUnion_pct','Precision_pct','IoU_pct', ...
        'VoxelsPoreCouverts','VoxelsHorsPoreCouverts','TempsGroupeTC_s'};
    colonnesAvant={'SeuilTD_pct','BoulesAdmissibles','CouverturePore_pct', ...
        'DebordementUnion_pct','IoU_pct','VoxelsPoreCouverts', ...
        'VoxelsHorsPoreCouverts','TempsGroupeTC_s','Marge'};
    colonnesComp={'SeuilTD_pct','SeuilTC_pct','CouvertureAvantTC_pct', ...
        'CouvertureApresTC_pct','CouverturePerdue_points', ...
        'CouvertureConservee_pct','ReductionBoules_pct'};

    statsTD=[table(TDr(:,1),'VariableNames',{'Methode'}), ...
        array2table(cell2mat(TDr(:,2:end)),'VariableNames',colonnesTD)];
    statsTC=[table(TCr(:,1),'VariableNames',{'Methode'}), ...
        array2table(cell2mat(TCr(:,2:end)),'VariableNames',colonnesTC)];
    statsAvantTC=[table(Ar(:,1),'VariableNames',{'Methode'}), ...
        array2table(cell2mat(Ar(:,2:end)),'VariableNames',colonnesAvant)];
    statsComparaison=[table(Cr(:,1),'VariableNames',{'Methode'}), ...
        array2table(cell2mat(Cr(:,2:end)),'VariableNames',colonnesComp)];

    assert(height(statsTD)==21 && height(statsTC)==84 && ...
           height(statsAvantTC)==21 && height(statsComparaison)==84, ...
           'Nombre inattendu de configurations.');

    writetable(statsTD,fullfile(rep,'statistiques_TD.csv'));
    writetable(statsTC,fullfile(rep,'statistiques_TC.csv'));
    writetable(statsAvantTC,fullfile(rep,'statistiques_avant_TC.csv'));
    writetable(statsComparaison,fullfile(rep,'comparaison_avant_apres_TC.csv'));
    save(fullfile(rep,'statistiques_completes.mat'), ...
        'statsTD','statsTC','statsAvantTC','statsComparaison','opts','-v7.3');

    fprintf('\n===== RESULTATS ENREGISTRES : %s =====\n',rep);
end

function lin=points_de_voxels_dans_sphere(c,r,taillePad,marge)
% Enumerer TOUS les centres entiers dans la sphere, y compris hors image.
% Coordonnees du voxel etendu p -> p-marge dans l'image originale.
    tol=1e-12*max(1,r*r);
    lo=max([1 1 1],ceil(c-r+marge-1e-12*max(1,r)));
    hi=min(taillePad,floor(c+r+marge+1e-12*max(1,r)));
    if any(lo>hi), lin=zeros(0,1); return; end
    x=(lo(1):hi(1))'-marge-c(1);
    y=(lo(2):hi(2))'-marge-c(2);
    z=(lo(3):hi(3))'-marge-c(3);
    d2=reshape(x.^2,[],1,1)+reshape(y.^2,1,[],1)+ ...
        reshape(z.^2,1,1,[]);
    ids=find(d2<=r*r+tol);
    if isempty(ids), lin=zeros(0,1); return; end
    [ii,jj,kk]=ind2sub(size(d2),ids);
    lin=sub2ind(taillePad,lo(1)-1+ii,lo(2)-1+jj,lo(3)-1+kk);
end

function ok=sphere_inscrite_dans_cubes(S,sz,c,r,tol)
% Test deterministe de penetration dans les cubes des voxels solides.
% Les contacts tangents sont autorises a la tolerance numerique donnee.
    marge=tol*max(1,r);
    if any(c-r<0.5-marge) || any(c+r>sz+0.5+marge)
        ok=false; return;
    end
    lo=max([1 1 1],ceil(c-r-0.5));
    hi=min(sz,floor(c+r+0.5));
    if any(lo>hi), ok=true; return; end
    [I,J,K]=ndgrid(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3));
    ids=sub2ind(sz,I(:),J(:),K(:));
    solide=~S(ids);
    if ~any(solide), ok=true; return; end
    Q=[double(I(solide)),double(J(solide)),double(K(solide))];
    d=max(abs(bsxfun(@minus,Q,c))-0.5,0);
    distance2=sum(d.^2,2);
    ok=~any(distance2<r*r-tol*max(1,r*r));
end

function opts=defaut_selection(opts,champ,valeur)
    if ~isfield(opts,champ) || isempty(opts.(champ))
        opts.(champ)=valeur;
    end
end