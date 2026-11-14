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
