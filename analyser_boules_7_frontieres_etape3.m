function statsB = analyser_boules_7_frontieres_etape3(img, dossier)
% ANALYSER_BOULES_7_FRONTIERES_ETAPE3
% Etape 3 du benchmark : spheres circonscrites des tetraedres Delaunay.
% Ne realise NI filtrage TD, NI selection TC, NI evaluation de couverture.
%
% ENTREE
%   img     : image 3D binaire ; 255 (ou true / 1) = pore, 0 = solide.
%   dossier : dossier produit par analyser_delaunay_7_frontieres
%             (defaut : 'resultats_delaunay_64').
%
% SORTIE
%   statsB : statistiques par representation de frontiere.
%   Fichiers : dossier/boules_brutes_<methode>.mat
%              dossier/comparaison_boules_brutes.csv
%
% CONVENTIONS
%   * Coordonnees MATLAB des centres de voxels : (i,j,k) = entiers.
%   * Chaque voxel est un cube de cote 1, de centre (i,j,k).
%   * L'exterieur du tableau est considere solide (domaine borne).
%   * Les centres sur une interface pore-solide sont classes a part.
%   * Dedoublonnage geometrique par quantification a 1e-7 voxel.
%   * Centres et rayons sauvegardes en precision double, sans arrondi.
%
% EXEMPLE
%   statsB = analyser_boules_7_frontieres_etape3(img, ...
%                              'resultats_delaunay_64');

    if nargin < 2 || isempty(dossier)
        dossier = 'resultats_delaunay_64';
    end
    if ~isnumeric(img) && ~islogical(img)
        error('img doit etre une image numerique ou logique.');
    end
    if ndims(img) ~= 3 || ~isreal(img) || isempty(img)
        error('img doit etre un volume 3D reel non vide.');
    end
    if any(img(:) ~= 0 & img(:) ~= 1 & img(:) ~= 255)
        error('img doit contenir uniquement 0 et 255 (ou 0 et 1).');
    end
    S = (img ~= 0);  % true = pore
    sz = size(S);

    noms = {'I6','I18','I26','E6','E18','E26','F'};
    n = numel(noms);
    tailleBloc = 50000;
    tolDoublon = 1e-7;  % tolerance de dedoublonnage en voxels
    tolInterface = 1e-8;  % tolerance numerique en voxels

    NombreTetraedres = zeros(n,1);
    SpheresFinies = zeros(n,1);
    SpheresDistinctes = zeros(n,1);
    CentresDansPore = zeros(n,1);
    CentresSurInterface = zeros(n,1);
    CentresDansSolide = zeros(n,1);
    CentresHorsDomaine = zeros(n,1);
    RayonMedian = nan(n,1);
    RayonMax = nan(n,1);
    Temps_s = zeros(n,1);

    fprintf('\n===== ETAPE 3 : SPHERES CIRCONSCRITES =====\n');
    fprintf('Masque de pore : %d voxels\n',nnz(S));
    fprintf('Tol. dedoublonnage : %.1g voxel\n',tolDoublon);

    for m = 1:n
        nom = noms{m};
        f = fullfile(dossier,['delaunay_' nom '.mat']);
        if ~isfile(f)
            error('Fichier introuvable : %s',f);
        end
        data = load(f,'points','tetra');
        points = double(data.points);
        tetra = double(data.tetra);
        N = size(tetra,1);
        NombreTetraedres(m) = N;

        fprintf('\n[%d/%d] %s : %d tetraedres\n',m,n,nom,N);
        chrono = tic;
        % La fonction circumcenter fonctionne aussi sur un objet triangulation
        % construit a partir des tetraedres deja sauvegardes ; nous ne
        % recalculons pas Delaunay, pour conserver exactement les memes tetra.
        TR = triangulation(tetra,points);
        C = nan(N,3);
        R = nan(N,1);

        for debut = 1:tailleBloc:N
            ids = (debut:min(N,debut+tailleBloc-1)).';
            [Cb,Rb] = circumcenter(TR,ids);
            C(ids,:) = Cb;
            R(ids) = Rb;
        end
        clear TR data

        valide = all(isfinite(C),2) & isfinite(R) & R > 0;
        SpheresFinies(m) = nnz(valide);
        centresValides = C(valide,:);
        rayonsValides = R(valide);
        idsValides = find(valide);
        clear C R

        if ~isempty(rayonsValides)
            % Trouver les spheres geometriquement equivalentes, sans
            % degrader la precision des spheres qui seront sauvegardees.
            cle = round([centresValides,rayonsValides] / tolDoublon);
            [~,ia] = unique(cle,'rows','stable');
            clear cle
            centres = centresValides(ia,:);
            rayons = rayonsValides(ia);
            tetraID = idsValides(ia);
            clear centresValides rayonsValides idsValides

            % classe : 1 = pore, 2 = interface, 0 = solide, 3 = hors domaine
            classe = classer_centres_voxels(S,centres,tolInterface);
            SpheresDistinctes(m) = numel(rayons);
            CentresDansPore(m) = nnz(classe == 1);
            CentresSurInterface(m) = nnz(classe == 2);
            CentresDansSolide(m) = nnz(classe == 0);
            CentresHorsDomaine(m) = nnz(classe == 3);
            RayonMedian(m) = median(rayons);
            RayonMax(m) = max(rayons);
        else
            centres = zeros(0,3);
            rayons = zeros(0,1);
            tetraID = zeros(0,1);
            classe = zeros(0,1,'uint8');
        end

        Temps_s(m) = toc(chrono);
        metadata = struct('Methode',nom, ...
            'ToleranceDedoublonnage',tolDoublon, ...
            'ToleranceInterface',tolInterface, ...
            'Note',['Sphere circonscrite, sans verification de ' ...
                    'contenance de la sphere dans le pore.']);

        save(fullfile(dossier,['boules_brutes_' nom '.mat']), ...
            'centres','rayons','tetraID','classe','metadata','-v7.3');

        fprintf('  Spheres finies          : %d\n',SpheresFinies(m));
        fprintf('  Spheres distinctes      : %d\n',SpheresDistinctes(m));
        fprintf('  Centres dans pore       : %d\n',CentresDansPore(m));
        fprintf('  Centres sur interface   : %d\n',CentresSurInterface(m));
        fprintf('  Centres dans solide     : %d\n',CentresDansSolide(m));
        fprintf('  Centres hors domaine    : %d\n',CentresHorsDomaine(m));
        fprintf('  Rayon median / max      : %.4f / %.4f voxels\n', ...
            RayonMedian(m),RayonMax(m));
        fprintf('  Temps                  : %.3f s\n',Temps_s(m));
    end

    Methode = noms(:);
    statsB = table(Methode,NombreTetraedres,SpheresFinies, ...
        SpheresDistinctes,CentresDansPore,CentresSurInterface, ...
        CentresDansSolide,CentresHorsDomaine,RayonMedian,RayonMax,Temps_s);

    disp(statsB);
    writetable(statsB,fullfile(dossier,'comparaison_boules_brutes.csv'));
    save(fullfile(dossier,'comparaison_boules_brutes.mat'),'statsB');
    fprintf('\nResultats enregistres dans : %s\n',dossier);
end

function classe = classer_centres_voxels(S,C,tol)
% Classe les points relativement a l'UNION DES CUBES poreux.
% En presence d'un centre sur une face/arete/sommet voxelique,
% il examine tous les voxels adjacents (jusqu'a 8).
% Les voxels virtuels au-dela du volume sont consideres solides.

    n = size(C,1);
    sz = size(S);
    classe = zeros(n,1,'uint8');
    if n == 0, return; end

    borneMin = repmat(0.5,1,3);
    borneMax = double(sz) + 0.5;
    hors = any(bsxfun(@lt,C,borneMin-tol) | ...
               bsxfun(@gt,C,borneMax+tol),2);

    partieFractionnaire = C - floor(C);
    surPlan = abs(partieFractionnaire - 0.5) <= tol;

    % Lorsque la coordonnee coincide avec une demi-valeur, il y a
    % potentiellement deux voxels adjacents selon cette direction.
    haut = floor(C + 0.5);
    haut(surPlan) = floor(C(surPlan)) + 1;
    bas = haut - double(surPlan);

    dansPore = false(n,1);
    dansSolide = false(n,1);
    for a = 0:1
        for b = 0:1
            for c = 0:1
                ii = bas(:,1) + a*(haut(:,1)-bas(:,1));
                jj = bas(:,2) + b*(haut(:,2)-bas(:,2));
                kk = bas(:,3) + c*(haut(:,3)-bas(:,3));

                dans = ii>=1 & ii<=sz(1) & ...
                       jj>=1 & jj<=sz(2) & ...
                       kk>=1 & kk<=sz(3);
                valeurs = false(n,1);
                if any(dans)
                    idx = sub2ind(sz,ii(dans),jj(dans),kk(dans));
                    valeurs(dans) = S(idx);
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
