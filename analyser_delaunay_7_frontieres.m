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
