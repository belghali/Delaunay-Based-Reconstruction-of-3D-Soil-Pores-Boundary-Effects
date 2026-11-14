function resultats = main_benchmark_7_frontieres_64(img, dossier, etapes)
% MAIN_BENCHMARK_7_FRONTIERES_64
% Orchestre le protocole FINAL du benchmark 64^3 :
%   1. Extraction des 7 frontieres I6 I18 I26 E6 E18 E26 F
%   2. Delaunay et mesure de qualite des tetraedres
%   3. Nettoyage des tetraedres (qualite, degenerescence, centroide)
%   4. Spheres circonscrites, dedoublonnage et classement des centres
%   5. TD VOLUMIQUE (grille exhaustive k=4) + TC glouton (112 essais)
%      Pour TD=0, test geometrique deterministe d'inscription.
%   6. Union de TOUTES les boules admissibles avant TC (28 essais)
%      et comparaison aux 112 essais apres TC.
%   7. Controle de l'equivalence geometrique E6/E18/E26.
%
% Entree : img binaire 3D 64x64x64, 255 ou 1 = PORE, 0 = SOLIDE.
% ATTENTION : si votre image d'origine a 0 = PORE, inversez le masque
% AVANT d'appeler cette fonction : img = uint8(imgOriginal==0)*255;
%
% Exemples :
%   resultats = main_benchmark_7_frontieres_64(img);
%   resultats = main_benchmark_7_frontieres_64(img, 'resultats_delaunay_64');
%
% Pour REPRENDRE apres le calcul des boules sans les recalculer :
%   e = struct('Extraction',false,'Delaunay',false, ...
%       'Nettoyage',false,'Boules',false);
%   resultats = main_benchmark_7_frontieres_64( ...
%       img, 'resultats_delaunay_64', e);
%
% Ce MAIN reprend les fonctions existantes sans modifier la selection
% gloutonne. Il est destine a 64^3 ; la version 512^3 exige un moteur
% volumique par blocs (k=4 conserve, mais plus de grille globale RAM).

    if nargin < 2 || isempty(dossier)
        dossier = 'resultats_delaunay_64';
    end
    if nargin < 3 || isempty(etapes)
        etapes = struct();
    end

    % Valeur par defaut : toutes les etapes actives.
    noms = {'Extraction','Delaunay','Nettoyage','Boules', ...
            'TD_TC','AvantTC','VerifierE'};
    for i = 1:numel(noms)
        if ~isfield(etapes,noms{i})
            etapes.(noms{i}) = true;
        end
        validateattributes(etapes.(noms{i}),{'logical','numeric'}, ...
            {'scalar'},mfilename,['etapes.' noms{i}]);
        etapes.(noms{i}) = logical(etapes.(noms{i}));
    end

    % Les fonctions associees doivent etre dans le meme dossier que ce main.
    dossierScripts = fileparts(mfilename('fullpath'));
    addpath(dossierScripts);

    validateattributes(img,{'numeric','logical'}, ...
        {'nonempty','real','3d'},mfilename,'img');
    if ~isequal(size(img),[64 64 64])
        error(['Ce MAIN reproduit le benchmark 64x64x64. ' ...
               'Ne pas l''utiliser directement sur 512^3 : ' ...
               'la grille k=4 du moteur TC actuel serait trop volumineuse.']);
    end
    if any(img(:) ~= 0 & img(:) ~= 1 & img(:) ~= 255)
        error('Convention attendue : pore=255 ou 1, solide=0.');
    end
    if nnz(img)==0 || nnz(img)==numel(img)
        error('Le volume doit contenir a la fois du pore et du solide.');
    end

    if ~exist(dossier,'dir')
        mkdir(dossier);
    end

    dependances = {};
    if etapes.Extraction
        dependances{end+1} = 'extraire_7_frontieres';
    end
    if etapes.Delaunay
        dependances{end+1} = 'analyser_delaunay_7_frontieres';
    end
    if etapes.Nettoyage
        dependances{end+1} = 'nettoyer_tetra_7_frontieres';
    end
    if etapes.Boules
        dependances{end+1} = 'calculer_boules_7_frontieres';
    end
    if etapes.TD_TC
        dependances{end+1} = 'selection_TC_volumique_7_frontieres_exhaustif';
    end
    if etapes.AvantTC
        dependances{end+1} = 'couverture_union_avant_TC_7_frontieres';
    end
    for i = 1:numel(dependances)
        if isempty(which(dependances{i}))
            error('Fonction MATLAB introuvable : %s.m',dependances{i});
        end
    end

    % Parametres figes du protocole 64^3.
    optsNet = struct('SeuilQualite',1e-3, ...
                     'SeuilDetRelatif',1e-12, ...
                     'ToleranceFace',1e-10, ...
                     'GarderInterface',true, ...
                     'TailleBloc',100000);

    optsBoules = struct('TailleBloc',50000, ...
                        'ToleranceDoublon',1e-7, ...
                        'ToleranceInterface',1e-9, ...
                        'ToleranceResiduelle',1e-8);

    optsTC = struct('SeuilsTD',[0 0.05 0.10 0.20], ...
                    'SeuilsTC',[0 0.05 0.10 0.20], ...
                    'SousVoxelsParAxe',4, ...
                    'ToleranceGeometrique',1e-10, ...
                    'VerboseTousLes',5000, ...
                    'SauverSelections',true, ...
                    'MaxPointsGrille',8e7);

    optsAvant = struct('SeuilsTD',optsTC.SeuilsTD, ...
                       'SousVoxelsParAxe',4, ...
                       'VerboseTousLes',5000, ...
                       'SauverMasques',false, ...
                       'MaxPointsGrille',8e7);

    % Metadonnees, sans stocker tout le volume dans les statistiques.
    resultats = struct();
    resultats.Dossier = dossier;
    resultats.Dimensions = size(img);
    resultats.VoxelsPoreux = nnz(img ~= 0);
    resultats.Porosite_pct = 100 * nnz(img ~= 0)/numel(img);
    resultats.Etapes = etapes;
    resultats.OptionsNettoyage = optsNet;
    resultats.OptionsBoules = optsBoules;
    resultats.OptionsTDTC = optsTC;
    resultats.OptionsAvantTC = optsAvant;

    fprintf('\n====================================================\n');
    fprintf('       BENCHMARK COMPLET : 7 FRONTIERES / 64^3\n');
    fprintf('====================================================\n');
    fprintf('Pore : %d voxels (%.3f %%)\n', ...
        resultats.VoxelsPoreux,resultats.Porosite_pct);
    fprintf('Dossier : %s\n',dossier);
    fprintf('TD volumique seulement ; k=4 ; TC glouton inchange.\n');
    tMain = tic;

    %% ETAPE 1 : SEPT FRONTIERES
    if etapes.Extraction
        fprintf('\n===== 1/7 EXTRACTION DES FRONTIERES =====\n');
        [P,statsFrontieres] = extraire_7_frontieres(img);
        save(fullfile(dossier,'frontieres_7_methodes.mat'), ...
            'P','statsFrontieres','-v7.3');
        writetable(statsFrontieres, ...
            fullfile(dossier,'statistiques_frontieres.csv'));
        resultats.Frontieres = statsFrontieres;
    end

    %% ETAPE 2 : TRIANGULATIONS DE DELAUNAY
    if etapes.Delaunay
        fprintf('\n===== 2/7 DELAUNAY =====\n');
        if ~exist('P','var')
            fP = fullfile(dossier,'frontieres_7_methodes.mat');
            if ~isfile(fP)
                error(['Points de frontiere manquants. ' ...
                       'Relancer l''extraction ou fournir %s.'],fP);
            end
            D = load(fP,'P');
            P = D.P;
        end
        resultats.Delaunay = analyser_delaunay_7_frontieres(P,dossier);
        clear P;
    end

    %% ETAPE 3 : NETTOYAGE DES TETRAEDRES
    if etapes.Nettoyage
        fprintf('\n===== 3/7 NETTOYAGE DES TETRAEDRES =====\n');
        resultats.Nettoyage = nettoyer_tetra_7_frontieres( ...
            img,dossier,optsNet);
    end

    %% ETAPE 4 : SPHERES CIRCONSCRITES DES TETRAEDRES NETTOYES
    if etapes.Boules
        fprintf('\n===== 4/7 CALCUL DES SPHERES =====\n');
        resultats.Boules = calculer_boules_7_frontieres( ...
            img,dossier,optsBoules);
    end

    %% ETAPE 5 : TD VOLUMIQUE ET SELECTION GLOUTONNE TC
    if etapes.TD_TC
        fprintf('\n===== 5/7 TD VOLUMIQUE ET TC GLOUTON =====\n');
        [statsTC,statsTD] = ...
            selection_TC_volumique_7_frontieres_exhaustif( ...
                img,dossier,optsTC);
        resultats.TD = statsTD;                 % 28 lignes
        resultats.TC = statsTC;                 % 112 lignes
    end

    %% ETAPE 6 : UNION AVANT TC ET COMPARAISON AVANT/APRES
    if etapes.AvantTC
        fprintf('\n===== 6/7 COUVERTURE AVANT TC =====\n');
        [statsAvantTC,statsComparaison] = ...
            couverture_union_avant_TC_7_frontieres( ...
                img,dossier,optsAvant);
        resultats.CouvertureAvantTC = statsAvantTC;  % 28 lignes
        resultats.Comparaison = statsComparaison;     % 112 lignes
    end

    %% ETAPE 7 : VERIFICATION E6 / E18 / E26
    if etapes.VerifierE
        fprintf('\n===== 7/7 EQUIVALENCE DES BANQUES E =====\n');
        methodesE = {'E6','E18','E26'};
        K = cell(3,1);
        for i = 1:3
            f = fullfile(dossier,'boules_tetra_nettoyes', ...
                ['boules_' methodesE{i} '.mat']);
            if ~isfile(f)
                error('Banque de boules introuvable : %s',f);
            end
            B = load(f,'centres','rayons','classeCentres');
            id = (B.classeCentres(:) == 1);
            C = double(B.centres(id,:));
            R = double(B.rayons(id));
            K{i} = sortrows(round([C R]*1e7));
            fprintf('%s : %d boules avec centre pore\n', ...
                methodesE{i},size(K{i},1));
        end
        resultatE = struct( ...
            'E6_E18',isequal(K{1},K{2}), ...
            'E6_E26',isequal(K{1},K{3}), ...
            'E18_E26',isequal(K{2},K{3}), ...
            'ToleranceQuantification_voxel',1e-7);
        resultats.EquivalenceE = resultatE;
        fprintf('E6 = E18 : %d\n',resultatE.E6_E18);
        fprintf('E6 = E26 : %d\n',resultatE.E6_E26);
        fprintf('E18 = E26 : %d\n',resultatE.E18_E26);
        save(fullfile(dossier,'verification_banques_E.mat'), ...
            'resultatE');
        clear K;
    end

    %% BILAN FINAL
    resultats.TempsTotal_s = toc(tMain);
    save(fullfile(dossier,'bilan_main_benchmark.mat'), ...
        'resultats','-v7.3');

    fprintf('\n====================================================\n');
    fprintf('FIN DU MAIN : %.2f secondes\n',resultats.TempsTotal_s);
    fprintf('Resultats : %s\n',dossier);
    fprintf('====================================================\n');
end
