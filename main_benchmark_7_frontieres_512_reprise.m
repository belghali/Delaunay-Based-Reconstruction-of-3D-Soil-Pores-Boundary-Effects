function resultats = main_benchmark_7_frontieres_512_reprise(img,dossier,opts)
% MAIN_BENCHMARK_7_FRONTIERES_512_REPRISE
% Reprise automatique des fichiers deja crees par le MAIN 512 precedent.
% Reprend les 7 methodes, TD=0/5/10/20, TC=0/5/10/20, k=4.
% La selection gloutonne, les tests et les conventions sont inchanges.
%
% REPRISE : saute les fichiers de methode complets pour les etapes 1-4 ;
% TD et TC reutilisent aussi les sorties de l'ancien pipeline, puis ecrivent
% de nouveaux checkpoints periodiques pour les calculs inacheves.
% Les vieux programmes ne possedaient pas de checkpoints internes : on
% NE PEUT PAS reconstruire leurs iterations non sauvegardees.
%
% img : 512x512x512 ; PORE=255 ou 1 ; SOLIDE=0.
% Exemple : resultats=main_benchmark_7_frontieres_512_reprise(img, ...
%     'resultats_delaunay_512');

if nargin<2 || isempty(dossier),dossier='resultats_delaunay_512';end
if nargin<3,opts=struct();end
opts=def(opts,'Reprendre',true);
opts=def(opts,'Extraction',true);
opts=def(opts,'Delaunay',true);
opts=def(opts,'Nettoyage',true);
opts=def(opts,'Boules',true);
opts=def(opts,'TD_TC',true);
opts=def(opts,'AvantTC',true);
opts=def(opts,'VerifierE',true);
opts=def(opts,'CheckpointChaque',10000);
opts=def(opts,'VerboseTousLes',10000);
opts=def(opts,'Methodes',{'I6','I18','I26','E6','E18','E26','F'});

addpath(fileparts(mfilename('fullpath')));
validateattributes(img,{'numeric','logical'},{'3d','real','nonempty'});
if ~isequal(size(img),[512 512 512])
    error('Image 512x512x512 attendue.');
end
if any(img(:)~=0 & img(:)~=1 & img(:)~=255)
    error('Convention : PORE=255/1, SOLIDE=0.');
end
if ~isfolder(dossier),mkdir(dossier);end
S=(img~=0);
if nnz(S)==0 || nnz(S)==numel(S),error('Image uniforme.');end
fprintf('\n====== MAIN 512 REPRISE AUTOMATIQUE ======\n');
fprintf('Dossier : %s\nPore : %d voxels\n',dossier,nnz(S));
fprintf('7 frontieres ; TD=0/5/10/20 ; TC=0/5/10/20 ; k=4.\n');
fprintf('ATTENTION : les donnees preexistantes doivent correspondre a la MEME image.\n');
fprintf('Une methode en cours dans l''ancien code sans checkpoint sera relancee.\n');
resultats=struct('Dossier',dossier,'Taille',size(img),'Pore',nnz(S));
methodes=opts.Methodes;
base=fullfile(dossier,'frontieres_7_methodes.mat');
t= tic;

%% 1. Extraction, une seule fois, si manquante
if opts.Extraction
    if opts.Reprendre && complet(base,{'P','statsFrontieres'})
        fprintf('[SKIP] Frontieres existantes.\n');
    else
        fprintf('\n[1] Extraction des sept frontieres...\n');
        [P,statsFrontieres]=extraire_7_frontieres(img);
        save(base,'P','statsFrontieres','-v7.3');
        writetable(statsFrontieres,fullfile(dossier,'statistiques_frontieres.csv'));
        clear P
    end
end

%% 2. Delaunay methode par methode
if opts.Delaunay
    missing={};
    for m=1:numel(methodes)
        nom=methodes{m};
        ff=fullfile(dossier,['delaunay_' nom '.mat']);
        if opts.Reprendre && complet(ff,{'points','tetra','qualite','metadata'})
            fprintf('[SKIP] Delaunay %s deja sauvegarde.\n',nom);
        else
            missing{end+1}=nom; %#ok<AGROW>
        end
    end
    if ~isempty(missing)
        if ~complet(base,{'P'})
            error('Fichier %s manquant. Activer Extraction.',base);
        end
        L=load(base,'P'); P=L.P; clear L
        for m=1:numel(missing)
            nom=missing{m};
            fprintf('\n[2] Delaunay %s...\n',nom);
            analyser_delaunay_une_methode(P,dossier,nom);
        end
        clear P
    end
end

%% 3. Nettoyage methode par methode
if opts.Nettoyage
    o=struct('SeuilQualite',1e-3,'SeuilDetRelatif',1e-12, ...
        'ToleranceFace',1e-10,'GarderInterface',true,'TailleBloc',100000);
    for m=1:numel(methodes)
        nom=methodes{m};
        ff=fullfile(dossier,'nettoyes',['delaunay_nettoye_' nom '.mat']);
        if opts.Reprendre && complet(ff,{'points','tetra','qualite', ...
                'tetra_ids_origine','classe_centroides','metadata'})
            fprintf('[SKIP] Nettoyage %s deja sauvegarde.\n',nom);
        else
            fprintf('\n[3] Nettoyage %s...\n',nom);
            nettoyer_tetra_une_methode(img,dossier,o,nom);
        end
    end
end

%% 4. Boules methode par methode
if opts.Boules
    o=struct('TailleBloc',50000,'ToleranceDoublon',1e-7, ...
        'ToleranceInterface',1e-9,'ToleranceResiduelle',1e-8);
    for m=1:numel(methodes)
        nom=methodes{m};
        ff=fullfile(dossier,'boules_tetra_nettoyes',['boules_' nom '.mat']);
        if opts.Reprendre && complet(ff,{'centres','rayons','classeCentres','metadata'})
            fprintf('[SKIP] Boules %s deja sauvegardees.\n',nom);
        else
            fprintf('\n[4] Calcul des boules %s...\n',nom);
            calculer_boules_une_methode(img,dossier,o,nom);
        end
    end
end

%% 5. TD volumique et quatre TC simultanes, avec checkpoints
if opts.TD_TC
    fprintf('\n[5] TD volumique + TC glouton, reprise...\n');
    q=struct('SeuilsTD',[0 .05 .10 .20], ...
        'SeuilsTC',[0 .05 .10 .20], ...
        'SousVoxelsParAxe',4,'ToleranceGeometrique',1e-10, ...
        'CheckpointChaque',opts.CheckpointChaque, ...
        'VerboseTousLes',opts.VerboseTousLes, ...
        'Reprendre',opts.Reprendre,'Methodes',{methodes});
    [statsTC,statsTD]=selection_TC_512_reprise(img,dossier,q);
    resultats.TC=statsTC;resultats.TD=statsTD;
end

%% 6. Union de toutes les boules avant TC
if opts.AvantTC
    fprintf('\n[6] Union avant TC, reprise...\n');
    q=struct('SousVoxelsParAxe',4, ...
        'SeuilsTD',[0 .05 .10 .20], ...
        'VerboseTousLes',opts.VerboseTousLes, ...
        'Reprendre',opts.Reprendre,'Methodes',{methodes});
    [a,c]=couverture_avant_TC_512_reprise(img,dossier,q);
    resultats.AvantTC=a;resultats.Comparaison=c;
end

%% 7. Verifier l'identite numerique des banques E
if opts.VerifierE
    fdir=fullfile(dossier,'boules_tetra_nettoyes');
    Ms={'E6','E18','E26'};
    cles=cell(1,3);
    for i=1:3
        B=load(fullfile(fdir,['boules_' Ms{i} '.mat']), ...
            'centres','rayons','classeCentres');
        id=B.classeCentres(:)==1;
        cles{i}=sortrows(round([double(B.centres(id,:)), ...
            double(B.rayons(id))]*1e7));
    end
    resultats.EquivalenceE=[isequal(cles{1},cles{2}), ...
        isequal(cles{1},cles{3}),isequal(cles{2},cles{3})];
    fprintf('\nE6==E18 : %d ; E6==E26 : %d ; E18==E26 : %d\n', ...
        resultats.EquivalenceE);
end
resultats.TempsCetteExecution_s=toc(t);
save(fullfile(dossier,'bilan_main_reprise_512.mat'),'resultats','opts','-v7.3');
fprintf('\nMAIN TERMINE : %.1f s (temps de cette execution).\n', ...
    resultats.TempsCetteExecution_s);
end

function yes=complet(fname,vars)
yes=false;
if ~isfile(fname),return;end
try
    info=whos('-file',fname);
    present={info.name};
    yes=all(ismember(vars,present));
catch
    yes=false;
end
end
function o=def(o,n,v)
if ~isfield(o,n)||isempty(o.(n)),o.(n)=v;end
end
