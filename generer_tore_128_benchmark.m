function [imgTore, resultatsTore] = generer_tore_128_benchmark(lancerBenchmark)
% GENERER_TORE_128_BENCHMARK
% Tore analytique voxelise dans un volume 128 x 128 x 128.
% Convention STRICTEMENT identique aux benchmarks BD12 et DRP-372 :
%   1 = domaine d'interet (pore), 0 = exterieur (solide).
%
% Geometrie (centre des voxels, coordonnees MATLAB 1..128) :
%   (sqrt((x-c)^2 + (y-c)^2) - R)^2 + (z-c)^2 <= r^2
% avec c=64.5, R=36 et r=16 (unites : voxels).
%
% Appel :
%   [imgTore, resultatsTore] = generer_tore_128_benchmark(true);
% Ou uniquement generer le tore, sans calculs :
%   [imgTore, ~] = generer_tore_128_benchmark(false);
%
% Le programme de benchmark utilise sans aucune modification :
%   main_benchmark_7_frontieres_step1_sequentiel_audit(img,dossier)
%
% Dependances : Image Processing Toolbox, programme de benchmark existant.

    if nargin < 1 || isempty(lancerBenchmark)
        lancerBenchmark = true;
    end
    validateattributes(lancerBenchmark,{'logical','numeric'}, {'scalar'}, mfilename);
    lancerBenchmark = logical(lancerBenchmark);

    N = 128;
    R = 36;        % grand rayon : axe du tore
    r = 16;        % petit rayon : rayon du tube
    c = (N+1)/2;   % 64.5 : centre geom. du volume
    tailleVoxel = 1; % distances exprimees en voxels

    assert(R > r && R+r < (N-1)/2, ...
        'Dimensions du tore incompatibles avec le volume.');

    dossier = fullfile(pwd,'TORE_SYNTHETIQUE_128');
    dossierBench = fullfile(dossier,'BENCHMARK_TORE_128_NEUF');

    if lancerBenchmark && isfolder(dossierBench)
        contenu = dir(dossierBench);
        if any(~ismember({contenu.name},{'.','..'}))
            error(['Le benchmark refuse de reprendre un dossier non vide : %s\n' ...
                 'Renommer le dossier de sortie ou le deplacer pour un NOUVEL essai.'], ...
                 dossierBench);
        end
    end
    if ~isfolder(dossier), mkdir(dossier); end

    %% 1. Image binaire echantillonnee au CENTRE des voxels
    [x,y,z] = ndgrid(1:N,1:N,1:N);
    distanceAxe = sqrt((x-c).^2 + (y-c).^2);
    distanceSection = sqrt((distanceAxe-R).^2 + (z-c).^2);
    signedDistance = single(distanceSection-r); % <0 dedans, >0 dehors
    imgTore = uint8(signedDistance <= 0);
    clear x y z distanceAxe distanceSection

    %% 2. Grandeurs de reference, utiles pour l'article
    volAnalytique = 2*pi^2*R*r^2;
    surfAnalytique = 4*pi^2*R*r;
    nVoxelsPore = nnz(imgTore);
    fractionPore = nVoxelsPore/numel(imgTore);

    S = (imgTore ~= 0);
    CC6 = bwconncomp(S,6);
    nombreFaces = nnz(S(1:end-1,:,:) ~= S(2:end,:,:)) + ...
                  nnz(S(:,1:end-1,:) ~= S(:,2:end,:)) + ...
                  nnz(S(:,:,1:end-1) ~= S(:,:,2:end));
    toucheBord = any(S(1,:,:),'all') || any(S(end,:,:),'all') || ...
                 any(S(:,1,:),'all') || any(S(:,end,:),'all') || ...
                 any(S(:,:,1),'all') || any(S(:,:,end),'all');

    torus = struct('N',N,'R_vox',R,'r_vox',r, ...
                   'centre_vox',c,'tailleVoxel',tailleVoxel, ...
                   'equation', ...
                   '(sqrt((x-c)^2+(y-c)^2)-R)^2+(z-c)^2<=r^2');
    mesures = struct('voxelsPore',nVoxelsPore, ...
        'porosite_pct',100*fractionPore, ...
        'volumeAnalytique_vox3',volAnalytique, ...
        'surfaceAnalytique_vox2',surfAnalytique, ...
        'surfaceVoxelisee_faces',nombreFaces, ...
        'nbComposantes6',CC6.NumObjects, ...
        'toucheBord',toucheBord);

    assert(CC6.NumObjects == 1,'Le tore voxelise est fragmente en connectivite 6.');
    assert(~toucheBord,'Le tore touche les limites du volume.');

    % "img" est fourni en plus de "imgTore" pour faciliter le chargement.
    img = imgTore; %#ok<NASGU>
    save(fullfile(dossier,'tore_analytique_128.mat'), ...
        'img','imgTore','signedDistance','torus','mesures','-v7.3');

    fprintf('\n========== TORE SYNTHETIQUE 128^3 ==========\n');
    fprintf('Grand rayon R = %g voxels ; petit rayon r = %g voxels\n',R,r);
    fprintf('Centre = (%g, %g, %g)\n',c,c,c);
    fprintf('Domaine (1=pore) : %d voxels (%.6f %%)\n',nVoxelsPore,100*fractionPore);
    fprintf('Volume analytique : %.3f voxels^3\n',volAnalytique);
    fprintf('Surface analytique : %.3f voxels^2\n',surfAnalytique);
    fprintf('Surface par faces voxelisees : %d faces (pas une aire lisse)\n',nombreFaces);
    fprintf('Composantes connexes (6) : %d | touche le bord : %d\n', ...
        CC6.NumObjects,toucheBord);
    fprintf('Image sauvegardee : %s\n',fullfile(dossier,'tore_analytique_128.mat'));

    %% 3. Representation tridimensionnelle de l'union des voxels
    fig3D = figure('Color','w','Name','Tore synthetique 128^3');
    p = patch(isosurface(double(imgTore),0.5));
    p.FaceColor = [0.16,0.54,0.82];
    p.EdgeColor = 'none';
    axis equal; axis tight; box on; grid on;
    xlabel('X (voxels)'); ylabel('Y (voxels)'); zlabel('Z (voxels)');
    title(sprintf('Tore analytique voxelise - R=%g, r=%g, 128^3',R,r));
    view(35,25); camlight('headlight'); lighting gouraud;
    try
        exportgraphics(fig3D,fullfile(dossier,'tore_128_3D.png'),'Resolution',300);
    catch
        saveas(fig3D,fullfile(dossier,'tore_128_3D.png'));
    end

    %% 4. Coupes orthogonales (verifier que le trou est bien present)
    nCoupe = N/2;
    figCoupes = figure('Color','w','Name','Coupes du tore 128^3');
    images = {imgTore(:,:,nCoupe), squeeze(imgTore(:,nCoupe,:)), ...
              squeeze(imgTore(nCoupe,:,:))};
    titres = {'XY (z=64)','XZ (y=64)','YZ (x=64)'};
    for a = 1:3
        subplot(1,3,a);
        imagesc(images{a}); axis image; colormap(gray(2));
        title(titres{a}); xlabel('voxels'); ylabel('voxels');
    end
    try
        exportgraphics(figCoupes,fullfile(dossier,'tore_128_coupes.png'),'Resolution',250);
    catch
        saveas(figCoupes,fullfile(dossier,'tore_128_coupes.png'));
    end

    %% 5. Lancer exactement le meme programme que BD12 et DRP-372
    resultatsTore = [];
    if lancerBenchmark
        if exist('main_benchmark_7_frontieres_step1_sequentiel_audit','file') ~= 2
            error(['Programme main_benchmark_7_frontieres_step1_sequentiel_audit.m introuvable.\n' ...
                   'Placer les deux fichiers .m dans le meme dossier MATLAB,\n' ...
                   'ou ajouter le dossier du benchmark au path.']);
        end
        logFile = fullfile(dossier,'journal_benchmark_tore_128.txt');
        diary(logFile);
        nettoyageDiary = onCleanup(@() diary('off')); %#ok<NASGU>
        fprintf('\n>>> LANCEMENT DU BENCHMARK TORE : 7 frontieres, TD/TC 0/5/10\n');
        resultatsTore = main_benchmark_7_frontieres_step1_sequentiel_audit( ...
            imgTore,dossierBench);
        fprintf('\n>>> BENCHMARK TORE TERMINE. Dossier : %s\n',dossierBench);
        diary('off');
    else
        fprintf('Benchmark non lance. Pour le lancer :\n');
        fprintf('[imgTore, resultatsTore] = generer_tore_128_benchmark(true);\n');
    end
end
