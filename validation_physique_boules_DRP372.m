function T = validation_physique_boules_DRP372(dossierBenchmark, fichierImage, fichierLBM, fichierLBMcsv, TDs, TCs, sauverMasques)
% VALIDATION_PHYSIQUE_BOULES_DRP372
% Post-traitement des spheres selectionnees par
% main_benchmark_7_frontieres_step1_sequentiel_audit.m.
%
% Aucun nouveau Delaunay, aucun nouvel LBM.
% Chaque union de boules est voxelisee exactement selon la convention
% du benchmark : un centre de voxel entier (x,y,z) est dans une boule si
% ||(x,y,z)-c||^2 <= r^2 + 1e-12*max(1,r^2).
%
% Les mesures de porosite, surface specifique et connectivite utilisent
% le volume 256^3 INTERIEUR uniquement. L'IoU de reference du benchmark,
% qui inclut aussi le debordement hors image, est lu depuis son CSV.
%
% MODE D'EMPLOI (a executer depuis le dossier des donnees) :
%   T = validation_physique_boules_DRP372;
%
% Ou explicitement, pour les sept configurations TD10/TC10 :
%   T = validation_physique_boules_DRP372( ...
%      'BENCHMARK_DRP372_256_NEUF', ...
%      '374_02_024_256.mat', 'LBM.mat', 'LBM.csv', 10, 10, false);
%
% Pour toutes les 63 configurations :
%   T = validation_physique_boules_DRP372( ...
%      'BENCHMARK_DRP372_256_NEUF', ...
%      '374_02_024_256.mat', 'LBM.mat', 'LBM.csv', ...
%      [0 5 10], [0 5 10], false);
%
% SauverMasques=true produit un fichier .mat logique par reconstruction.
%
% Limitation importante : la permeabilite KC est un PROXY de permeabilite,
% calibre sur le LBM de l'original ; ce n'est PAS la permeabilite physique
% des reconstructions. L'indicateur de flux conserve n'est PAS une DNS
% des reconstructions.

    if nargin < 1 || isempty(dossierBenchmark)
        dossierBenchmark = 'BENCHMARK_DRP372_256_NEUF';
    end
    if nargin < 2 || isempty(fichierImage)
        fichierImage = '374_02_024_256.mat';
    end
    if nargin < 3 || isempty(fichierLBM)
        fichierLBM = 'LBM.mat';
    end
    if nargin < 4 || isempty(fichierLBMcsv)
        fichierLBMcsv = 'LBM.csv';
    end
    if nargin < 5 || isempty(TDs)
        TDs = 10;
    end
    if nargin < 6 || isempty(TCs)
        TCs = 10;
    end
    if nargin < 7 || isempty(sauverMasques)
        sauverMasques = false;
    end
    assert(all(ismember(TDs,[0 5 10])) && all(ismember(TCs,[0 5 10])), ...
        'Seuils autorises : 0, 5, 10.');
    assert(isfile(fichierImage), 'Image manquante : %s', fichierImage);
    assert(isfile(fichierLBM), 'Champ LBM manquant : %s', fichierLBM);
    assert(isfile(fichierLBMcsv), 'CSV LBM manquant : %s', fichierLBMcsv);

    repertoireSelection = fullfile(dossierBenchmark,'TD_TC_step1_sequentiel_audit_k1');
    assert(isfolder(repertoireSelection), ...
        'Dossier selections introuvable : %s', repertoireSelection);

    % 0=pore et 1=solide dans le fichier DRP-372 d'origine.
    D = load(fichierImage,'bin');
    assert(isfield(D,'bin'),'Variable bin absente de %s.', fichierImage);
    assert(all(D.bin(:)==0 | D.bin(:)==1), ...
        'Le volume d''origine n''est pas binaire 0/1.');
    P = (D.bin == 0);  % Convention du benchmark : true=pore.
    clear D
    sz = size(P); sz = sz(1:3);
    nPore = nnz(P);
    nVox = numel(P);
    phi0 = nPore / nVox;
    Sv0 = surface_interne_voxel(P) / nVox;  % en (voxel)^(-1)
    fprintf('\nIMAGE REFERENCE : [%d %d %d], porosite = %.6f %%\n', ...
        sz(1),sz(2),sz(3),100*phi0);
    fprintf('Surface pore-solide : %.0f faces internes ; Sv=%.8g /voxel\n', ...
        surface_interne_voxel(P),Sv0);

    % Le champ 'uz' est aligne directement avec 'bin' quand les deux
    % fichiers sont lus avec LOAD dans MATLAB (pas avec un lecteur HDF5 brut).
    D = load(fichierLBM,'uz');
    assert(isfield(D,'uz'),'Variable uz absente de %s.', fichierLBM);
    assert(isequal(size(D.uz),sz), 'Dimensions de bin et uz incompatibles.');
    uz = double(D.uz);
    clear D
    amplitudeMax = max(abs(uz(:)));
    if amplitudeMax <= 0
        error('Le champ uz de reference est nul partout.');
    end
    nVitessesHorsPore = nnz(abs(uz(~P)) > 1e-8*amplitudeMax);
    assert(nVitessesHorsPore == 0, ...
        ['Champ uz mal aligne avec le pore (vitesses dans le solide). ' ...
         'Verifier l''orientation, ne pas permuter arbitrairement les axes.']);
    vitesseAbs = abs(uz);
    vitesseAbs(~P) = 0;
    sommeVitesse = sum(vitesseAbs(:));
    clear uz
    assert(sommeVitesse>0,'Somme des vitesses de reference nulle.');

    Kref = permeabilite_depuis_csv(fichierLBMcsv);
    fprintf('K_LBM reference = %.12g [lus : unite du fichier LBM.csv]\n',Kref);
    fprintf('Somme |uz| sur le pore original = %.9g\n',sommeVitesse);
    fprintf('Note : la permeabilite reconstruite est un estimateur KC calibre.\n\n');

    % Reference : les composantes poreuses originales.
    [nC0,perco0,frac0,grand0] = connectivite_6(P);
    ref = struct('Dimensions',sz,'NPore',nPore,'Porosite',phi0, ...
        'Sv_par_voxel',Sv0,'Kref_LBM_lus',Kref, ...
        'NComposantes',nC0,'PercolationXYZ',perco0, ...
        'FracVolumeTraversantXYZ',frac0, ...
        'PartPlusGrandeComposante',grand0);
    fprintf('Connectivite originale (6) : [%d %d %d] ; %d composantes.\n\n', ...
        perco0(1),perco0(2),perco0(3),nC0);

    % Valeurs officielles de couverture et IoU (incluant les voxels
    % virtuels solides hors image), telles que publiees par le benchmark.
    fichierStats = fullfile(repertoireSelection,'statistiques_TC.csv');
    if isfile(fichierStats)
        stats = readtable(fichierStats);
    else
        warning('statistiques_TC.csv absent : colonnes benchmark = NaN.');
        stats = table();
    end

    noms = {'I6','I18','I26','E6','E18','E26','F'};
    sorties = struct([]);
    repResultats = fullfile(dossierBenchmark,'VALIDATION_PHYSIQUE_SANS_LBM');
    if ~isfolder(repResultats),mkdir(repResultats);end

    nLigne = 0;
    for m = 1:numel(noms)
        nom = noms{m};
        for it = 1:numel(TDs)
            td = TDs(it);
            for ic = 1:numel(TCs)
                tc = TCs(ic);
                f = fullfile(repertoireSelection, ...
                    sprintf('selection_%s_TD%02d_TC%02d.mat',nom,td,tc));
                assert(isfile(f),'Selection manquante : %s',f);
                D = load(f,'centresSelectionnes','rayonsSelectionnes');
                assert(isfield(D,'centresSelectionnes') && ...
                       isfield(D,'rayonsSelectionnes'), ...
                       'Variables inattendues dans %s',f);
                C = double(D.centresSelectionnes);
                R = double(D.rayonsSelectionnes(:));
                clear D
                assert(size(C,1)==numel(R) && size(C,2)==3, ...
                    'Centres/rayons de dimensions incompatibles dans %s',f);
                assert(all(isfinite(C(:))) && all(isfinite(R)) && all(R>0), ...
                    'Coordonnees ou rayons invalides dans %s',f);
                nbBoules = numel(R);
                t0 = tic;
                U = union_boules_centres(C,R,sz);
                duree = toc(t0);
                clear C R

                % Les voxels hors image ne figurent pas dans U : les
                % indicateurs ci-dessous sont limites au meme domaine image.
                TP = nnz(U & P);
                FPinterne = nnz(U & ~P);
                nB = nnz(U);
                phiB = nB / nVox;
                surfB = surface_interne_voxel(U);
                SvB = surfB / nVox;
                couv = 100*TP/max(1,nPore);
                debInterne = 100*FPinterne/max(1,nB);
                iouInterne = 100*TP/max(1,nPore+FPinterne);
                fluxRetenu = 100*sum(vitesseAbs(U & P))/sommeVitesse;
                [nC,perco,fracPerco,fracGrand] = connectivite_6(U);

                % Modele KC : meme coefficient que l'image de reference,
                % obtenu en calant le modele sur sa permeabilite LBM.
                % K_B/K_orig = (phi_B/phi_orig)^3 * (Sv_orig/Sv_B)^2.
                if nB == 0
                    Kkc = 0;
                elseif Sv0>0 && SvB>0 && phi0>0
                    Kkc = Kref * (phiB/phi0)^3 * (Sv0/SvB)^2;
                else
                    Kkc = NaN;
                end
                erreurK = 100*abs(Kkc-Kref)/Kref;
                % Variante diagnostic conditionnee a la connectivite Z.
                % Ce n'est pas une loi de Darcy ni une nouvelle DNS.
                KkcAvecConnexion = Kkc;
                if ~perco(3), KkcAvecConnexion=0; end
                erreurKConnexion = 100*abs(KkcAvecConnexion-Kref)/Kref;

                couvBenchmark=NaN; debBenchmark=NaN; iouBenchmark=NaN;
                if ~isempty(stats)
                    trouve = strcmp(string(stats.Methode),string(nom)) & ...
                        stats.SeuilTD_pct==td & stats.SeuilTC_pct==tc;
                    assert(nnz(trouve)==1,'Ligne TD/TC manquante ou multiple dans stats_TC.');
                    couvBenchmark=stats.CouverturePore_pct(trouve);
                    debBenchmark=stats.DebordementUnion_pct(trouve);
                    iouBenchmark=stats.IoU_pct(trouve);
                    assert(abs(couvBenchmark-couv)<1e-5, ...
                        ['Couverture recalculee differente du benchmark. ' ...
                         'Verifier la convention des centres et les selections.']);
                end

                nLigne=nLigne+1;
                sorties(nLigne).Methode = string(nom); %#ok<AGROW>
                sorties(nLigne).TD_pct = td;
                sorties(nLigne).TC_pct = tc;
                sorties(nLigne).NBoules = nbBoules;
                sorties(nLigne).PorositeOriginale_pct = 100*phi0;
                sorties(nLigne).PorositeBoules_pct = 100*phiB;
                sorties(nLigne).ErreurPorosite_points = 100*(phiB-phi0);
                sorties(nLigne).SurfaceInterface_faces = surfB;
                sorties(nLigne).SvOriginal_par_voxel = Sv0;
                sorties(nLigne).SvBoules_par_voxel = SvB;
                sorties(nLigne).NComposantes_6 = nC;
                sorties(nLigne).PartGrandeComposante_pct = 100*fracGrand;
                sorties(nLigne).TraverseX = logical(perco(1));
                sorties(nLigne).TraverseY = logical(perco(2));
                sorties(nLigne).TraverseZ = logical(perco(3));
                sorties(nLigne).VolumeComposantesTraversantesZ_pct = 100*fracPerco(3);
                sorties(nLigne).CouverturePoreInterne_pct = couv;
                sorties(nLigne).DebordementInterne_pct = debInterne;
                sorties(nLigne).IoUInterne_pct = iouInterne;
                sorties(nLigne).FluxAbsUzOriginalConserve_pct = fluxRetenu;
                sorties(nLigne).Kref_LBM_lus = Kref;
                sorties(nLigne).K_KC_proxy_lus = Kkc;
                sorties(nLigne).Erreur_K_KC_proxy_pct = erreurK;
                sorties(nLigne).K_KC_connexionZ_proxy_lus = KkcAvecConnexion;
                sorties(nLigne).Erreur_K_KC_connexionZ_pct = erreurKConnexion;
                sorties(nLigne).CouvertureBenchmark_pct = couvBenchmark;
                sorties(nLigne).DebordementBenchmark_pct = debBenchmark;
                sorties(nLigne).IoUBenchmark_pct = iouBenchmark;
                sorties(nLigne).TempsVoxelisation_s = duree;

                fprintf('%-4s TD%02d TC%02d | boules=%6d | phi=%.4f%% | IoU=%.2f%% | flux=%.2f%% | Z=%d | Khat/Kref=%.4f | %.1f s\n', ...
                    nom,td,tc,sorties(nLigne).NBoules,100*phiB, ...
                    iouBenchmark,fluxRetenu,perco(3),Kkc/Kref,duree);

                if sauverMasques
                    save(fullfile(repResultats, ...
                        sprintf('union_%s_TD%02d_TC%02d.mat',nom,td,tc)), ...
                        'U','nom','td','tc','-v7.3');
                end
                clear U
            end
        end
    end

    T = struct2table(sorties);
    % Les lignes sont sauvegardees dans l'ordre standard des 7 methodes.
    if numel(TDs)==1 && numel(TCs)==1
        suffixe=sprintf('TD%02d_TC%02d',TDs,TCs);
    else
        suffixe='configurations_demandees';
    end
    sortieCSV=fullfile(repResultats,['proprietes_physiques_' suffixe '.csv']);
    sortieMAT=fullfile(repResultats,['proprietes_physiques_' suffixe '.mat']);
    writetable(T,sortieCSV);
    save(sortieMAT,'T','ref','-v7.3');

    % Classements differents : geom (IoU), flux de reference et erreur KC.
    fprintf('\nCLASSEMENT PAR IoU benchmark (plus grand = meilleur) :\n');
    afficheColonnes = {'Methode','TD_pct','TC_pct','NBoules', ...
        'IoUBenchmark_pct','FluxAbsUzOriginalConserve_pct', ...
        'TraverseZ','Erreur_K_KC_proxy_pct'};
    disp(sortrows(T(:,afficheColonnes),'IoUBenchmark_pct','descend'));
    fprintf('\nCLASSEMENT PAR ERREUR DU PROXY Kozeny-Carman :\n');
    disp(sortrows(T(:,afficheColonnes),'Erreur_K_KC_proxy_pct','ascend'));
    fprintf('\nCSV : %s\nMAT : %s\n',sortieCSV,sortieMAT);
    fprintf(['ATTENTION : K_KC_proxy_lus n''est pas une permeabilite DNS des boules.\n' ...
             'FluxAbsUzOriginalConserve est une couverture ponderee par le champ DNS original.\n']);
end

function U = union_boules_centres(C,R,sz)
% Reproduit la condition discrete du benchmark ; ne cree pas de grille
% virtuelle hors image car les indicateurs physiques sont intra-image.
    U=false(sz);
    for b=1:numel(R)
        c=C(b,:); r=R(b);
        a=1e-12*max(1,r);
        lo=max([1 1 1], ceil(c-r-a));
        hi=min(sz, floor(c+r+a));
        if any(lo>hi),continue;end
        xx=double(lo(1):hi(1));
        yy=double(lo(2):hi(2));
        [X,Y]=ndgrid(xx,yy);
        dxy=(X-c(1)).^2+(Y-c(2)).^2;
        tol=1e-12*max(1,r*r);
        for k=lo(3):hi(3)
            masque = dxy+(k-c(3))^2 <= r*r+tol;
            if any(masque(:))
                slab=U(lo(1):hi(1),lo(2):hi(2),k);
                slab(masque)=true;
                U(lo(1):hi(1),lo(2):hi(2),k)=slab;
            end
        end
    end
end

function A = surface_interne_voxel(M)
% Surface pore-solide INTERNE au volume en nombre de faces de voxel.
% Les faces exterieures du cube d'imagerie ne sont pas comptees.
    A = nnz(xor(M(1:end-1,:,:),M(2:end,:,:))) + ...
        nnz(xor(M(:,1:end-1,:),M(:,2:end,:))) + ...
        nnz(xor(M(:,:,1:end-1),M(:,:,2:end)));
end

function [nComp,perco,fracPerco,fracMax] = connectivite_6(M)
% Composantes a 6 voisins et ponts entre faces opposees X, Y, Z.
    CC=bwconncomp(M,6);  % Image Processing Toolbox
    nComp=CC.NumObjects;
    perco=false(1,3);
    fracPerco=zeros(1,3);
    fracMax=0;
    if nComp==0,return;end
    tailles=cellfun(@numel,CC.PixelIdxList);
    nb=sum(tailles);
    fracMax=max(tailles)/nb;
    L=labelmatrix(CC);
    for d=1:3
        switch d
            case 1
                a=unique(reshape(L(1,:,:),[],1));
                b=unique(reshape(L(end,:,:),[],1));
            case 2
                a=unique(reshape(L(:,1,:),[],1));
                b=unique(reshape(L(:,end,:),[],1));
            otherwise
                a=unique(reshape(L(:,:,1),[],1));
                b=unique(reshape(L(:,:,end),[],1));
        end
        ab=intersect(a,b);
        ab=double(ab(ab>0));
        perco(d)=~isempty(ab);
        if perco(d)
            fracPerco(d)=sum(tailles(ab))/nb;
        end
    end
end

function K = permeabilite_depuis_csv(f)
% Lecture robuste de la ligne 'permeability (lus), ...'.
    fid=fopen(f,'rt');
    assert(fid~=-1,'Impossible de lire %s',f);
    c=onCleanup(@() fclose(fid));
    K=NaN;
    while true
        line=fgetl(fid);
        if ~ischar(line),break;end
        parts=strsplit(line,',');
        if numel(parts)>=2 && startsWith(lower(strtrim(parts{1})),'permeability')
            K=str2double(strtrim(parts{2}));
            break;
        end
    end
    assert(isfinite(K) && K>0, ...
        'Impossible de trouver la permeabilite de reference dans %s',f);
end
