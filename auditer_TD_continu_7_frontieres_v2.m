function [statsAudit, details] = auditer_TD_continu_7_frontieres_v2(img, dossier, opts)
% AUDITER_TD_CONTINU_7_FRONTIERES_V2
% Version avec selection ALEATOIRE REPRODUCTIBLE des spheres, apres tri
% canonique par rayon et coordonnees, sans dependance a l ordre en fichier.
% Les memes banques donnent alors le meme audit Monte Carlo.
% Audit independant du TD discret par une estimation volumique Monte Carlo.
%
% Le pore est l'union des cubes unitaires centres sur les voxels non nuls.
% Les cubes hors de l'image sont solides. Pour chaque sphere, des points
% sont tires UNIFORMEMENT dans son volume continu, puis attribues aux
% voxels par floor(coordonnee+0.5). L'estimateur obtenu approche :
%       TD_volume = volume(sphere \ pore_voxelise) / volume(sphere)
%
% Cet audit ne remplace pas le TD discret et ne modifie aucune banque.
% N.B. Precision statistique limitee par opts.NEchantillons ; pour les
% resultats finaux, augmenter NEchantillons et controler la convergence.
%
% UTILISATION :
% opts = struct('MaxBoules',1000,'NEchantillons',8192,'Graine',42);
% [statsAudit,details] = auditer_TD_continu_7_frontieres_v2( ...
%     img,'resultats_delaunay_64',opts);
%
% Entree : dossier/filtrage_TD/boules_TD_<methode>.mat
% Sortie : dossier/audit_TD_continu_v2/*.csv et *.mat

    if nargin < 2 || isempty(dossier)
        dossier = 'resultats_delaunay_64';
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts,'MaxBoules') || isempty(opts.MaxBoules)
        opts.MaxBoules = 1000;
    end
    if ~isfield(opts,'NEchantillons') || isempty(opts.NEchantillons)
        opts.NEchantillons = 8192;
    end
    if ~isfield(opts,'Graine') || isempty(opts.Graine)
        opts.Graine = 42;
    end

    validateattributes(img,{'numeric','logical'},{'real','nonempty','3d'});
    if any(img(:) ~= 0 & img(:) ~= 1 & img(:) ~= 255)
        error('Valeurs de img attendues : 0 et 255 (ou 0 et 1).');
    end
    validateattributes(opts.MaxBoules,{'numeric'}, ...
        {'scalar','integer','positive'});
    validateattributes(opts.NEchantillons,{'numeric'}, ...
        {'scalar','integer','>=',100});
    validateattributes(opts.Graine,{'numeric'}, ...
        {'scalar','integer','nonnegative'});

    S = img ~= 0;
    sz = size(S);
    entree = fullfile(dossier,'filtrage_TD');
    sortie = fullfile(dossier,'audit_TD_continu_v2');
    if ~exist(sortie,'dir'), mkdir(sortie); end

    % Echantillons communs aux methodes : directions uniformes sur S^2,
    % rayons selon U^(1/3), ce qui donne un tirage uniforme dans B(0,1).
    rng(opts.Graine,'twister');
    D = randn(opts.NEchantillons,3);
    D = D ./ sqrt(sum(D.^2,2));
    U = D .* (rand(opts.NEchantillons,1).^(1/3));

    methodes = {'I6','I18','I26','E6','E18','E26','F'};
    nM = numel(methodes);
    Methode = cell(nM,1);
    N_BoulesTestees = zeros(nM,1);
    MoyTD_Discret_pct = nan(nM,1);
    MoyTD_Volume_pct = nan(nM,1);
    MedTD_Volume_pct = nan(nM,1);
    BoulesVolumeTDPositif_pct = nan(nM,1);
    BoulesVolumeTD5_pct = nan(nM,1);
    BoulesVolumeTD20_pct = nan(nM,1);
    TempsAudit_s = zeros(nM,1);

    toutesMethodes = {};
    tousIndices = [];
    tousRayons = [];
    tousTDDiscret = [];
    tousTDVolume = [];
    banques = cell(nM,1); % cles geometriques pour comparaison directe

    fprintf('\n===== AUDIT TD DISCRET vs VOLUME =====\n');
    fprintf('Tirages par sphere : %d | Max spheres/methode : %d\n', ...
        opts.NEchantillons,opts.MaxBoules);
    fprintf('Pore : %d voxels | Dimensions : %s\n',nnz(S),mat2str(sz));

    for m = 1:nM
        nom = methodes{m};
        f = fullfile(entree,['boules_TD_' nom '.mat']);
        if ~isfile(f)
            error('Fichier TD introuvable : %s',f);
        end
        A = load(f,'centres','rayons','TD');
        centres = double(A.centres);
        rayons = double(A.rayons(:));
        TDdis = double(A.TD(:));
        if numel(rayons) ~= size(centres,1) || numel(TDdis) ~= numel(rayons)
            error('Dimensions incoherentes dans %s',f);
        end
        % Canonisation de la geometrie : ordre independant du fichier.
        % La quantification est coherente avec la tolerance des doublons.
        cles = round([rayons,centres] / 1e-7);
        [banques{m},ordre] = sortrows(cles,[1 2 3 4]);

        % Echantillon uniforme SANS remplacement dans la banque.
        % Meme graine par methode : banques identiques => memes spheres.
        nSel = min(opts.MaxBoules,numel(rayons));
        if nSel == 0
            error('%s : aucune sphere admissible.',nom);
        end
        rng(opts.Graine + 1,'twister');
        choix = randperm(numel(rayons),nSel);
        sel = ordre(choix(:));
        TDvol = zeros(nSel,1);
        t0 = tic;
        for b = 1:nSel
            c = centres(sel(b),:);
            r = rayons(sel(b));
            X = c + r*U;
            I = floor(X+0.5);  % voxel de centre entier contenant X
            valide = I(:,1)>=1 & I(:,1)<=sz(1) & ...
                     I(:,2)>=1 & I(:,2)<=sz(2) & ...
                     I(:,3)>=1 & I(:,3)<=sz(3);
            nPore = 0;
            if any(valide)
                ind = sub2ind(sz,I(valide,1),I(valide,2),I(valide,3));
                nPore = nnz(S(ind));
            end
            TDvol(b) = 1 - nPore/opts.NEchantillons;
        end
        TempsAudit_s(m) = toc(t0);
        Methode{m} = nom;
        N_BoulesTestees(m) = nSel;
        MoyTD_Discret_pct(m) = 100*mean(TDdis(sel));
        MoyTD_Volume_pct(m) = 100*mean(TDvol);
        MedTD_Volume_pct(m) = 100*median(TDvol);
        BoulesVolumeTDPositif_pct(m) = 100*mean(TDvol > 0);
        BoulesVolumeTD5_pct(m) = 100*mean(TDvol <= 0.05);
        BoulesVolumeTD20_pct(m) = 100*mean(TDvol <= 0.20);

        toutesMethodes = [toutesMethodes;repmat({nom},nSel,1)]; %#ok<AGROW>
        tousIndices = [tousIndices;sel(:)]; %#ok<AGROW>
        tousRayons = [tousRayons;rayons(sel)]; %#ok<AGROW>
        tousTDDiscret = [tousTDDiscret;TDdis(sel)]; %#ok<AGROW>
        tousTDVolume = [tousTDVolume;TDvol]; %#ok<AGROW>

        fprintf('%-3s : N=%4d | TD discret moyen=%7.3f%% | TD volume moyen=%7.3f%% | TD volume>0=%5.1f%% | %.2fs\n', ...
            nom,nSel,MoyTD_Discret_pct(m),MoyTD_Volume_pct(m), ...
            BoulesVolumeTDPositif_pct(m),TempsAudit_s(m));
    end

    % Test : E6, E18, E26 sont-elles geometriquement equivalentes
    % apres selection des centres poreux ?
    fprintf('\n===== TEST GEOMETRIQUE DES BANQUES EXTERIEURES =====\n');
    fprintf('E6 vs E18 identiques : %d\n',isequal(banques{4},banques{5}));
    fprintf('E18 vs E26 identiques : %d\n',isequal(banques{5},banques{6}));
    if isequal(banques{5},banques{6})
        fprintf('E18/E26 : ecart TDvol moyen = %.12g points de pourcentage\n', ...
            abs(MoyTD_Volume_pct(5)-MoyTD_Volume_pct(6)));
    end

    statsAudit = table(Methode,N_BoulesTestees,MoyTD_Discret_pct, ...
        MoyTD_Volume_pct,MedTD_Volume_pct, ...
        BoulesVolumeTDPositif_pct,BoulesVolumeTD5_pct, ...
        BoulesVolumeTD20_pct,TempsAudit_s);
    details = table(toutesMethodes,tousIndices,tousRayons, ...
        tousTDDiscret,tousTDVolume, ...
        'VariableNames',{'Methode','IndiceBoule','Rayon', ...
        'TDDiscret','TDVolumeEstime'});
    disp(statsAudit);
    writetable(statsAudit,fullfile(sortie,'stats_audit_TD.csv'));
    writetable(details,fullfile(sortie,'details_audit_TD.csv'));
    save(fullfile(sortie,'audit_TD.mat'),'statsAudit','details','opts','-v7.3');
    fprintf('Resultats enregistres dans %s\n',sortie);
end
