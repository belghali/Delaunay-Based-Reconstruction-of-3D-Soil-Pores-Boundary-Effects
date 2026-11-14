function [statsAudit, statsSeuils] = auditer_TD_volume_boules_retenues(img, dossier, opts)
% AUDITER_TD_VOLUME_BOULES_RETENUES
% Audit Monte Carlo du TD volumique UNIQUEMENT sur les boules retenues
% a l'etape 5 : centre strictement dans le pore ET TD discret <= seuil.
%
% Entrees :
%   img       : image 3D, 255 = pore, 0 = solide (ou 1/0)
%   dossier   : 'resultats_delaunay_64' par defaut
%   opts.SeuilTDDiscretMax : 0.20 par defaut (fraction, pas pourcentage)
%   opts.NEchantillons      : 8192 par sphere
%   opts.Graine             : 42
%   opts.SeuilsTDVol        : [0 0.05 0.10 0.20]
%   opts.VerboseTousLes     : 5000
%
% Source : dossier/filtrage_TD/boules_TD_<methode>.mat
%    Le fichier contient deja seulement les centres admis par le filtrage
%    de l'etape 5. Cette fonction exige ici des centres STRICTEMENT poreux.
%    Elle retient en outre les boules dont TD <= SeuilTDDiscretMax.
%
% Sortie : dossier/audit_TD_volume_retenues/
%   audit_retenues_<methode>.mat (une estimation TD par boule + IC de Wilson)
%   statistiques_audit_retenues.csv
%   statistiques_seuils_TDvol.csv
%
% Remarque : le TD volumique est estime par Monte Carlo dans la sphere
% continue intersectee avec l'union des cubes poreux de cote 1.
% Un TD estime a zero n'est PAS une preuve de contenance exacte.

    if nargin < 2 || isempty(dossier)
        dossier = 'resultats_delaunay_64';
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = defaut(opts,'SeuilTDDiscretMax',0.20);
    opts = defaut(opts,'NEchantillons',8192);
    opts = defaut(opts,'Graine',42);
    opts = defaut(opts,'SeuilsTDVol',[0 0.05 0.10 0.20]);
    opts = defaut(opts,'VerboseTousLes',5000);

    validateattributes(img,{'numeric','logical'},{'real','nonempty','3d'},mfilename,'img');
    if any(~isfinite(double(img(:)))) || any(img(:)~=0 & img(:)~=1 & img(:)~=255)
        error('Image attendue : 255/1 = pore, 0 = solide.');
    end
    validateattributes(opts.SeuilTDDiscretMax,{'numeric'}, ...
        {'scalar','real','finite','>=',0,'<=',1});
    validateattributes(opts.NEchantillons,{'numeric'}, ...
        {'scalar','integer','>=',100});
    validateattributes(opts.Graine,{'numeric'}, ...
        {'scalar','integer','nonnegative'});
    validateattributes(opts.SeuilsTDVol,{'numeric'}, ...
        {'vector','nonempty','finite','>=',0,'<=',1});
    validateattributes(opts.VerboseTousLes,{'numeric'}, ...
        {'scalar','integer','nonnegative'});

    opts.SeuilsTDVol = unique(double(opts.SeuilsTDVol(:)'),'sorted');
    S = (img ~= 0);
    sz = size(S);
    cheminEntree = fullfile(dossier,'filtrage_TD');
    cheminSortie = fullfile(dossier,'audit_TD_volume_retenues');
    if ~exist(cheminSortie,'dir'), mkdir(cheminSortie); end

    methodes = {'I6','I18','I26','E6','E18','E26','F'};
    nMethodes = numel(methodes);
    nSeuils = numel(opts.SeuilsTDVol);
    Nmc = double(opts.NEchantillons);
    z95 = 1.959963984540054;

    % Tirages partages par toutes les boules et toutes les methodes.
    % Meme sphere => meme TD Monte Carlo independamment de l'ordre du fichier.
    etatRng = rng;
    cleanup = onCleanup(@() rng(etatRng)); %#ok<NASGU>
    rng(opts.Graine,'twister');
    a = rand(Nmc,3);
    zz = 2*a(:,1)-1;
    theta = 2*pi*a(:,2);
    rho = a(:,3).^(1/3);
    uxy = rho.*sqrt(max(0,1-zz.^2));
    U = [uxy.*cos(theta), uxy.*sin(theta), rho.*zz];
    clear a zz theta rho uxy

    Methode = methodes(:);
    N_CentresAdmis = zeros(nMethodes,1);
    N_RetenuesTDDiscret = zeros(nMethodes,1);
    N_ExcluesTDDiscret = zeros(nMethodes,1);
    MoyTDDiscret_pct = nan(nMethodes,1);
    MoyTDVolume_pct = nan(nMethodes,1);
    MedTDVolume_pct = nan(nMethodes,1);
    PartTDVolumePositif_pct = nan(nMethodes,1);
    TempsAudit_s = zeros(nMethodes,1);

    MethodeSeuil = cell(nMethodes*nSeuils,1);
    SeuilTDVolume_pct = zeros(nMethodes*nSeuils,1);
    N_Evaluees = zeros(nMethodes*nSeuils,1);
    N_EstimeesSousSeuil = zeros(nMethodes*nSeuils,1);
    N_IC95EntierSousSeuil = zeros(nMethodes*nSeuils,1);
    N_IC95EntierSurSeuil = zeros(nMethodes*nSeuils,1);
    N_IC95Indecises = zeros(nMethodes*nSeuils,1);

    fprintf('\n===== AUDIT TD VOLUMIQUE : BOULES RETENUES UNIQUEMENT =====\n');
    fprintf('Volume : %s, %d voxels poreux\n',mat2str(sz),nnz(S));
    fprintf('Selection : centre pore et TD discret <= %.1f %%\n', ...
        100*opts.SeuilTDDiscretMax);
    fprintf('Tirages Monte Carlo par boule : %d\n',Nmc);

    for m = 1:nMethodes
        nom = methodes{m};
        fichier = fullfile(cheminEntree,['boules_TD_' nom '.mat']);
        if ~isfile(fichier), error('Fichier introuvable : %s',fichier); end
        D = load(fichier,'centres','rayons','TD', ...
                 'classeCentres','indicesOriginaux');
        if ~all(isfield(D,{'centres','rayons','TD','classeCentres','indicesOriginaux'}))
            error('%s : fichier TD incomplet.',nom);
        end
        Nb = numel(D.rayons);
        if size(D.centres,1) ~= Nb || size(D.centres,2) ~= 3 || ...
                numel(D.TD) ~= Nb || numel(D.classeCentres) ~= Nb || ...
                numel(D.indicesOriginaux) ~= Nb
            error('%s : dimensions incoherentes.',nom);
        end
        if any(D.classeCentres(:) ~= 1)
            error(['%s : banque TD contenant des centres non strictement ' ...
                   'poreux. Reexecuter filtrer_boules_TD_7_frontieres ' ...
                   'avec opts.CentresAdmis = ''pore''.'],nom);
        end
        N_CentresAdmis(m) = Nb;
        idx = find(D.TD(:) <= opts.SeuilTDDiscretMax + 1e-12);
        centres = double(D.centres(idx,:));
        rayons = double(D.rayons(idx));
        TDDiscret = double(D.TD(idx));
        indicesOriginaux = D.indicesOriginaux(idx);
        indicesDansTD = uint32(idx);
        N = numel(rayons);
        if any(~isfinite(centres(:))) || any(~isfinite(rayons)) || any(rayons <= 0)
            error('%s : boules non finies ou rayons invalides.',nom);
        end
        N_RetenuesTDDiscret(m) = N;
        N_ExcluesTDDiscret(m) = Nb-N;
        fprintf('\n[%d/%d] %s : %d / %d boules a auditer\n', ...
            m,nMethodes,nom,N,Nb);

        nTiragesHorsPore = zeros(N,1,'uint32');
        TDVolume = zeros(N,1);
        t = tic;
        for b = 1:N
            X = U*rayons(b) + centres(b,:);
            V = floor(X + 0.5); % cube centre sur voxel [i-.5,i+.5]
            dans = V(:,1)>=1 & V(:,1)<=sz(1) & ...
                   V(:,2)>=1 & V(:,2)<=sz(2) & ...
                   V(:,3)>=1 & V(:,3)<=sz(3);
            nPore = 0;
            if any(dans)
                indicesVoxels = sub2ind(sz,V(dans,1),V(dans,2),V(dans,3));
                nPore = nnz(S(indicesVoxels));
            end
            nTiragesHorsPore(b) = uint32(Nmc-nPore);
            TDVolume(b) = double(nTiragesHorsPore(b))/Nmc;
            if opts.VerboseTousLes > 0 && mod(b,opts.VerboseTousLes)==0
                fprintf('  %d/%d | %.1f s\n',b,N,toc(t));
            end
        end
        TempsAudit_s(m) = toc(t);

        den = 1+z95^2/Nmc;
        p0 = (TDVolume+z95^2/(2*Nmc))/den;
        dp = z95/den * sqrt(TDVolume.*(1-TDVolume)/Nmc + z95^2/(4*Nmc^2));
        IC95_Bas = max(0,p0-dp);
        IC95_Haut = min(1,p0+dp);

        if N > 0
            MoyTDDiscret_pct(m) = 100*mean(TDDiscret);
            MoyTDVolume_pct(m) = 100*mean(TDVolume);
            MedTDVolume_pct(m) = 100*median(TDVolume);
            PartTDVolumePositif_pct(m) = 100*mean(TDVolume>0);
        end

        for s = 1:nSeuils
            pos = (m-1)*nSeuils+s;
            th = opts.SeuilsTDVol(s);
            MethodeSeuil{pos} = nom;
            SeuilTDVolume_pct(pos) = th*100;
            N_Evaluees(pos) = N;
            N_EstimeesSousSeuil(pos) = nnz(TDVolume<=th);
            N_IC95EntierSousSeuil(pos) = nnz(IC95_Haut<=th);
            N_IC95EntierSurSeuil(pos) = nnz(IC95_Bas>th);
            N_IC95Indecises(pos) = N - N_IC95EntierSousSeuil(pos) ...
                - N_IC95EntierSurSeuil(pos);
        end

        meta = struct('Methode',nom, ...
            'SeuilTDDiscretMax',opts.SeuilTDDiscretMax, ...
            'NEchantillons',Nmc,'Graine',opts.Graine, ...
            'Origine',fichier, ...
            'Modele','Union des cubes poreux unitaires ; hors image = solide', ...
            'IC','Wilson 95%% : incertitude statistique, non certification');
        save(fullfile(cheminSortie,['audit_retenues_' nom '.mat']), ...
            'centres','rayons','TDDiscret','TDVolume', ...
            'nTiragesHorsPore','IC95_Bas','IC95_Haut', ...
            'indicesOriginaux','indicesDansTD','meta','-v7.3');

        fprintf('  TD discret moyen     : %.5f %%\n',MoyTDDiscret_pct(m));
        fprintf('  TD volumique moyen   : %.5f %%\n',MoyTDVolume_pct(m));
        fprintf('  TD volumique > 0     : %.2f %% des boules\n', ...
            PartTDVolumePositif_pct(m));
        fprintf('  Duree                : %.2f s\n',TempsAudit_s(m));
    end

    statsAudit = table(Methode,N_CentresAdmis,N_RetenuesTDDiscret, ...
        N_ExcluesTDDiscret,MoyTDDiscret_pct,MoyTDVolume_pct, ...
        MedTDVolume_pct,PartTDVolumePositif_pct,TempsAudit_s);
    statsSeuils = table(MethodeSeuil,SeuilTDVolume_pct,N_Evaluees, ...
        N_EstimeesSousSeuil,N_IC95EntierSousSeuil, ...
        N_IC95EntierSurSeuil,N_IC95Indecises);
    writetable(statsAudit,fullfile(cheminSortie,'statistiques_audit_retenues.csv'));
    writetable(statsSeuils,fullfile(cheminSortie,'statistiques_seuils_TDvol.csv'));
    save(fullfile(cheminSortie,'statistiques_audit_retenues.mat'), ...
        'statsAudit','statsSeuils','opts','-v7.3');
    fprintf('\n===== RESULTATS GLOBAUX =====\n');
    disp(statsAudit);
    fprintf('Resultats enregistres dans : %s\n',cheminSortie);
end

function opts = defaut(opts,champ,val)
    if ~isfield(opts,champ) || isempty(opts.(champ))
        opts.(champ) = val;
    end
end
