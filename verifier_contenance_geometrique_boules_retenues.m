function [statsGeo, statsSeuils] = verifier_contenance_geometrique_boules_retenues(img, dossier, opts)
% VERIFIER_CONTENANCE_GEOMETRIQUE_BOULES_RETENUES
% Controle deterministe de penetration des boules dans le SOLIDE et hors
% de la boite image, pour les SEPT methodes deja auditees (TD discret <=20%).
%
% Modele continu : chaque voxel pore (255/1) est un cube ferme de cote 1
% centre en (i,j,k), avec i,j,k indices MATLAB. Le solide est l'union des
% cubes complementaires + l'exterieur de [0.5,nx+0.5]x[...]x[...].
%
% Critere : une sphere de centre c et rayon r penetre dans un cube solide
% de centre v si sum(max(abs(c-v)-0.5,0).^2) < r^2, a tolerance pres.
% Un contact tangent n'est pas une penetration volumique.
%
% NOTE : Cette fonction DETECTE la penetration, mais ne calcule PAS le
% volume d'intersection. Le TD volumique reste estime dans l'audit MC.
%
% ENTREE : dossier/audit_TD_volume_retenues/audit_retenues_<M>.mat
%         contient les boules centre-pore et TD-discret <= 20%.
% opts.SeuilTDDiscretMax : 0.20 (peut etre diminue sans recalculer)
% opts.ToleranceRelative: 1e-10 (tolerance en unite de voxel)
% opts.VerboseTousLes   : 5000
%
% SORTIES : dossier/controle_geometrique_TD0/
%  controle_<M>.mat et statistiques_controle_geometrique.csv/.mat
%  statsSeuils : detail du controle par seuil TD discret.
%
% UTILISATION :
%   opts=struct('SeuilTDDiscretMax',0.20,'ToleranceRelative',1e-10);
%   [statsGeo,statsSeuils] = ...
%       verifier_contenance_geometrique_boules_retenues(img,'resultats_delaunay_64',opts);

    if nargin < 2 || isempty(dossier), dossier = 'resultats_delaunay_64'; end
    if nargin < 3 || isempty(opts), opts = struct(); end
    opts = option_defaut(opts,'SeuilTDDiscretMax',0.20);
    opts = option_defaut(opts,'ToleranceRelative',1e-10);
    opts = option_defaut(opts,'VerboseTousLes',5000);
    opts = option_defaut(opts,'SeuilsTDDiscret',[0 0.05 0.10 0.20]);

    validateattributes(img,{'logical','numeric'},{'3d','real','nonempty'},mfilename,'img');
    if any(img(:) ~= 0 & img(:) ~= 1 & img(:) ~= 255)
        error('Image attendue : 255 (ou 1) = pore, 0 = solide.');
    end
    validateattributes(opts.SeuilTDDiscretMax,{'numeric'}, ...
        {'scalar','real','finite','>=',0,'<=',0.20});
    validateattributes(opts.ToleranceRelative,{'numeric'}, ...
        {'scalar','real','finite','positive','<=',1e-5});
    validateattributes(opts.VerboseTousLes,{'numeric'}, ...
        {'scalar','integer','nonnegative'});
    validateattributes(opts.SeuilsTDDiscret,{'numeric'}, ...
        {'vector','finite','>=',0,'<=',0.20});

    S = (img ~= 0);
    sz = size(S);
    entre = fullfile(dossier,'audit_TD_volume_retenues');
    sortie = fullfile(dossier,'controle_geometrique_TD0');
    if ~isfolder(entre)
        error('Dossier introuvable : %s. Lancer audit TD volume retenues.',entre);
    end
    if ~exist(sortie,'dir'), mkdir(sortie); end
    noms = {'I6','I18','I26','E6','E18','E26','F'};
    seuils = unique(double(opts.SeuilsTDDiscret(:)'),'sorted');
    nm = numel(noms); ns = numel(seuils);

    Methode = noms(:);
    N_Boules = zeros(nm,1);
    N_Inscrites = zeros(nm,1);
    N_DebordentSolide = zeros(nm,1);
    N_DebordentDomaine = zeros(nm,1);
    N_DebordentDeux = zeros(nm,1);
    N_ContactNumerique = zeros(nm,1);
    PartInscrites_pct = zeros(nm,1);
    TempsControle_s = zeros(nm,1);

    ms = cell(nm*ns,1);
    ss = zeros(nm*ns,1);
    ne = zeros(nm*ns,1);
    ni = zeros(nm*ns,1);
    nd = zeros(nm*ns,1);

    fprintf('\n===== CONTROLE GEOMETRIQUE DETERMINISTE : BOULES RETENUES =====\n');
    fprintf('Pore=%d voxels | volume=%s | tolerance relative=%.1e\n', ...
        nnz(S),mat2str(sz),opts.ToleranceRelative);

    for m = 1:nm
        nom = noms{m};
        fichier = fullfile(entre,['audit_retenues_' nom '.mat']);
        if ~isfile(fichier)
            error('Fichier introuvable : %s',fichier);
        end
        D = load(fichier,'centres','rayons','TDDiscret','TDVolume', ...
            'indicesOriginaux','indicesDansTD','meta');
        besoins = {'centres','rayons','TDDiscret','indicesOriginaux','indicesDansTD'};
        if ~all(isfield(D,besoins))
            error('Fichier audit incomplet : %s',fichier);
        end
        idx = find(D.TDDiscret(:) <= opts.SeuilTDDiscretMax + 1e-12);
        centres = double(D.centres(idx,:));
        rayons = double(D.rayons(idx));
        TDDiscret = double(D.TDDiscret(idx));
        indicesOriginaux = D.indicesOriginaux(idx);
        indicesDansTD = D.indicesDansTD(idx);
        if isfield(D,'TDVolume'), TDVolume = D.TDVolume(idx); else, TDVolume=[]; end
        N = numel(rayons);
        if size(centres,1) ~= N || size(centres,2) ~= 3 || ...
                any(~isfinite(centres(:))) || any(~isfinite(rayons)) || any(rayons<=0)
            error('Coordonnees ou rayons invalides pour %s.',nom);
        end
        if isfield(D,'meta') && isfield(D.meta,'SeuilTDDiscretMax') && ...
                D.meta.SeuilTDDiscretMax + 1e-12 < opts.SeuilTDDiscretMax
            error('%s : le fichier audit ne couvre pas le seuil exige.',nom);
        end

        DebordeSolide = false(N,1);
        DebordeDomaine = false(N,1);
        ContactNumerique = false(N,1);
        DistanceMinSolide = inf(N,1);
        MargeMinSolide = inf(N,1);
        fprintf('\n[%d/7] %s : %d boules a controler\n',m,nom,N);
        t = tic;
        for b = 1:N
            c = centres(b,:); r = rayons(b);
            td = opts.ToleranceRelative * max(1,r);
            td2 = opts.ToleranceRelative * max(1,r*r);

            % Penetration hors de la boite image.
            bas = c-r; haut = c+r;
            DebordeDomaine(b) = any(bas < 0.5-td | haut > double(sz)+0.5+td);
            toucheBoite = any(abs(bas-0.5)<=td | ...
                              abs(haut-(double(sz)+0.5))<=td);

            % Seuls les voxels solides avec cube susceptible de toucher
            % la sphere sont testes : bornes des centres entiers.
            lo = max([1 1 1], ceil(c-r-0.5));
            hi = min(double(sz), floor(c+r+0.5));
            toucheSolide = false;
            if all(lo <= hi)
                [I,J,K] = ndgrid(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3));
                id = sub2ind(sz,I(:),J(:),K(:));
                estSolide = ~S(id);
                if any(estSolide)
                    Q = double([I(estSolide),J(estSolide),K(estSolide)]);
                    Dxyz = max(abs(Q-c)-0.5,0);
                    d2 = sum(Dxyz.^2,2);
                    dm2 = min(d2);
                    DistanceMinSolide(b) = sqrt(dm2);
                    MargeMinSolide(b) = sqrt(dm2)-r;
                    DebordeSolide(b) = (dm2 < r*r-td2);
                    toucheSolide = any(abs(d2-r*r)<=td2);
                end
            end
            ContactNumerique(b) = toucheSolide || toucheBoite;
            if opts.VerboseTousLes > 0 && mod(b,opts.VerboseTousLes)==0
                fprintf('  %d/%d | %.1f s\n',b,N,toc(t));
            end
        end
        TempsControle_s(m) = toc(t);
        EstInscrite = ~DebordeSolide & ~DebordeDomaine;
        N_Boules(m) = N;
        N_Inscrites(m) = nnz(EstInscrite);
        N_DebordentSolide(m) = nnz(DebordeSolide);
        N_DebordentDomaine(m) = nnz(DebordeDomaine);
        N_DebordentDeux(m) = nnz(DebordeSolide & DebordeDomaine);
        N_ContactNumerique(m) = nnz(ContactNumerique);
        PartInscrites_pct(m) = 100 * nnz(EstInscrite)/max(1,N);

        for s=1:ns
            row = (m-1)*ns+s;
            test = TDDiscret <= seuils(s)+1e-12;
            ms{row} = nom;
            ss(row) = 100*seuils(s);
            ne(row) = nnz(test);
            ni(row) = nnz(test & EstInscrite);
            nd(row) = nnz(test & ~EstInscrite);
        end
        metadata = struct('Methode',nom,'SourceAudit',fichier, ...
            'ToleranceRelative',opts.ToleranceRelative, ...
            'SeuilTDDiscretMax',opts.SeuilTDDiscretMax, ...
            'Convention','Voxel pore = cube unite. Hors volume = solide.', ...
            'Nature','Test de penetration deterministe avec tolerance numerique; pas de volume TD.');
        save(fullfile(sortie,['controle_' nom '.mat']), ...
            'centres','rayons','TDDiscret','TDVolume', ...
            'DebordeSolide','DebordeDomaine','ContactNumerique', ...
            'DistanceMinSolide','MargeMinSolide','EstInscrite', ...
            'indicesOriginaux','indicesDansTD','metadata','-v7.3');
        fprintf('  Sans penetration detectee : %d / %d (%.2f%%)\n', ...
            N_Inscrites(m),N,PartInscrites_pct(m));
        fprintf('  Penetration solide=%d | hors domaine=%d | contacts proches=%d\n', ...
            N_DebordentSolide(m),N_DebordentDomaine(m),N_ContactNumerique(m));
        fprintf('  Temps : %.2f s\n',TempsControle_s(m));
    end

    statsGeo = table(Methode,N_Boules,N_Inscrites,PartInscrites_pct, ...
        N_DebordentSolide,N_DebordentDomaine,N_DebordentDeux, ...
        N_ContactNumerique,TempsControle_s);
    statsSeuils = table(ms,ss,ne,ni,nd, 'VariableNames', ...
        {'Methode','SeuilTDDiscret_pct','N_Evaluees','N_Inscrites','N_AvecPenetration'});
    writetable(statsGeo,fullfile(sortie,'statistiques_controle_geometrique.csv'));
    writetable(statsSeuils,fullfile(sortie,'statistiques_controle_par_seuil.csv'));
    save(fullfile(sortie,'statistiques_controle_geometrique.mat'), ...
        'statsGeo','statsSeuils','opts','-v7.3');
    fprintf('\n===== RECAPITULATIF CONTROLE GEOMETRIQUE =====\n');
    disp(statsGeo);
    fprintf('Resultats sauves dans : %s\n',sortie);
end

function opts=option_defaut(opts,nom,valeur)
    if ~isfield(opts,nom) || isempty(opts.(nom)), opts.(nom)=valeur; end
end
