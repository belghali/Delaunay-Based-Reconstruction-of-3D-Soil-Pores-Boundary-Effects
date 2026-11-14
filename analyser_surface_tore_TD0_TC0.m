function T = analyser_surface_tore_TD0_TC0(dossierTore)
% ANALYSER_SURFACE_TORE_TD0_TC0
% Calcule les memes distances de surface que dans
% analyser_tore_128_3_points.m, mais UNIQUEMENT pour TD=0% et TC=0%.
%
% Utilisation :
%   T = analyser_surface_tore_TD0_TC0('TORE_SYNTHETIQUE_128');
%
% Fichiers requis :
%   TORE_SYNTHETIQUE_128/tore_analytique_128.mat
%       variables : imgTore, torus (R_vox, r_vox, centre_vox)
%   TORE_SYNTHETIQUE_128/BENCHMARK_TORE_128_NEUF/...
%       TD_TC_step1_sequentiel_audit_k1/selection_I6_TD00_TC00.mat
%       ... selection_I18_TD00_TC00.mat, ..., selection_F_TD00_TC00.mat
%       variables : centresSelectionnes, rayonsSelectionnes
%
% Les méthodes exterieures E6/E18/E26 n'ont aucune boule a TD0/TC0 :
% leurs distances sont NaN (non definies), et non zero.
%
% IMPORTANT : conserver le meme echantillonnage du tore analytique
% (240 x 120, uniformise en aire) et le meme domaine etendu (+32 voxels)
% que dans le calcul TD10/TC10 original.

    if nargin < 1 || isempty(dossierTore)
        dossierTore = fullfile(pwd,'TORE_SYNTHETIQUE_128');
    end
    fRef = fullfile(dossierTore,'tore_analytique_128.mat');
    assert(isfile(fRef),'Fichier de reference absent : %s',fRef);
    A = load(fRef,'imgTore','torus');
    assert(isfield(A,'imgTore') && isfield(A,'torus'), ...
        'imgTore et torus sont necessaires dans %s',fRef);
    ref128 = logical(A.imgTore);
    assert(isequal(size(ref128),[128 128 128]),'Tore attendu : 128 x 128 x 128');

    R = double(A.torus.R_vox);
    r = double(A.torus.r_vox);
    c0 = double(A.torus.centre_vox);
    assert(isscalar(c0),'Le centre du tore du script original doit etre scalaire.');
    marge = 32;
    ref = false(size(ref128)+2*marge);
    ref((1:128)+marge,(1:128)+marge,(1:128)+marge) = ref128;
    sz = size(ref);
    c = c0 + marge;

    dossierSel = fullfile(dossierTore,'BENCHMARK_TORE_128_NEUF', ...
        'TD_TC_step1_sequentiel_audit_k1');
    dossierOut = fullfile(dossierTore,'ANALYSES_GEOMETRIQUES_IMPROVE');
    if ~isfolder(dossierOut), mkdir(dossierOut); end

    % Echantillonnage analytique IDENTIQUE au script TD10/TC10.
    ptsAna = echantillons_surface_tore(R,r,c,240,120);
    methodes = {'I6','I18','I26','E6','E18','E26','F'};

    Methode = strings(7,1);
    NBoules = zeros(7,1);
    IoU_pct = zeros(7,1);
    Statut = strings(7,1);
    D_faces_vers_tore_moy_vox = nan(7,1);
    D_faces_vers_tore_RMSE_vox = nan(7,1);
    D_faces_vers_tore_P95_vox = nan(7,1);
    D_tore_vers_faces_moy_APPROX_vox = nan(7,1);
    D_sym_moy_APPROX_vox = nan(7,1);

    for k = 1:numel(methodes)
        nom = methodes{k};
        Methode(k) = string(nom);
        fSel = fullfile(dossierSel, ...
            sprintf('selection_%s_TD00_TC00.mat',nom));
        assert(isfile(fSel),'Fichier de selection absent : %s',fSel);
        S = load(fSel,'centresSelectionnes','rayonsSelectionnes');
        assert(isfield(S,'centresSelectionnes') && ...
               isfield(S,'rayonsSelectionnes'), ...
               'Variables centresSelectionnes/rayonsSelectionnes manquantes : %s',fSel);

        RR = double(S.rayonsSelectionnes(:));
        NBoules(k) = numel(RR);

        if isempty(RR)
            % Reconstruction vide : aucune surface reconstruite, mesures N/A.
            IoU_pct(k) = 0;
            Statut(k) = "vide";
            fprintf('%-4s : 0 sphere - distances non definies (NaN).\n',nom);
            continue;
        end

        C = double(S.centresSelectionnes);
        assert(size(C,2)==3 && size(C,1)==numel(RR), ...
               'Centres et rayons incompatibles : %s',fSel);
        C = C + marge;
        verifier_marge(C,RR,sz,nom);

        U = union_boules(C,RR,sz);
        tp = nnz(U & ref);
        IoU_pct(k) = 100 * tp / max(1,nnz(U | ref));
        faces = points_centres_faces(U);
        if isempty(faces)
            Statut(k) = "sans surface";
            continue;
        end

        % Distances reconstruit -> tore (SDF analytique exacte aux centres de faces).
        dRT = abs(sdf_tore(faces,R,r,c));
        % Distances tore -> reconstruit (echantillonnage analytique en aire).
        dTR = distances_surface_vers_faces(ptsAna,faces,sz);

        D_faces_vers_tore_moy_vox(k) = mean(dRT);
        D_faces_vers_tore_RMSE_vox(k) = sqrt(mean(dRT.^2));
        D_faces_vers_tore_P95_vox(k) = percentile_local(dRT,95);
        D_tore_vers_faces_moy_APPROX_vox(k) = mean(dTR);
        D_sym_moy_APPROX_vox(k) = ...
            (D_faces_vers_tore_moy_vox(k) + ...
             D_tore_vers_faces_moy_APPROX_vox(k))/2;
        Statut(k) = "calcule";
        fprintf('%-4s : %5d spheres | IoU=%7.3f%% | ',nom,NBoules(k),IoU_pct(k));
        fprintf('mean R->T=%7.3f | RMSE=%7.3f | P95=%7.3f | ', ...
            D_faces_vers_tore_moy_vox(k), ...
            D_faces_vers_tore_RMSE_vox(k), ...
            D_faces_vers_tore_P95_vox(k));
        fprintf('mean T->R~=%7.3f | sym~=%7.3f vox\n', ...
            D_tore_vers_faces_moy_APPROX_vox(k), ...
            D_sym_moy_APPROX_vox(k));
    end

    T = table(Methode,NBoules,IoU_pct,Statut, ...
        D_faces_vers_tore_moy_vox, ...
        D_faces_vers_tore_RMSE_vox, ...
        D_faces_vers_tore_P95_vox, ...
        D_tore_vers_faces_moy_APPROX_vox, ...
        D_sym_moy_APPROX_vox);

    disp(T);

    % Controle de coherence avec les nombres publies pour TD0/TC0.
    attendusN = [481;2112;2572;0;0;0;1399];
    attendusIoU = [67.651;72.070;71.419;0;0;0;74.623];
    if any(NBoules ~= attendusN)
        warning('Nombre de spheres different du tableau TD0/TC0 : verifier la banque.');
    end
    if any(abs(IoU_pct-attendusIoU) > 0.02)
        warning('IoU differente du tableau TD0/TC0 : verifier le domaine ou les selections.');
    end

    fichierCsv = fullfile(dossierOut,'03_surface_tore_TD00_TC00.csv');
    fichierMat = fullfile(dossierOut,'03_surface_tore_TD00_TC00.mat');
    writetable(T,fichierCsv);
    save(fichierMat,'T');
    fprintf('\nResultats TD0/TC0 :\n%s\n%s\n',fichierCsv,fichierMat);
end

function verifier_marge(C,R,sz,nom)
    margeEps = 1e-12*max(1,R);
    bas = ceil(C-R-margeEps);
    haut = floor(C+R+margeEps);
    assert(~any(bas(:)<1) && ...
           ~any(haut(:,1)>sz(1)) && ...
           ~any(haut(:,2)>sz(2)) && ...
           ~any(haut(:,3)>sz(3)), ...
           'Marge insuffisante pour %s : augmenter marge.',nom);
end

function U = union_boules(C,R,sz)
    U = false(sz);
    for j = 1:numel(R)
        ind = indices_sphere(C(j,:),R(j),sz);
        U(ind) = true;
    end
end

function ind = indices_sphere(c,r,sz)
    % Meme test d'appartenance des CENTRES de voxels que le benchmark.
    tol = 1e-12*max(1,r*r);
    lo = max([1 1 1],ceil(c-r-1e-12*max(1,r)));
    hi = min(sz,floor(c+r+1e-12*max(1,r)));
    if any(lo>hi), ind=zeros(0,1); return; end
    x=(lo(1):hi(1))-c(1);
    y=(lo(2):hi(2))-c(2);
    z=(lo(3):hi(3))-c(3);
    if prod(hi-lo+1)<=350000
        d2=reshape(x.^2,[],1,1)+reshape(y.^2,1,[],1)+reshape(z.^2,1,1,[]);
        [ii,jj,kk]=ind2sub(hi-lo+1,find(d2<=r*r+tol));
        ind=sub2ind(sz,lo(1)-1+ii,lo(2)-1+jj,lo(3)-1+kk);
    else
        [xx,yy]=ndgrid(lo(1):hi(1),lo(2):hi(2));
        dxy=(xx-c(1)).^2+(yy-c(2)).^2;
        ind=zeros(0,1);
        for kk=lo(3):hi(3)
            ok=dxy+(kk-c(3))^2<=r*r+tol;
            [ii,jj]=ind2sub(size(dxy),find(ok));
            if ~isempty(ii)
                ind=[ind;sub2ind(sz,lo(1)-1+ii,lo(2)-1+jj, ...
                    kk*ones(numel(ii),1))]; %#ok<AGROW>
            end
        end
    end
end

function p = points_centres_faces(U)
    % Centres de faces aux interfaces U / fond, sur les 3 axes.
    sz=size(U); a=sz(1);b=sz(2);c=sz(3);
    ax=cat(1,false(1,b,c),U); bx=cat(1,U,false(1,b,c));
    [i,j,k]=ind2sub(size(ax),find(xor(ax,bx)));
    px=[double(i)-0.5,double(j),double(k)];
    ay=cat(2,false(a,1,c),U); by=cat(2,U,false(a,1,c));
    [i,j,k]=ind2sub(size(ay),find(xor(ay,by)));
    py=[double(i),double(j)-0.5,double(k)];
    az=cat(3,false(a,b,1),U); bz=cat(3,U,false(a,b,1));
    [i,j,k]=ind2sub(size(az),find(xor(az,bz)));
    pz=[double(i),double(j),double(k)-0.5];
    p=[px;py;pz];
end

function d = sdf_tore(P,R,r,c)
    rho=hypot(P(:,1)-c,P(:,2)-c);
    d=hypot(rho-R,P(:,3)-c)-r;
end

function pts = echantillons_surface_tore(R,r,c,nu,nv)
    % Original TD10/TC10 : echantillonnage uniforme en AIRE.
    u=2*pi*((0:nu-1)+0.5)/nu;
    q=((0:nv-1)+0.5)/nv;
    v=2*pi*q;
    for k=1:12
        v=v-(R*v+r*sin(v)-2*pi*R*q)./(R+r*cos(v));
    end
    [uu,vv]=ndgrid(u,v);
    pts=[c+(R+r*cos(vv(:))).*cos(uu(:)), ...
         c+(R+r*cos(vv(:))).*sin(uu(:)), ...
         c+r*sin(vv(:))];
end

function d = distances_surface_vers_faces(P,faces,sz)
    % Identique a l'original : voisins centres de faces les plus proches.
    if isempty(faces),d=inf(size(P,1),1);return;end
    if exist('knnsearch','file')==2
        [~,d]=knnsearch(faces,P,'K',1);
    else
        % Repli du code original sans Statistics Toolbox : EDT demi-voxel.
        dim=2*double(sz)+2;
        B=false(dim);
        pos=round(2*faces)+1;
        B(sub2ind(dim,pos(:,1),pos(:,2),pos(:,3)))=true;
        D=0.5*bwdist(B);
        p=2*P+1;
        d=double(interpn(D,p(:,1),p(:,2),p(:,3),'linear',NaN));
        assert(all(isfinite(d)), ...
            'Points d''echantillonnage hors grille.');
    end
end

function p = percentile_local(x,pc)
    % Meme interpolation de percentile que l'analyse TD10/TC10.
    x=sort(double(x(isfinite(x))));
    if isempty(x),p=NaN;return;end
    i=1+(numel(x)-1)*pc/100;
    a=floor(i);b=ceil(i);
    p=x(a)+(i-a)*(x(b)-x(a));
end
