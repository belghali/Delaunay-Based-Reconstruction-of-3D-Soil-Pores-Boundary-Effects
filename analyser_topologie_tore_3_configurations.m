function T = analyser_topologie_tore_3_configurations(dossierTore)
% ANALYSER_TOPOLOGIE_TORE_3_CONFIGURATIONS
% Reutilise exactement la methode de "analyser_tore_128_3_points.m":
%  - union des spheres sur centres de voxels (grille avec marge de 32),
%  - beta0 = composantes connexes 26 du premier plan,
%  - beta2 = cavites de fond 6 qui ne touchent pas la boite,
%  - chi = V - E + F - C du complexe cubique de voxels fermes,
%  - beta1 = beta0 + beta2 - chi.
%
% AUCUN recalcul de Delaunay ou de selection des spheres.
% Les trois jeux TD/TC sont traites : 0/0, 5/5 et 10/10.
%
% Exemple :
%   T = analyser_topologie_tore_3_configurations('TORE_SYNTHETIQUE_128');
%
% Sorties :
%  ANALYSES_GEOMETRIQUES_IMPROVE/
%    04_topologie_tore_TD00_TD05_TD10.csv
%    04_topologie_tore_TD00_TD05_TD10.mat
%
% Requires Image Processing Toolbox (bwconncomp).

    if nargin < 1 || isempty(dossierTore)
        dossierTore = fullfile(pwd,'TORE_SYNTHETIQUE_128');
    end

    fRef = fullfile(dossierTore,'tore_analytique_128.mat');
    assert(isfile(fRef),'Fichier de reference introuvable : %s',fRef);
    A = load(fRef,'imgTore');
    assert(isfield(A,'imgTore'),'Variable imgTore introuvable : %s',fRef);
    ref128 = logical(A.imgTore);
    assert(isequal(size(ref128),[128 128 128]), ...
           'Le benchmark attend un tore de taille 128 x 128 x 128.');

    marge = 32;
    ref = false(size(ref128)+2*marge);
    ref((1:128)+marge,(1:128)+marge,(1:128)+marge) = ref128;
    sz = size(ref);

    repSel = fullfile(dossierTore,'BENCHMARK_TORE_128_NEUF', ...
        'TD_TC_step1_sequentiel_audit_k1');
    repOut = fullfile(dossierTore,'ANALYSES_GEOMETRIQUES_IMPROVE');
    assert(isfolder(repSel),'Dossier de selections introuvable : %s',repSel);
    if ~isfolder(repOut), mkdir(repOut); end

    autotest_topologie();
    refTop = topologie_cubique(ref);
    assert(refTop.Beta0_26 == 1 && refTop.Beta1 == 1 && ...
           refTop.Beta2 == 0 && refTop.Chi == 0, ...
           'Le tore de reference ne presente pas la topologie (1,1,0), chi=0.');

    methodes = {'I6','I18','I26','E6','E18','E26','F'};
    seuils = [0 5 10];

    % Ajouter le tore de reference en premier (TD/TC non applicable).
    lignes = struct([]);
    ligne = struct('Configuration',"Reference", ...
        'Methode',"Reference torus",'TD_pct',NaN,'TC_pct',NaN, ...
        'NBoules',NaN,'IoU_pct',100, ...
        'Beta0_26',refTop.Beta0_26,'Beta1',refTop.Beta1, ...
        'Beta2',refTop.Beta2,'Euler',refTop.Chi, ...
        'TopologiePreservee',true);
    lignes = ligne;

    fprintf('\nTORE ORIGINAL : beta0=%d beta1=%d beta2=%d chi=%d\n', ...
        refTop.Beta0_26,refTop.Beta1,refTop.Beta2,refTop.Chi);

    for s = 1:numel(seuils)
        v = seuils(s);
        fprintf('\n===== TD%02d / TC%02d =====\n',v,v);
        for j = 1:numel(methodes)
            nom = methodes{j};
            fSel = fullfile(repSel, ...
                sprintf('selection_%s_TD%02d_TC%02d.mat',nom,v,v));
            assert(isfile(fSel),'Selection introuvable : %s',fSel);
            S = load(fSel,'centresSelectionnes','rayonsSelectionnes');
            assert(isfield(S,'centresSelectionnes') && ...
                   isfield(S,'rayonsSelectionnes'), ...
                   'Variables centresSelectionnes ou rayonsSelectionnes absentes.');
            RR = double(S.rayonsSelectionnes(:));
            C = double(S.centresSelectionnes);

            if isempty(RR)
                % Topologie de l'ensemble vide : beta0=beta1=beta2=chi=0.
                % Contrairement aux distances de surface, ce n'est PAS NaN.
                topo = struct('Beta0_26',0,'Beta1',0,'Beta2',0,'Chi',0);
                IoU_pct = 0;
            else
                assert(size(C,2)==3 && size(C,1)==numel(RR), ...
                    'Centres/rayons incoherents dans %s',fSel);
                C = C + marge;
                verifier_marge(C,RR,sz,nom);
                U = union_boules(C,RR,sz);
                tp = nnz(U & ref);
                fp = nnz(U & ~ref);
                IoU_pct = 100 * tp / (nnz(ref)+fp);
                topo = topologie_cubique(U);
                clear U
            end

            ok = topo.Beta0_26==1 && topo.Beta1==1 && topo.Beta2==0;
            ligne = struct('Configuration',string(sprintf('TD%d/TC%d',v,v)), ...
                'Methode',string(nom),'TD_pct',v,'TC_pct',v, ...
                'NBoules',numel(RR),'IoU_pct',IoU_pct, ...
                'Beta0_26',topo.Beta0_26,'Beta1',topo.Beta1, ...
                'Beta2',topo.Beta2,'Euler',topo.Chi, ...
                'TopologiePreservee',ok);
            lignes(end+1) = ligne; %#ok<AGROW>

            fprintf('%-4s : %5d spheres | IoU=%7.3f%% | beta=(%d,%d,%d) | chi=%d\n', ...
                nom,numel(RR),IoU_pct, ...
                topo.Beta0_26,topo.Beta1,topo.Beta2,topo.Chi);
        end
    end

    T = struct2table(lignes);
    writetable(T,fullfile(repOut,'04_topologie_tore_TD00_TD05_TD10.csv'));
    save(fullfile(repOut,'04_topologie_tore_TD00_TD05_TD10.mat'), ...
        'T','refTop','-v7.3');
    fprintf('\nTableaux enregistres dans :\n%s\n',repOut);
end

function verifier_marge(C,R,sz,nom)
    if isempty(R),return;end
    epsR=1e-12*max(1,R);
    bas=ceil(C-R-epsR);
    haut=floor(C+R+epsR);
    assert(~any(bas(:)<1) && ~any(haut(:,1)>sz(1)) && ...
           ~any(haut(:,2)>sz(2)) && ~any(haut(:,3)>sz(3)), ...
       'Marge insuffisante pour %s : augmenter marge dans le script.',nom);
end

function ind=indices_sphere(c,r,sz)
    % Exactement le meme test voxel-centre du programme de benchmark,
    % sur une grille etendue, exterieur au tore = solide virtuel.
    tol=1e-12*max(1,r*r);
    lo=max([1 1 1],ceil(c-r-1e-12*max(1,r)));
    hi=min(sz,floor(c+r+1e-12*max(1,r)));
    if any(lo>hi), ind=zeros(0,1);return;end
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
        for k=lo(3):hi(3)
            ok=dxy+(k-c(3))^2<=r*r+tol;
            [ii,jj]=ind2sub(size(dxy),find(ok));
            if ~isempty(ii)
                ind=[ind;sub2ind(sz,lo(1)-1+ii,lo(2)-1+jj, ...
                    k*ones(numel(ii),1))]; %#ok<AGROW>
            end
        end
    end
end

function U=union_boules(C,R,sz)
    U=false(sz);
    for b=1:numel(R)
        ind=indices_sphere(C(b,:),R(b),sz);
        U(ind)=true;
    end
end


function t=topologie_cubique(U)
    % Cubes fermes (voxels) : notion de connexite par faces/aretes/sommets.
    % Beta0 cubique = composantes 26. Cavites = composantes du complement
    % en voisinage 6 qui NE touchent pas le bord de la boite.
    % Beta1 = Beta0 + Beta2 - chi (identite d'Euler-Poincare).
    sz=size(U); a=sz(1);b=sz(2);c=sz(3);
    t.Beta0_6=bwconncomp(U,6).NumObjects;
    t.Beta0_26=bwconncomp(U,26).NumObjects;
    CC=bwconncomp(~U,6);
    cav=0;
    for j=1:CC.NumObjects
        [x,y,z]=ind2sub(sz,CC.PixelIdxList{j});
        if ~any(x==1 | x==a | y==1 | y==b | z==1 | z==c)
            cav=cav+1;
        end
    end
    t.Beta2=cav;

    % Nombre de cellules d'un complexe cubique (cubes fermes).
    v=false(a+1,b+1,c+1);
    for i=0:1,for j=0:1,for k=0:1
        v((1:a)+i,(1:b)+j,(1:c)+k)= ...
            v((1:a)+i,(1:b)+j,(1:c)+k) | U;
    end,end,end
    V=nnz(v); clear v

    ex=false(a,b+1,c+1);
    for j=0:1,for k=0:1
        ex(:,(1:b)+j,(1:c)+k)=ex(:,(1:b)+j,(1:c)+k) | U;
    end,end
    ey=false(a+1,b,c+1);
    for i=0:1,for k=0:1
        ey((1:a)+i,:,(1:c)+k)=ey((1:a)+i,:,(1:c)+k) | U;
    end,end
    ez=false(a+1,b+1,c);
    for i=0:1,for j=0:1
        ez((1:a)+i,(1:b)+j,:)=ez((1:a)+i,(1:b)+j,:) | U;
    end,end
    E=nnz(ex)+nnz(ey)+nnz(ez);clear ex ey ez

    fx=false(a+1,b,c);
    fx(1:a,:,:)=fx(1:a,:,:) | U;
    fx(2:a+1,:,:)=fx(2:a+1,:,:) | U;
    fy=false(a,b+1,c);
    fy(:,1:b,:)=fy(:,1:b,:) | U;
    fy(:,2:b+1,:)=fy(:,2:b+1,:) | U;
    fz=false(a,b,c+1);
    fz(:,:,1:c)=fz(:,:,1:c) | U;
    fz(:,:,2:c+1)=fz(:,:,2:c+1) | U;
    F=nnz(fx)+nnz(fy)+nnz(fz);
    t.Chi=double(V)-double(E)+double(F)-double(nnz(U));
    t.Beta1=t.Beta0_26+t.Beta2-t.Chi;
    if t.Beta1<0
        warning('Beta1 negatif : verifier les conventions de topologie digitale.');
    end
end


function autotest_topologie()
    n=48;c=(n+1)/2;
    [x,y,z]=ndgrid(1:n,1:n,1:n);
    sphere=(x-c).^2+(y-c).^2+(z-c).^2<=9^2;
    bs=topologie_cubique(sphere);
    assert(bs.Beta0_26==1 && bs.Beta1==0 && bs.Beta2==0 && bs.Chi==1, ...
        'Autotest sphere topologique echoue.');
    tore=(hypot(x-c,y-c)-14).^2+(z-c).^2<=6^2;
    bt=topologie_cubique(tore);
    assert(bt.Beta0_26==1 && bt.Beta1==1 && bt.Beta2==0 && bt.Chi==0, ...
        'Autotest tore topologique echoue.');
    fprintf('Auto-tests topologiques OK : sphere (1,0,0), tore (1,1,0).\n');
end
