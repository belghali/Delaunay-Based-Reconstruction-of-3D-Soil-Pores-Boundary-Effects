function [Tfinal, Tbudget] = analyser_tore_128_3_points(dossierTore, tdPct, tcPct, budgets)
% ANALYSER_TORE_128_3_POINTS - etude pour IMPROVE 2027
%
% Trois experiences complementaires, sans recalculer Delaunay ni TD/TC :
%  (1) Topologie cubique : beta0, beta1, beta2 et Euler.
%  (2) Distance des faces reconstruites au tore analytique exact ;
%      distance retour echantillonnee (approximation vers centres des faces).
%  (3) IoU pour un MEME NOMBRE DE BOULES pour les 7 methodes,
%      depuis les banques APRES TD, AVANT TC, selon un protocole neuf
%      (rayon decroissant + contribution poreuse positive, sans TC).
%
% APPEL :
%  [Tfinal,Tbudget] = analyser_tore_128_3_points;
%  [Tfinal,Tbudget] = analyser_tore_128_3_points( ...
%       'TORE_SYNTHETIQUE_128',10,10,[25 50 100 250 500]);
%
% Necessite Image Processing Toolbox : bwconncomp, bwdist.
% Les distances sont en voxels ; un voxel est un cube de cote 1.
% Attention : un test TD=0 aux centres de voxels ne prouve PAS
% l'inscription continue de chaque sphere.

    if nargin<1 || isempty(dossierTore), dossierTore=fullfile(pwd,'TORE_SYNTHETIQUE_128'); end
    if nargin<2 || isempty(tdPct), tdPct=10; end
    if nargin<3 || isempty(tcPct), tcPct=10; end
    if nargin<4 || isempty(budgets), budgets=[25 50 100 250 500]; end

    validateattributes(tdPct,{'numeric'},{'scalar','nonnegative','finite'});
    validateattributes(tcPct,{'numeric'},{'scalar','nonnegative','finite'});
    assert(any(tdPct==[0 5 10]) && any(tcPct==[0 5 10]), ...
        'TD et TC doivent appartenir a {0,5,10}.');
    budgets=unique(double(budgets(:).'));
    assert(~isempty(budgets) && all(budgets>0 & mod(budgets,1)==0), ...
        'Les budgets doivent etre des entiers strictement positifs.');

    fichierTore=fullfile(dossierTore,'tore_analytique_128.mat');
    assert(isfile(fichierTore),'Image du tore introuvable : %s',fichierTore);
    A=load(fichierTore,'imgTore','torus');
    assert(isfield(A,'imgTore') && isfield(A,'torus'),'Variables imgTore / torus absentes.');
    ref128=logical(A.imgTore);
    sz128=size(ref128);
    assert(isequal(sz128,[128 128 128]),'Le programme a ete concu pour le tore 128^3.');
    R=double(A.torus.R_vox); r=double(A.torus.r_vox);
    c0=double(A.torus.centre_vox);
    % Reprendre la convention EXTERIEUR=SOLIDE du benchmark : les boules
    % peuvent couvrir des centres de voxels virtuels hors de 128^3.
    % Une marge de 32 voxels evite de tronquer ces portions.
    marge=32;
    ref=false(sz128+2*marge);
    ref((1:128)+marge,(1:128)+marge,(1:128)+marge)=ref128;
    sz=size(ref);
    c=c0+marge;

    rep=fullfile(dossierTore,'BENCHMARK_TORE_128_NEUF');
    repTC=fullfile(rep,'TD_TC_step1_sequentiel_audit_k1');
    repBoules=fullfile(rep,'boules_tetra_nettoyes');
    sortie=fullfile(dossierTore,'ANALYSES_GEOMETRIQUES_IMPROVE');
    if ~isfolder(sortie), mkdir(sortie); end
    methodes={'I6','I18','I26','E6','E18','E26','F'};

    % Auto-test cubique sur un tore plus petit et une boule :
    autotest_topologie();

    % Domaine analytique echantillonne avec densite d'aire constante.
    ptsAna=echantillons_surface_tore(R,r,c,240,120);
    bRef=topologie_cubique(ref);
    assert(bRef.Beta0_26==1 && bRef.Beta1==1 && bRef.Beta2==0, ...
        'Le tore de reference ne presente pas la topologie attendue (1,1,0).');
    fprintf('TORE ORIGINAL : beta0=%d, beta1=%d, beta2=%d, chi=%d\n', ...
        bRef.Beta0_26,bRef.Beta1,bRef.Beta2,bRef.Chi);

    % ----------------- 1 et 2 : selections finales ------------------
    enreg=struct([]);
    Uex=cell(1,7);
    fprintf('\n=== TOPOLOGIE ET DISTANCES A LA SURFACE EXACTE ===\n');
    for m=1:numel(methodes)
        nom=methodes{m};
        fichier=fullfile(repTC,sprintf('selection_%s_TD%02d_TC%02d.mat', ...
            nom,tdPct,tcPct));
        assert(isfile(fichier),'Selection introuvable : %s',fichier);
        S=load(fichier,'centresSelectionnes','rayonsSelectionnes');
        C=double(S.centresSelectionnes)+marge;
        RR=double(S.rayonsSelectionnes(:));
        assert(size(C,2)==3 && size(C,1)==numel(RR), ...
            'Donnees centres/rayons incoherentes dans %s',fichier);
        verifier_marge(C,RR,sz,nom);
        U=union_boules(C,RR,sz);
        Uex{m}=U;
        g=geometrie_voxelique(U,ref);
        top=topologie_cubique(U);
        p=points_centres_faces(U);
        dVers=abs(sdf_tore(p,R,r,c));
        dRetour=distances_surface_vers_faces(ptsAna,p,sz);
        ligneFinale=struct( ...
            'Methode',string(nom),'TD_pct',tdPct,'TC_pct',tcPct, ...
            'NBoules',numel(RR),'IoU_pct',g.IoU_pct, ...
            'Couverture_pct',g.Couverture_pct,'Debordement_pct',g.Debordement_pct, ...
            'Beta0_6',top.Beta0_6,'Beta0_26',top.Beta0_26, ...
            'Beta1',top.Beta1,'Beta2',top.Beta2,'Euler',top.Chi, ...
            'ToreTopologiePreservee',logical(top.Beta0_26==1 && top.Beta1==1 && top.Beta2==0), ...
            'D_faces_vers_tore_moy_vox',mean(dVers), ...
            'D_faces_vers_tore_P95_vox',percentile_local(dVers,95), ...
            'D_faces_vers_tore_RMSE_vox',sqrt(mean(dVers.^2)), ...
            'D_tore_vers_faces_moy_APPROX_vox',mean(dRetour), ...
            'D_tore_vers_faces_P95_APPROX_vox',percentile_local(dRetour,95), ...
            'D_sym_moy_APPROX_vox',(mean(dVers)+mean(dRetour))/2, ...
            'D_sym_P95_max_APPROX_vox',max(percentile_local(dVers,95), ...
                                             percentile_local(dRetour,95)) );
        % IMPORTANT : struct([]) est sans champs ; MATLAB refuse d'y
        % affecter directement un struct ayant des champs (erreur ligne 93).
        % On initialise au premier enregistrement avant d'agrandir le tableau.
        if isempty(enreg)
            enreg=ligneFinale;
        else
            enreg(end+1)=ligneFinale; %#ok<AGROW>
        end
        fprintf('%-4s : %5d boules | IoU=%7.3f%% | Betti=(%d,%d,%d) | dA->T=%6.3f | dT->A~=%6.3f vox\n', ...
            nom,numel(RR),g.IoU_pct,top.Beta0_26,top.Beta1,top.Beta2, ...
            mean(dVers),mean(dRetour));
    end
    Tfinal=struct2table(enreg);
    writetable(Tfinal,fullfile(sortie,'01_topologie_distances_7_frontieres.csv'));
    save(fullfile(sortie,'01_topologie_distances_7_frontieres.mat'), ...
         'Tfinal','bRef','-v7.3');

    % ---------------- 3 : budgets strictement egaux -----------------
    fprintf('\n=== BUDGETS IDENTIQUES APRES TD / AVANT TC ===\n');
    lignes=struct([]);
    for m=1:numel(methodes)
        nom=methodes{m};
        fB=fullfile(repBoules,['boules_' nom '.mat']);
        fAudit=fullfile(repTC,['audit_TD_' nom '.mat']);
        assert(isfile(fB) && isfile(fAudit), ...
            'Banque ou audit TD manquant pour %s',nom);
        B=load(fB,'centres','rayons','classeCentres');
        D=load(fAudit,'TD','idsBruts');
        identifiants=find(B.classeCentres(:)==1);
        assert(isequal(double(identifiants(:)),double(D.idsBruts(:))), ...
            'Les identifiants de banque et d''audit TD divergent pour %s',nom);
        C=double(B.centres(identifiants,:))+marge;
        RR=double(B.rayons(identifiants)); RR=RR(:);
        td=double(D.TD(:));
        assert(numel(td)==numel(RR),'Tailles TD/rayons incompatibles pour %s',nom);
        bons=find(isfinite(td) & td<=tdPct/100+1e-12);
        [~,o]=sortrows([-RR(bons),C(bons,:)], [1 2 3 4]);
        bons=bons(o);
        verifier_marge(C(bons,:),RR(bons),sz,[nom ' TD']);

        UB=false(sz);
        nChoisi=0; kSuivant=1;
        for j=1:numel(bons)
            b=bons(j);
            ind=indices_sphere(C(b,:),RR(b),sz);
            if isempty(ind),continue;end
            % Critere commun : au moins un NOUVEAU voxel de la forme d'origine.
            % PAS de TC ici : c'est une experience additionnelle a budget fixe.
            if ~any(ref(ind) & ~UB(ind)),continue;end
            UB(ind)=true;
            nChoisi=nChoisi+1;
            if nChoisi==budgets(kSuivant)
                g=geometrie_voxelique(UB,ref);
                ligneBudget=struct('Methode',string(nom),'TD_pct',tdPct, ...
                    'Budget',budgets(kSuivant),'NBoules',nChoisi, ...
                    'IoU_pct',g.IoU_pct,'Couverture_pct',g.Couverture_pct, ...
                    'Debordement_pct',g.Debordement_pct, ...
                    'BudgetAtteint',true);
                if isempty(lignes)
                    lignes=ligneBudget;
                else
                    lignes(end+1)=ligneBudget; %#ok<AGROW>
                end
                fprintf('%-4s N=%4d : IoU=%7.3f%% | couv=%7.3f%%\n', ...
                    nom,nChoisi,g.IoU_pct,g.Couverture_pct);
                kSuivant=kSuivant+1;
                if kSuivant>numel(budgets),break;end
            end
        end
        % Les budgets impossibles sont NOTEs (pas imputes a leur maximum).
        while kSuivant<=numel(budgets)
            ligneBudget=struct('Methode',string(nom),'TD_pct',tdPct, ...
                'Budget',budgets(kSuivant),'NBoules',nChoisi, ...
                'IoU_pct',NaN,'Couverture_pct',NaN, ...
                'Debordement_pct',NaN,'BudgetAtteint',false);
            if isempty(lignes)
                lignes=ligneBudget;
            else
                lignes(end+1)=ligneBudget; %#ok<AGROW>
            end
            warning('%s : budget %d non atteint (%d boules contributives).', ...
                nom,budgets(kSuivant),nChoisi);
            kSuivant=kSuivant+1;
        end
    end
    Tbudget=struct2table(lignes);
    writetable(Tbudget,fullfile(sortie,'02_comparaison_budget_identique.csv'));
    save(fullfile(sortie,'02_comparaison_budget_identique.mat'),'Tbudget','-v7.3');

    % Figures destinees a l'article.
    % Regressions : l'IoU doit correspondre au benchmark k=1.
    cheminStat=fullfile(repTC,'statistiques_TC.csv');
    if isfile(cheminStat)
        Sref=readtable(cheminStat);
        for m=1:height(Tfinal)
            q=strcmp(string(Sref.Methode),Tfinal.Methode(m)) & ...
              Sref.SeuilTD_pct==tdPct & Sref.SeuilTC_pct==tcPct;
            if nnz(q)==1
                delta=abs(Sref.IoU_pct(q)-Tfinal.IoU_pct(m));
                assert(delta<0.02, ...
                    'IoU non reproduit pour %s (ecart %.4f point).', ...
                    Tfinal.Methode(m),delta);
            end
        end
        fprintf('Controle de regression IoU : OK (ecart < 0.02 point).\n');
    end
    creer_figures(Tfinal,Tbudget,Uex,ref,sortie,marge,tdPct,tcPct);
    fprintf('\nFIN. Tableaux et figures dans :\n%s\n',sortie);
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

function g=geometrie_voxelique(U,ref)
    tp=nnz(U & ref); fp=nnz(U & ~ref);
    g.IoU_pct=100*tp/max(1,nnz(ref)+fp);
    g.Couverture_pct=100*tp/max(1,nnz(ref));
    g.Debordement_pct=100*fp/max(1,nnz(U));
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

function pts=points_centres_faces(U)
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
    pts=[px;py;pz];
end

function d=sdf_tore(P,R,r,c)
    if isempty(P),d=zeros(0,1);return;end
    rho=hypot(P(:,1)-c,P(:,2)-c);
    d=hypot(rho-R,P(:,3)-c)-r;
end

function pts=echantillons_surface_tore(R,r,c,nu,nv)
    % Parametrage torique uniforme en aire (et non uniforme en angle v).
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

function d=distances_surface_vers_faces(P,faces,sz)
    if isempty(faces),d=inf(size(P,1),1);return;end
    % Retour approxime au CENTRE de face le plus proche, pas a la face carree.
    % Si knnsearch est disponible, voisins exacts des centres de faces.
    if exist('knnsearch','file')==2
        [~,d]=knnsearch(faces,P,'K',1);
    else
        % Repli sans Statistics Toolbox : EDT sur grille a demi-voxel.
        dim=2*double(sz)+2;  % 258^3 pour 128^3.
        B=false(dim);
        pos=round(2*faces)+1;
        assert(all(pos(:)>=1) && all(pos(:,1)<=dim(1)) && ...
               all(pos(:,2)<=dim(2)) && all(pos(:,3)<=dim(3)));
        B(sub2ind(dim,pos(:,1),pos(:,2),pos(:,3)))=true;
        D=0.5*bwdist(B);
        p=2*P+1;
        d=double(interpn(D,p(:,1),p(:,2),p(:,3),'linear',NaN));
        assert(all(isfinite(d)),'Echantillonnage de distance hors grille.');
        clear D B
    end
end

function p=percentile_local(x,pc)
    x=sort(double(x(isfinite(x))));
    if isempty(x),p=NaN;return;end
    i=1+(numel(x)-1)*pc/100;
    a=floor(i);b=ceil(i);
    p=x(a)+(i-a)*(x(b)-x(a));
end

function creer_figures(T,Tbudget,U,ref,sortie,marge,tdPct,tcPct)
    f=figure('Visible','off','Color','w','Position',[100 100 1050 460]);
    subplot(1,2,1)
    bar(categorical(cellstr(T.Methode)),T.IoU_pct);
    ylabel('IoU (%)');ylim([0 100]);grid on;
    title(sprintf('Fidelite geometrique - TD%d / TC%d',tdPct,tcPct));
    subplot(1,2,2)
    bar(categorical(cellstr(T.Methode)),T.NBoules);
    ylabel('Nombre de boules');grid on;
    title('Compacite des reconstructions');
    exporter_figure(f,fullfile(sortie,'figure_01_iou_et_nb_boules.png'));close(f);

    f=figure('Visible','off','Color','w','Position',[100 100 850 510]);hold on;
    noms=unique(Tbudget.Methode,'stable');
    for i=1:numel(noms)
        x=Tbudget(Tbudget.Methode==noms(i) & Tbudget.BudgetAtteint,:);
        plot(x.Budget,x.IoU_pct,'-o','LineWidth',1.5,'DisplayName',char(noms(i)));
    end
    xlabel('Budget de spheres identique (N)');ylabel('IoU (%)');
    title('Comparaison a budget fixe - memes banques TD, sans TC');
    grid on;legend('Location','best');
    exporter_figure(f,fullfile(sortie,'figure_02_iou_budget_fixe.png'));close(f);

    f=figure('Visible','off','Color','w','Position',[100 100 1000 390]);
    a={ref,U{7},U{6}}; titres={'Tore original','F','E26'};
    for i=1:3
        subplot(1,3,i);
        imagesc(a{i}(:,:,64+marge)); axis image off;
        colormap(gray(2));title(titres{i});
    end
    sgtitle(sprintf('Coupe XY - original / F / E26 (TD%d TC%d)',tdPct,tcPct));
    exporter_figure(f,fullfile(sortie,'figure_03_comparaison_coupes.png'));close(f);
end

function exporter_figure(h,path)
    try
        exportgraphics(h,path,'Resolution',250);
    catch
        saveas(h,path);
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
