function T = calculer_MSD_RSAE_7_frontieres(ref, dossierBenchmark, nomDataset, modeMSD, tore, configurations)
% CALCULER_MSD_RSAE_7_FRONTIERES
% Calcule 2 mesures uniquement :
%   - MSD : mean surface distance, en voxels, reconstruction -> reference
%   - RSAE : erreur relative absolue de l'aire voxelisee, en %
%
% ENTREE : ref = logical 3D (true = objet / pore) ; tableau deja converti.
% Lit les selections du benchmark existant (pas de triangulation Delaunay).
%   .../TD_TC_step1_sequentiel_audit_k1/selection_I6_TD10_TC10.mat
%   variables centresSelectionnes, rayonsSelectionnes
%
% modeMSD = 'voxel' : distance EUCLIDIENNE depuis chaque voxel de frontiere
%                       reconstruite au plus proche voxel de frontiere
%                       originale (bwperim 6, puis bwdist 3D).
% modeMSD = 'analytique_tore': depuis les CENTRES DE FACES reconstruites
%                       vers la surface torique analytique (exacte).
%   NB : les MSD 'voxel' et 'analytique_tore' utilisent une definition de
%   reference differente : ne pas comparer leurs valeurs directement.
%   tore.R, tore.r, tore.centre sont requis en mode analytique.
%
% RSAE : aire = nombre de faces voxelaires exposees (surface de cubes
% unitaires), dans TOUS les cas ; reference = volume original voxelise.
% RSAE_pct = 100*abs(Arec-Aref)/Aref.
% Hors-image : fond virtuel, non coupe par le volume de reference.
%
% Par defaut : configurations [0 0; 5 5; 10 10].
% Pour une seule configuration, passer [10 10].
% L'absence de spheres donne NaN et statut 'vide' (pas MSD=0).
%
% Necessite Image Processing Toolbox : bwperim, bwdist (mode 'voxel').

    if nargin<5 || isempty(tore), tore=struct(); end
    if nargin<6 || isempty(configurations)
        configurations=[0 0;5 5;10 10];
    end
    assert(islogical(ref) && ndims(ref)==3 && nnz(ref)>0 && nnz(ref)<numel(ref), ...
        'ref doit etre un volume logique 3D non trivial, true=objet.');
    assert(ischar(modeMSD) || isstring(modeMSD),'modeMSD incorrect.');
    modeMSD=char(modeMSD);
    assert(ismember(modeMSD,{'voxel','analytique_tore'}), ...
        'Mode MSD : voxel ou analytique_tore.');
    assert(size(configurations,2)==2 && all(ismember(configurations(:),[0 5 10])), ...
        'Configurations : couples de pourcentages [TD TC] parmi {0,5,10}.');
    if strcmp(modeMSD,'analytique_tore')
        assert(isfield(tore,'R') && isfield(tore,'r') && isfield(tore,'centre'), ...
            'Le tore necessite R, r, centre dans un struct.');
        assert(isscalar(tore.R) && isscalar(tore.r) && tore.R>tore.r && tore.r>0);
        centre=double(tore.centre(:).');
        if numel(centre)==1,centre=repmat(centre,1,3);end
        assert(numel(centre)==3,'Le centre doit etre un scalaire ou [x y z].');
    end
    repSel=fullfile(dossierBenchmark,'TD_TC_step1_sequentiel_audit_k1');
    assert(isfolder(repSel),'Dossier de selections absent : %s',repSel);
    methodes={'I6','I18','I26','E6','E18','E26','F'};
    sz0=double(size(ref));

    % Calcul de l'extension necessaire pour eviter de tronquer les boules.
    padBas=[0 0 0];padHaut=[0 0 0];
    fprintf('\n[%s] Analyse des limites des sphères selectionnees...\n',nomDataset);
    for t=1:size(configurations,1)
        td=configurations(t,1); tc=configurations(t,2);
        for k=1:numel(methodes)
            f=selection_path(repSel,methodes{k},td,tc);
            assert(isfile(f),'Selection absente : %s',f);
            S=load(f,'centresSelectionnes','rayonsSelectionnes');
            assert(isfield(S,'centresSelectionnes') && isfield(S,'rayonsSelectionnes'), ...
                'Variables centresSelectionnes / rayonsSelectionnes absentes de %s',f);
            R=double(S.rayonsSelectionnes(:));
            C=double(S.centresSelectionnes);
            assert(size(C,2)==3 && size(C,1)==numel(R), ...
                'Centres / rayons incoherents : %s',f);
            assert(all(isfinite(C(:))) && all(isfinite(R)) && all(R>=0), ...
                'Valeurs de sphere non finies ou rayons negatifs : %s',f);
            if isempty(R),continue;end
            delta=1e-12*max(1,R);
            mini=min(ceil(C-(R+delta)),[],1);
            maxi=max(floor(C+(R+delta)),[],1);
            padBas=max(padBas,max(0,1-mini));
            padHaut=max(padHaut,max(0,maxi-sz0));
        end
    end
    padBas=ceil(padBas); padHaut=ceil(padHaut);
    sz=sz0+padBas+padHaut;
    fprintf('Padding bas=%s, haut=%s, grille=%s\n', ...
        mat2str(padBas),mat2str(padHaut),mat2str(sz));
    if prod(sz)>1.8*prod(sz0)
        warning('Le padding multiplie le volume par %.2f : RAM importante.',prod(sz)/prod(sz0));
    end
    refPad=false(sz);
    ii=(1:sz0(1))+padBas(1);
    jj=(1:sz0(2))+padBas(2);
    kk=(1:sz0(3))+padBas(3);
    refPad(ii,jj,kk)=ref;
    clear ref

    Aref=aire_faces_voxels(refPad);
    fprintf('Aire reference voxelisee = %.0f faces unite\n',Aref);

    if strcmp(modeMSD,'voxel')
        fprintf('Calcul EDT 3D depuis voxels de frontiere originaux (RAM importante pour 512^3)...\n');
        bRef=bwperim(refPad,6);
        Dref=bwdist(bRef);  % Distances euclidiennes en voxel-centres.
        clear bRef
    else
        Dref=[];
        c=tore.centre(:).';
        if numel(c)==1,c=repmat(c,1,3);end
        c=double(c)+padBas;
        % Baseline simple : MSD de la frontiere voxelisee originale au tore.
        base=msd_analytique_tore(refPad,double(tore.R),double(tore.r),c);
        fprintf('MSD tore voxelise original vers analytique = %.6f voxels\n',base);
    end

    statFile=fullfile(repSel,'statistiques_TC.csv');
    stat=table();
    if isfile(statFile)
        stat=readtable(statFile);
    else
        warning('statistiques_TC.csv absent : verification des IoU impossible.');
    end

    N=numel(methodes)*size(configurations,1);
    Dataset=repmat(string(nomDataset),N,1);
    Boundary=strings(N,1);
    TD_pct=zeros(N,1);TC_pct=zeros(N,1);
    NBoules=zeros(N,1);
    MSD_vox=nan(N,1);RSAE_pct=nan(N,1);
    AireRef_faces=repmat(Aref,N,1);
    AireRec_faces=nan(N,1);
    IoU_controle_pct=nan(N,1);
    Statut=strings(N,1);
    kL=0;
    for t=1:size(configurations,1)
        td=configurations(t,1);tc=configurations(t,2);
        fprintf('\n====== TD%02d TC%02d / %s ======\n',td,tc,nomDataset);
        for k=1:numel(methodes)
            nom=methodes{k};kL=kL+1;
            Boundary(kL)=string(nom);
            TD_pct(kL)=td;TC_pct(kL)=tc;
            f=selection_path(repSel,nom,td,tc);
            S=load(f,'centresSelectionnes','rayonsSelectionnes');
            C=double(S.centresSelectionnes)+padBas;
            R=double(S.rayonsSelectionnes(:));
            NBoules(kL)=numel(R);
            if isempty(R)
                Statut(kL)="vide";
                fprintf('%-4s : pas de spheres, MSD et RSAE non definies.\n',nom);
                continue;
            end
            ticM=tic;
            U=union_spheres(C,R,sz);
            Au=aire_faces_voxels(U);
            AireRec_faces(kL)=Au;
            RSAE_pct(kL)=100*abs(Au-Aref)/Aref;
            if strcmp(modeMSD,'voxel')
                bRec=bwperim(U,6);
                assert(any(bRec(:)),'Reconstruction sans frontiere non vide : %s',nom);
                % MSD : moyenne des distances au plus proche voxel de
                % frontiere de l'image originale (centre de voxel).
                MSD_vox(kL)=mean(double(Dref(bRec)));
                clear bRec
            else
                % MSD vers la surface ANALYTIQUE exacte du tore.
                % Meme definition que dans analyser_tore_128_3_points.m.
                MSD_vox(kL)=msd_analytique_tore(U,double(tore.R),double(tore.r),c);
            end
            % IoU de controle : pas une 3eme metrique de surface.
            tp=nnz(U & refPad);fp=nnz(U & ~refPad);
            IoU_controle_pct(kL)=100*tp/(nnz(refPad)+fp);
            if ~isempty(stat) && all(ismember( ...
                    {'Methode','SeuilTD_pct','SeuilTC_pct','IoU_pct'},stat.Properties.VariableNames))
                q=strcmp(string(stat.Methode),string(nom)) & ...
                    abs(stat.SeuilTD_pct-td)<1e-8 & abs(stat.SeuilTC_pct-tc)<1e-8;
                if nnz(q)==1
                    deltaIoU=abs(IoU_controle_pct(kL)-stat.IoU_pct(q));
                    assert(deltaIoU<0.05, ...
                        ['La reconstruction %s TD%d TC%d ne reproduit pas ', ...
                         'l IoU du benchmark (ecart %.3f point).'],nom,td,tc,deltaIoU);
                end
            end
            Statut(kL)="OK";
            fprintf('%-4s : %7d spheres | MSD=%8.4f vox | RSAE=%8.3f %% | ', ...
                nom,numel(R),MSD_vox(kL),RSAE_pct(kL));
            fprintf('IoU controle=%.3f %% | %.1f s\n',IoU_controle_pct(kL),toc(ticM));
            clear U S C R
        end
    end
    T=table(Dataset,Boundary,TD_pct,TC_pct,NBoules, ...
        MSD_vox,RSAE_pct,AireRef_faces,AireRec_faces,IoU_controle_pct,Statut);
    sortie=fullfile(dossierBenchmark,'SURFACE_MSD_RSAE');
    if ~isfolder(sortie),mkdir(sortie);end
    pCsv=fullfile(sortie,sprintf('MSD_RSAE_%s.csv',nomDataset));
    pMat=fullfile(sortie,sprintf('MSD_RSAE_%s.mat',nomDataset));
    writetable(T,pCsv);
    meta=struct('Dataset',nomDataset,'ModeMSD',modeMSD, ...
        'ReferenceMSD',reference_description(modeMSD), ...
        'RSAE_reference','Aire des faces exposees des voxels de l image originale', ...
        'PaddingBas',padBas,'PaddingHaut',padHaut, ...
        'Configurations',configurations);
    save(pMat,'T','meta','-v7.3');
    fprintf('\nRESULTATS :\n%s\n%s\n',pCsv,pMat);
    disp(T(:,{'Boundary','TD_pct','TC_pct','NBoules','MSD_vox','RSAE_pct','Statut'}));
end

function p=selection_path(repSel,nom,td,tc)
    p=fullfile(repSel,sprintf('selection_%s_TD%02d_TC%02d.mat',nom,td,tc));
end

function A=aire_faces_voxels(U)
    % Aire totale des faces exposees de cubes unitaires, EXTERIEUR inclus.
    % 6 faces par voxel moins 2 pour chaque paire face-adjacente.
    n=double(nnz(U));
    intFaces=double(nnz(U(1:end-1,:,:) & U(2:end,:,:))) + ...
        double(nnz(U(:,1:end-1,:) & U(:,2:end,:))) + ...
        double(nnz(U(:,:,1:end-1) & U(:,:,2:end)));
    A=6*n-2*intFaces;
end

function U=union_spheres(C,R,sz)
    U=false(sz);
    for b=1:numel(R)
        ind=indices_sphere(C(b,:),R(b),sz);
        U(ind)=true;
        if mod(b,50000)==0
            fprintf('    union : %d / %d spheres\n',b,numel(R));
        end
    end
end

function ind=indices_sphere(c,r,sz)
    tol=1e-12*max(1,r*r);
    dl=1e-12*max(1,r);
    lo=max([1 1 1],ceil(c-r-dl));
    hi=min(sz,floor(c+r+dl));
    if any(lo>hi),ind=zeros(0,1);return;end
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

function msd=msd_analytique_tore(U,R,r,c)
    % Calcul face par face, sans construire un gros nuage 3D en memoire.
    a=size(U,1);b=size(U,2);d=size(U,3);
    somme=0;nFaces=0;
    ax=cat(1,false(1,b,d),U);bx=cat(1,U,false(1,b,d));
    [i,j,k]=ind2sub(size(ax),find(xor(ax,bx)));
    if ~isempty(i)
        p=[double(i)-0.5,double(j),double(k)];
        di=distance_tore(p,R,r,c);somme=somme+sum(di);nFaces=nFaces+numel(di);
    end
    clear ax bx i j k p di
    ay=cat(2,false(a,1,d),U);by=cat(2,U,false(a,1,d));
    [i,j,k]=ind2sub(size(ay),find(xor(ay,by)));
    if ~isempty(i)
        p=[double(i),double(j)-0.5,double(k)];
        di=distance_tore(p,R,r,c);somme=somme+sum(di);nFaces=nFaces+numel(di);
    end
    clear ay by i j k p di
    az=cat(3,false(a,b,1),U);bz=cat(3,U,false(a,b,1));
    [i,j,k]=ind2sub(size(az),find(xor(az,bz)));
    if ~isempty(i)
        p=[double(i),double(j),double(k)-0.5];
        di=distance_tore(p,R,r,c);somme=somme+sum(di);nFaces=nFaces+numel(di);
    end
    msd=somme/nFaces;
end

function d=distance_tore(p,R,r,c)
    rho=hypot(p(:,1)-c(1),p(:,2)-c(2));
    d=abs(hypot(rho-R,p(:,3)-c(3))-r);
end

function s=reference_description(modeMSD)
    if strcmp(modeMSD,'analytique_tore')
        s='Surface analytique continue du tore ; centres des faces reconstruites';
    else
        s='Voxels de frontiere 6-adjacents de l image originale ; voxels de frontiere reconstruite';
    end
end
