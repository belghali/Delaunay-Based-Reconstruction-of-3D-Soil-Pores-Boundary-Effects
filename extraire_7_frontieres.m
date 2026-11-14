function [P, stats] = extraire_7_frontieres(img)
%EXTRAIRE_7_FRONTIERES Sept representations de la frontiere d'un volume 3D.
%
% [P, stats] = extraire_7_frontieres(img)
%
% Entree :
%   img : image binaire 3D, avec 255 (ou 1/true) = PORE/OBJET
%         et 0 = SOLIDE/FOND.
%
% Sorties :
%   P.I6, P.I18, P.I26 : frontieres INTERIEURES (centres de voxels poreux)
%   P.E6, P.E18, P.E26 : frontieres EXTERIEURES (centres de voxels solides)
%   P.F              : centres des faces pore-solide, y compris les faces
%                      de decoupe du volume.
%   Chaque champ est une matrice N x 3 de coordonnees [i,j,k] (double),
%   dans le repere MATLAB de l'image originale (1-based).
%
% Convention de bord : en dehors du tableau, la phase est supposée solide.
% Pour l'implementer, une couche de solide virtuelle entoure temporairement
% l'image. Les points E peuvent donc avoir 0 ou dim+1 comme coordonnee;
% les centres de faces F peuvent avoir 0.5 ou dim+0.5.
%
% stats.PointsLiesAuBord : pour I, points ajoutes par la convention de bord;
% pour E, points situees dans la couche virtuelle; pour F, faces situees
% sur les six limites du volume.
%
% Requiert Image Processing Toolbox (imdilate).

    if ~(islogical(img) || isnumeric(img)) || ~isreal(img) || ...
            isempty(img) || ndims(img) ~= 3
        error('img doit etre une image binaire 3D non vide.');
    end
    if any(img(:) ~= 0 & img(:) ~= 1 & img(:) ~= 255)
        error('Les valeurs permises sont 0 et 255 (ou 0 et 1).');
    end

    S = (img ~= 0);              % true = pore, false = solide
    sz = size(S);

    % Une couche de voxels solides virtuels tout autour de l'image.
    Sp = false(sz + 2);
    Sp(2:end-1, 2:end-1, 2:end-1) = S;

    % Voisinages 6, 18, 26 (centre inclus, sans consequence sur le resultat).
    [di,dj,dk] = ndgrid(-1:1, -1:1, -1:1);
    d1 = abs(di) + abs(dj) + abs(dk);
    V = {d1 <= 1, d1 <= 2, true(3,3,3)};
    suffixes = {'6','18','26'};

    P = struct();
    bordI = zeros(3,1);
    bordE = zeros(3,1);

    for t = 1:3
        voisinage = V{t};

        % Tout point de pore adjacent a un solide (y compris hors-volume).
        MI = Sp & imdilate(~Sp, voisinage);
        % Tout point de solide adjacent a un pore (y compris couche virtuelle).
        ME = ~Sp & imdilate(Sp, voisinage);

        nomI = ['I' suffixes{t}];
        nomE = ['E' suffixes{t}];
        P.(nomI) = masque_vers_points(MI) - 1;
        P.(nomE) = masque_vers_points(ME) - 1;

        % Vrais nouveaux points interieurs dus aux limites du volume.
        MI_sans_bord = S & imdilate(~S, voisinage);
        MI_dans_img = MI(2:end-1, 2:end-1, 2:end-1);
        bordI(t) = nnz(MI_dans_img & ~MI_sans_bord);

        % Points exterieurs situes dans la couche virtuelle.
        Q = P.(nomE);
        hors = false(size(Q,1),1);
        for a = 1:3
            hors = hors | Q(:,a) < 1 | Q(:,a) > sz(a);
        end
        bordE(t) = nnz(hors);
    end

    % Points au centre de chaque face separant deux voxels de phases
    % differentes, y compris les faces entre pore et solide virtuel.
    faces = cell(3,1);
    for a = 1:3
        low = {':',':',':'};
        high = {':',':',':'};
        low{a} = 1:size(Sp,a)-1;
        high{a} = 2:size(Sp,a);
        A = Sp(low{:});
        B = Sp(high{:});
        [i,j,k] = ind2sub(size(A), find(A ~= B));
        Q = double([i,j,k]);
        Q(:,a) = Q(:,a) + 0.5;
        faces{a} = Q - 1;  % coordonnees dans l'image non rembourree
    end
    P.F = vertcat(faces{:});

    % Faces placees sur les limites physiques du domaine numerique.
    Q = P.F;
    sur_limite = false(size(Q,1),1);
    for a = 1:3
        sur_limite = sur_limite | Q(:,a) == 0.5 | Q(:,a) == sz(a)+0.5;
    end
    bordF = nnz(sur_limite);

    methodes = {'I6';'I18';'I26';'E6';'E18';'E26';'F'};
    NombrePoints = zeros(7,1);
    for t = 1:7
        NombrePoints(t) = size(P.(methodes{t}),1);
    end
    PointsLiesAuBord = [bordI; bordE; bordF];
    stats = table(methodes, NombrePoints, PointsLiesAuBord, ...
        'VariableNames', {'Methode','NombrePoints','PointsLiesAuBord'});

    fprintf('\nDimensions du volume : %s\n', mat2str(sz));
    fprintf('Voxels poreux        : %d\n', nnz(S));
    fprintf('Convention de bord   : exterieur = solide\n');
    disp(stats);
end

function Q = masque_vers_points(M)
    [i,j,k] = ind2sub(size(M), find(M));
    Q = double([i,j,k]);
end
