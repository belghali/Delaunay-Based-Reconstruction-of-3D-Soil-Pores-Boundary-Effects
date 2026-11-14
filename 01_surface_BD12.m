%% 01_surface_BD12.m - MSD et RSAE, BD12 (512^3)
% Doit etre lance independamment des autres scripts.
% Reconstruit les unions des selections TD/TC deja sauvegardees.
% Original = image de pore true. Zero=solide OU 0=pore selon la source.
% Image Processing Toolbox requis. RAM : bwdist 512^3 ~ >1 Go au total.

fichier='bd12.mat';
variable='img';          % Adapter si besoin : whos('-file','bd12.mat')
convention='zero';       % 'zero' pour 0=pore ; 'nonzero' pour valeur non nulle=pore
benchmark='BENCHMARK_512_NEUF';
configs=[0 0;5 5;10 10]; % Pour TD10/TC10 uniquement : [10 10]

S=load(fichier,variable);
assert(isfield(S,variable),'Variable introuvable dans le MAT.');
brut=S.(variable);clear S
assert(isequal(size(brut),[512 512 512]),'BD12 : dimensions 512^3 attendues');
if strcmp(convention,'zero')
    ref=(brut==0);
elseif strcmp(convention,'nonzero')
    ref=(brut~=0);
else
    error('convention : zero ou nonzero');
end
clear brut
fprintf('BD12 : pores=%d (attendu 23056832)\n',nnz(ref));
assert(nnz(ref)==23056832, ...
    'Erreur de convention ou de fichier : %d pores (attendu 23056832)',nnz(ref));

T_BD12=calculer_MSD_RSAE_7_frontieres( ...
    ref,benchmark,'BD12','voxel',struct(),configs);
