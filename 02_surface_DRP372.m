%% 02_surface_DRP372.m - MSD et RSAE, DRP-372 (256^3)
% Lancer independamment ; le fichier original code 0=pore, 1=solide.

fichier='374_02_024_256.mat';
benchmark='BENCHMARK_DRP372_256_NEUF';
configs=[0 0;5 5;10 10]; % Pour TD10/TC10 uniquement : [10 10]

S=load(fichier,'bin');
assert(isfield(S,'bin'),'La variable bin est absente');
assert(isequal(size(S.bin),[256 256 256]),'DRP-372 : dimensions 256^3 attendues');
ref=(S.bin==0);
clear S
fprintf('DRP-372 : pores=%d (attendu 688309)\n',nnz(ref));
assert(nnz(ref)==688309, ...
    'Erreur de convention ou de fichier : %d pores (attendu 688309)',nnz(ref));

T_DRP372=calculer_MSD_RSAE_7_frontieres( ...
    ref,benchmark,'DRP372','voxel',struct(),configs);
