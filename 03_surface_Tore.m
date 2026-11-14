%% 03_surface_Tore.m - MSD analytique et RSAE, tore (128^3)
% MSD : EXACTEMENT comme le script analyser_tore_128_3_points.m :
%   distances depuis centres des faces reconstruites vers tore analytique.
% RSAE : comparaison des faces exposees dans la reconstruction et le tore
%   voxelise original (PAS avec l'aire de la surface analytique 4*pi^2*R*r).

repertoireTore='TORE_SYNTHETIQUE_128';
fichier=fullfile(repertoireTore,'tore_analytique_128.mat');
benchmark=fullfile(repertoireTore,'BENCHMARK_TORE_128_NEUF');
configs=[0 0;5 5;10 10];  % Pour TD10/TC10 uniquement : [10 10]

S=load(fichier,'imgTore','torus');
assert(isfield(S,'imgTore') && isfield(S,'torus'), ...
    'Variables imgTore et torus requises');
ref=logical(S.imgTore);
assert(isequal(size(ref),[128 128 128]),'Tore : dimensions 128^3 attendues');
torus=struct('R',double(S.torus.R_vox), ...
             'r',double(S.torus.r_vox), ...
             'centre',double(S.torus.centre_vox));
clear S
fprintf('Tore : foreground=%d (attendu 182384)\n',nnz(ref));
assert(nnz(ref)==182384,'Tore : foreground inattendu : %d',nnz(ref));

T_Tore=calculer_MSD_RSAE_7_frontieres( ...
    ref,benchmark,'Tore','analytique_tore',torus,configs);
