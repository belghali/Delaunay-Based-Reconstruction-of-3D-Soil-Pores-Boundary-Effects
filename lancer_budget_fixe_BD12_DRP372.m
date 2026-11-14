% LANCER_BUDGET_FIXE_BD12_DRP372
% Appeler apres avoir place ce script et budget_fixe_7_frontieres.m
% dans le meme dossier MATLAB (ou ajouter ce dossier au path).
%
% ADAPTER LES CHEMINS CI-DESSOUS A VOTRE ORDINATEUR.
% Les benchmarks originaux doivent deja avoir termine TD, et leurs
% sous-dossiers boules_tetra_nettoyes et TD_TC_step1_sequentiel_audit_k1
% doivent etre presents.

% --- BD12 : 512^3, foreground (pore) = 1 ou 255 ---
benchmarkBD12 = 'BENCHMARK_512_NEUF';
fichierBD12   = 'bd12.mat';  % a adapter si le fichier porte un autre nom
budgetsBD12   = [1000 5000 10000 25000 50000];

% --- DRP-372 : 256^3, variable bin (0 = pore, 1 = solide) ---
benchmarkDRP  = 'BENCHMARK_DRP372_256_NEUF';
fichierDRP    = '374_02_024_256.mat';
budgetsDRP    = [100 250 500 1000 2500];

% --- 1/2 : BD12 ---
% Le fichier BD12 peut avoir plusieurs noms de variables selon sa version.
% La fonction ci-dessous selectionne uniquement un tableau binaire 512^3
% contenant exactement 23 056 832 voxels foreground : aucun choix aveugle.
refBD12=charger_BD12(fichierBD12);
T_BD12=budget_fixe_7_frontieres(refBD12,benchmarkBD12, ...
    budgetsBD12,fullfile(benchmarkBD12,'BUDGET_FIXE_TD10'),'BD12');
clear refBD12

% --- 2/2 : DRP-372 ---
assert(isfile(fichierDRP),'Fichier DRP introuvable : %s',fichierDRP);
D=load(fichierDRP,'bin');
assert(isfield(D,'bin'),'Variable bin absente de %s.',fichierDRP);
assert(isequal(size(D.bin),[256 256 256]),'DRP : dimensions inattendues.');
assert(all(D.bin(:)==0 | D.bin(:)==1),'DRP doit contenir uniquement 0/1.');
refDRP=(D.bin==0); % ! IMPORTANT : 0 = PORE POUR DRP-372 !
clear D
assert(nnz(refDRP)==688309, ...
    'DRP : nombre de pores inattendu (%d au lieu de 688309).',nnz(refDRP));
T_DRP=budget_fixe_7_frontieres(refDRP,benchmarkDRP, ...
    budgetsDRP,fullfile(benchmarkDRP,'BUDGET_FIXE_TD10'),'DRP-372');
clear refDRP

% --- Resume commun ---
T_ensemble=[T_BD12;T_DRP];
writetable(T_ensemble,'BUDGET_FIXE_BD12_DRP372_COMBINE.csv');
save('BUDGET_FIXE_BD12_DRP372_COMBINE.mat', ...
     'T_BD12','T_DRP','T_ensemble','-v7.3');
disp(T_ensemble(:,{'Dataset','Methode','Budget','IoU_pct','BudgetAtteint'}));

function ref=charger_BD12(fichier)
    assert(isfile(fichier),'Fichier BD12 introuvable : %s',fichier);
    info=whos('-file',fichier);
    candidat={};
    for i=1:numel(info)
        q=info(i);
        if ~isequal(q.size,[512 512 512]),continue;end
        if ~ismember(q.class,{'logical','uint8','int8','uint16','int16', ...
                'uint32','int32','single','double'}),continue;end
        v=load(fichier,q.name);
        x=v.(q.name);
        if all(x(:)==0 | x(:)==1 | x(:)==255) && nnz(x) == 23056832
            candidat{end+1}=q.name; %#ok<SAGROW>
            ref=logical(x~=0);
        end
        clear v x
    end
    assert(numel(candidat)==1, ...
        ['BD12 : impossible d''identifier de maniere unique la variable ', ...
         '512^3, 23 056 832 voxels poreux dans %s. ', ...
         'Verifiez avec whos(''-file'',''%s'') puis fournissez directement ', ...
         'votre tableau img a budget_fixe_7_frontieres.'],fichier,fichier);
    fprintf('BD12 : variable detectee = %s\n',candidat{1});
end
