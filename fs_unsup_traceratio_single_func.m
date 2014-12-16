function [FeaNumCandi,res_gs,res_aio, res_gs_ps] = fs_unsup_traceratio_single_func(dataset, exp_settings, algo_settings)
%use trace ratio to select features.

%======================setup===========================
FeaNumCandi = exp_settings.FeaNumCandi;
nKmeans = exp_settings.nKmeans;
prefix_mdcs = [];
if isfield(exp_settings, 'prefix_mdcs')
    prefix_mdcs = exp_settings.prefix_mdcs;
end
%======================================================

disp(['dataset:',dataset]);
[X, Y] = extractXY(dataset);
[nSmp,nDim] = size(X);
nClass = length(unique(Y));

%================setup======================
knnCandi = 5;
n1 = length(knnCandi);
nP = n1;
paramCell = cell(nP, 1);
idx = 0;
for i1 = 1:n1
    param = [];
    param.k = knnCandi(i1);
    idx = idx + 1;
    paramCell{idx} = param;
end
%===========================================

t_start = clock;
disp('trace ratio...');
feaSubsets = cell(length(paramCell), length(FeaNumCandi));
for i1 = 1:length(paramCell)
    fprintf('UDFS parameter search %d out of %d...\n', i1, length(paramCell));
    param = paramCell{i1};
    [Sb, Sw] = compute_SbSw_unsup(X, param.k);
    parfor i2 = 1:length(FeaNumCandi)
        feaSubsets{i1, i2} = fs_unsup_traceratio(Sb, Sw, FeaNumCandi(i2));
    end
end
t_end = clock;
t1 = etime(t_end,t_start);
disp(['exe time: ',num2str(t1)]);

t_start = clock;
disp('evaluation ...');
res_aio = cell(1, length(FeaNumCandi)); 
parfor i1 = 1:length(FeaNumCandi)
    idx = feaSubsets{i1};
    res_aio{1, i1} = evalUnSupFS(X, Y, idx(1:FeaNumCandi(i1)), struct('nKm', nKmeans));
end
[res_gs, res_gs_ps] = grid_search_fs(res_aio);
res_gs.feaset = FeaNumCandi;
t_end = clock;
t2 = etime(t_end,t_start);
disp(['exe time: ',num2str(t2)]);
res_gs.time = t1;
res_gs.time2 = t2;

save(fullfile(prefix_mdcs, [dataset, '_best_result_TraceRatio.mat']),'FeaNumCandi','res_gs','res_aio', 'res_gs_ps');
end