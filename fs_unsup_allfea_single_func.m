function [FeaNumCandi,res_gs,res_aio, res_gs_ps] = fs_unsup_allfea_single_func(dataset, exp_settings, algo_settings)
%Unsupervised feature selection using AllFea

%======================setup===========================
FeaNumCandi = exp_settings.FeaNumCandi;
nKmeans = exp_settings.nKmeans;
prefix_mdcs = [];
if isfield(exp_settings, 'prefix_mdcs')
    prefix_mdcs = exp_settings.prefix_mdcs;
end
%===============================================
[X, Y] = extractXY(dataset);
[nSmp,nDim] = size(X);

t_start = clock;
disp('get AllFea ...');
fs_res = evalUnSupFS(X, Y, [1:nDim], struct('nKm', nKmeans));
res_aio = cell(1, length(FeaNumCandi)); 
parfor feaIdx = 1:length(FeaNumCandi)
    res_aio{1, feaIdx} = fs_res;
end
[res_gs, res_gs_ps] = grid_search_fs(res_aio);
res_gs.feaset = FeaNumCandi;
t_end = clock;
t1 = etime(t_end,t_start);
disp(['exe time: ',num2str(t1)]);
res_gs.time = t1;
res_gs.time2 = t1;

save(fullfile(prefix_mdcs, [dataset, '_best_result_AllFea.mat']),'FeaNumCandi','res_gs','res_aio', 'res_gs_ps');
end