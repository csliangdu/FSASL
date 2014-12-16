function [FeaNumCandi,res_gs,res_aio, res_gs_ps] = fs_unsup_lapscore_single_func(dataset, exp_settings, algo_settings)
%Unsupervised feature selection using LapScore

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
weightCandi = {'HeatKernel'};%{'Binary','HeatKernel'};
s1 = optSigma(X);
weight_param_Candi = {2.^[-3:3] .* s1.^2};% {[], 2.^[-3:3] .* s1.^2};
paramCell = fs_unsup_lapscore_build_param(knnCandi, weightCandi, weight_param_Candi);
%===========================================


disp('LapScore ...');
t_start = clock;
feaSubsets = cell(length(paramCell), 1);
parfor i1 = 1:length(paramCell)
    fprintf(['LapScore parameter search %d out of %d...\n'], i1, length(paramCell));
    param = paramCell{i1};
    W = constructW(X, param);
    LS = fs_unsup_lapscore(X, W);
    [~, idx] = sort(-LS);
    feaSubsets{i1,1} = idx;
end
t_end = clock;
t1 = etime(t_end,t_start);
disp(['exe time: ',num2str(t1)]);

disp('evaluation....');
t_start = clock;
res_aio = cell(length(paramCell), length(FeaNumCandi));
for i2 = 1:length(FeaNumCandi)
    parfor i1 = 1:length(paramCell)
        fprintf('LapScore parameter evaluation %d outof %d  ... %d out of %d...\n', i2, length(FeaNumCandi), i1, length(paramCell));
        idx = feaSubsets{i1,1};    
        res_aio{i1, i2} = evalUnSupFS(X, Y, idx(1:FeaNumCandi(i2)), struct('nKm', nKmeans));
    end
end
[res_gs, res_gs_ps] = grid_search_fs(res_aio);
res_gs.feaset = FeaNumCandi;
t_end = clock;
t2 = etime(t_end,t_start);
disp(['exe time: ',num2str(t2)]);
res_gs.time = t1;
res_gs.time2 = t2;

save(fullfile(prefix_mdcs, [dataset, '_best_result_LapScore.mat']), 'FeaNumCandi','res_gs','res_aio', 'res_gs_ps');
end