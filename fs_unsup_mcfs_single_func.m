function [FeaNumCandi,res_gs,res_aio, res_gs_ps] = fs_unsup_mcfs_single_func(dataset, exp_settings, algo_settings)
%Unsupervised feature selection using MCFS

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

%======================setup===========================
knnCandi = 5;
weightCandi = {'Binary','HeatKernel'};
s1 = optSigma(X);
weight_param_Candi = {[], 2.^[-3:3] .* s1.^2};
paramCell = fs_unsup_mcfs_build_param(knnCandi, weightCandi, weight_param_Candi);
%======================================================

t_start = clock;
feaSubsets = cell(length(paramCell), 1);
valid_ids = zeros(length(paramCell), 1);
parfor i1 = 1:length(paramCell)
    fprintf(['MCFS parameter search %d out of %d...\n'], i1, length(paramCell));
    param = paramCell{i1};
    W = constructW(X, param);
    options = [];
    options.nUseEigenfunction = nClass;
    options.W = W;
    % some may failed  due to SR code
    try
        index = fs_unsup_mcfs(X,max(FeaNumCandi),options);
        feaSubsets{i1} = index{1};
    catch
        valid_ids(i1) = 1;
    end
end
t_end = clock;
t1 = etime(t_end,t_start);
disp(['exe time: ',num2str(t1)]);

t_start = clock;
disp('evaluation ...');
valid_ids = find(valid_ids == 0);
paramCell_old = paramCell;
feaSubsets_old = feaSubsets;
paramCell = cell(length(valid_ids), 1);
feaSubsets = cell(length(valid_ids), 1);
for i1=1:length(valid_ids)
    paramCell{i1} = paramCell_old{valid_ids(i1)};
    feaSubsets{i1} = feaSubsets_old{valid_ids(i1)};
end
res_aio = cell(length(paramCell), length(FeaNumCandi));
for i2 = 1:length(FeaNumCandi)
    parfor i1 = 1:length(paramCell)
        tmp = feaSubsets{i1, 1};
        fprintf('MCFS parameter evaluation %d outof %d  ... %d out of %d...\n', i2, length(FeaNumCandi), i1, length(paramCell));
        res_aio{i1, i2} = evalUnSupFS(X, Y, tmp(1:FeaNumCandi(i2)), struct('nKm', nKmeans));
    end
end
[res_gs,res_gs_ps] = grid_search_fs(res_aio);
res_gs.feaset = FeaNumCandi;
t_end = clock;
t2 = etime(t_end,t_start);
disp(['exe time: ',num2str(t2)]);
res_gs.time = t1;
res_gs.time2 = t2;

save(fullfile(prefix_mdcs, [dataset, '_best_result_MCFS.mat']),'FeaNumCandi','res_gs','res_aio', 'res_gs_ps');
end