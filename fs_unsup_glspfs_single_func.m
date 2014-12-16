function [FeaNumCandi,res_gs,res_aio, res_gs_ps] = fs_unsup_glspfs_single_func(dataset, exp_settings, algo_settings)
%Unsupervised feature selection using GLSPFS

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

%===================setup=======================
local_type_candi = {'LPP', 'LLE', 'LTSA'};
local_type_param_candi = {[], [], []};
knn_size_candi = 5;
lambda1_candi = 10.^[-5:5];
lambda2_candi = 10.^[-5:5];
s1 = optSigma(X);
global_kernel_cell_candi = buildParamKernel({'Gaussian'}, {sqrt(2.^[-1]) * s1}, {''});
local_type_param_candi{1} = [sqrt(2.^[-1]) * s1];
local_type_param_candi{3} = [nClass];
paramCell = fs_unsup_glspfs_build_param(local_type_candi, local_type_param_candi, knn_size_candi, ...
    lambda1_candi, lambda2_candi, global_kernel_cell_candi);
%===============================================

disp('GLSPFS ...');
t_start = clock;
feaSubsets = cell(length(paramCell), 1);
for i1 = 1:length(paramCell)
    fprintf('GLSPFS parameter search %d out of %d...\n', i1, length(paramCell));
    param = paramCell{i1};
    K = constructKernel(X, X, param.global_kernel_option);
    L = computeLocalStructure(X, param.local_type, param.local_k, param.local_lpp_sigma, param.local_ltsa_embedded_dim);
    feaSubsets{i1,1} = fs_unsup_glspfs(X, K, L, param.lambda1, param.lambda2, max(FeaNumCandi));
end
t_end = clock;
t1 = etime(t_end,t_start);
disp(['exe time: ',num2str(t1)]);

t_start = clock;
disp('evaluation ...');
res_aio = cell(length(paramCell), length(FeaNumCandi));
for i2 = 1:length(FeaNumCandi)
    for i1 = 1:length(paramCell)
        fprintf('GLSPFS parameter evaluation %d outof %d  ... %d out of %d...\n', i2, length(FeaNumCandi), i1, length(paramCell));
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

save(fullfile(prefix_mdcs, [dataset, '_best_result_GLSPFS.mat']),'FeaNumCandi','res_gs','res_aio', 'res_gs_ps');
end