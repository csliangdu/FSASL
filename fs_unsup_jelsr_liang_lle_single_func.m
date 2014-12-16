function [FeaNumCandi,res_gs,res_aio, res_gs_ps] = fs_unsup_jelsr_liang_lle_single_func(dataset, exp_settings, algo_settings)
%Unsupervised feature selection using JELSR_liang

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
r1Candi = 10.^[-5:5];
r2Candi = 10.^[-5:5];
knnCandi = 5;
weightCandi = {'lle'};
s1 = optSigma(X);
weight_param_Candi = {s1};
paramCell = fs_unsup_jelsr_build_param(knnCandi, weightCandi, weight_param_Candi, r1Candi, r2Candi);
%===============================================
disp('JELSR ...');
t_start = clock;
feaSubsets = cell(length(paramCell), 1);
parfor i1 = 1:length(paramCell)
    fprintf('JELSR_liang parameter search %d out of %d...\n', i1, length(paramCell));
    param = paramCell{i1};
    param.nClusters = nClass;
    [model_jelsr] = fs_unsup_jelsr_liang(X', param);
    [~,idx] = sort(model_jelsr.z, 'descend');
    feaSubsets{i1,1} = idx;
end
t_end = clock;
t1 = etime(t_end,t_start);
disp(['exe time: ',num2str(t1)]);

t_start = clock;
disp('evaluation....');
res_aio = cell(length(paramCell), length(FeaNumCandi));
for i2 = 1:length(FeaNumCandi)
    parfor i1 = 1:length(paramCell)
        fprintf('JELSR parameter evaluation %d outof %d  ... %d out of %d...\n', i2, length(FeaNumCandi), i1, length(paramCell));
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

save([prefix_mdcs, filesep, dataset, '_best_result_JELSR_liang_lle.mat'],'FeaNumCandi','res_gs','res_aio', 'res_gs_ps');
end