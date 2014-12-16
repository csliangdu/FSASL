function [FeaNumCandi,res_gs,res_aio, res_gs_ps] = fs_unsup_fsasl_11_5_5_single_func(dataset, exp_settings, algo_settings)
%feature selection by FSASL

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
alphaCandi = 10.^[-5:5];
betaCandi = [0.1:0.2:0.9];
gammaCandi = [0.001, 0.005, 0.01, 0.05, 0.1];
maxIter = 50;
nnCandi = 5;
paramCell = fs_unsup_fsasl_build_param({'SLEP'}, {gammaCandi}, nnCandi, ...
    alphaCandi, betaCandi, {'LS21'}, maxIter);
%===============================================

disp('FSSL ...');
t_start = clock;
feaSubsets = cell(length(paramCell), 1);
parfor i1 = 1:length(paramCell)
    fprintf(['FSSL parameter search %d out of %d...\n'], i1, length(paramCell));
    param = paramCell{i1};
    param.lambda1 = 1 - param.lambda2;
    W = FSASL(X', nClass, param);
    [~, idx] = sort(sum(W.^2,2),'descend');
    % save([dataset,'\','feaIdx_param_', num2str(i1), '.mat'],'idx');
    feaSubsets{i1,1} = idx;
end
t_end = clock;
t1 = etime(t_end,t_start);
disp(['exe time: ',num2str(t1)]);

t_start = clock;
disp('evaluation ...');
res_aio = cell(length(paramCell), length(FeaNumCandi));
for i2 = 1:length(FeaNumCandi)
    parfor i1 = 1:length(paramCell)
        fprintf('FSASL parameter evaluation %d outof %d  ... %d out of %d...\n', i2, length(FeaNumCandi), i1, length(paramCell));
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


save(fullfile(prefix_mdcs, [dataset, '_best_result_FSSL_11_5_5.mat']),'FeaNumCandi','res_gs','res_aio', 'res_gs_ps', 'paramCell', 'feaSubsets');
end