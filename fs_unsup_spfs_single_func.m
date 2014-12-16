function [FeaNumCandi,res_gs,res_aio, res_gs_ps] = fs_unsup_spfs_single_func(dataset, exp_settings, algo_settings)
%feature selection by SPFS

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
s1 = optSigma(X);
tCandi = 2.^[-3:3] * s1.^2;
spfs_typeCandi = {'SFS', 'NES'};
nP = length(tCandi) * length(spfs_typeCandi);
paramCell = cell(nP, 1) ;
idx = 0;
for i1 = 1:length(tCandi)
    for i2 = 1:length(spfs_typeCandi)
        param = [];
        param.t = tCandi(i1);
        param.spfs_type = spfs_typeCandi{i2};
        idx = idx + 1;
        paramCell{idx} = param;
    end
end
%===========================================

disp('get SPFS...');
t_start = clock;
Dist = EuDist2(X, X, 0);

feaSubsets = cell(length(paramCell), 1);
for i1 = 1:length(paramCell)
    param = paramCell{i1};
    param.nClass = nClass;
	K = exp( - Dist / param.t);
    index = fs_unsup_spfs(X, K, [], max(FeaNumCandi), param);
	feaSubsets{i1,1} = index;
end
t_end = clock;
t1 = etime(t_end,t_start);
disp(['exe time: ',num2str(t1)]);

t_start = clock;
disp('evaluation ...');
res_aio = cell(length(paramCell), length(FeaNumCandi));
for i2 = 1:length(FeaNumCandi)
    for i1 = 1:length(paramCell)
        fprintf('SPFS parameter evaluation %d outof %d  ... %d out of %d...\n', i2, length(FeaNumCandi), i1, length(paramCell));
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

save(fullfile(prefix_mdcs, [dataset, '_best_result_SPFS.mat']),'FeaNumCandi','res_gs','res_aio', 'res_gs_ps');
end