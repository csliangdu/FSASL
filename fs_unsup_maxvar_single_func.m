function [FeaNumCandi,res_gs,res_aio, res_gs_ps] = fs_unsup_maxvar_single_func(dataset, exp_settings, algo_settings)
%use laplacian score to select features.

%======================setup===========================
FeaNumCandi = exp_settings.FeaNumCandi;
nKmeans = exp_settings.nKmeans;
prefix_mdcs = [];
if isfield(exp_settings, 'prefix_mdcs')
    prefix_mdcs = exp_settings.prefix_mdcs;
end
%======================================================

%================setup======================

%===========================================
disp(['dataset:',dataset]);
[X, Y] = extractXY(dataset);
[nSmp,nDim] = size(X);

%get maxvar score
disp('get maxvar score...');
t_start = clock;
FeaScore = fs_unsup_maxvar(X);
[~, index] = sort(FeaScore, 'descend');
% save([dataset, filesep,'feaIdx.mat'],'index');
t_end = clock;
t1 = etime(t_end,t_start);
disp(['exe time: ',num2str(t1)]);

t_start = clock;
disp('evaluation ...');
res_aio = cell(1, length(FeaNumCandi)); 
parfor feaIdx = 1:length(FeaNumCandi)
    res_aio{1, feaIdx} = evalUnSupFS(X, Y, index(1:FeaNumCandi(feaIdx)), struct('nKm', nKmeans));
end
[res_gs, res_gs_ps] = grid_search_fs(res_aio);
res_gs.feaset = FeaNumCandi;
t_end = clock;
t2 = etime(t_end,t_start);
disp(['exe time: ',num2str(t2)]);
res_gs.time = t1;
res_gs.time2 = t2;

save([prefix_mdcs, filesep, dataset, '_best_result_MaxVar.mat'],'FeaNumCandi','res_gs','res_aio', 'res_gs_ps');
end