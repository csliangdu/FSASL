function [res_gs2, res_gs_ps2] = grid_search_fs(res_aio)
% input
%	nParam * nFea, cell
%
% for each feature subset,
%     for each evaluation measure,
%         choose the bset result
% Liang Du (csliangdu@gmail.com)
%

[nParam, nSubsets] = size(res_aio);
res_gs = cell(1, nSubsets);
res_gs_ps = res_gs;
fn1 = {'mean_acc', 'mean_nmi_sqrt', 'mean_nmi_max', 'mean_purity', 'mean_prec', 'mean_recall', 'mean_f1'};
fn2 = {'std_acc', 'std_nmi_sqrt', 'std_nmi_max', 'std_purity', 'std_prec', 'std_recall', 'std_f1'};
fn3 = {'best_obj_acc', 'best_obj_nmi_max', 'best_obj_nmi_sqrt', 'best_obj_purity',...
    'best_obj_prec', 'best_obj_recall', 'best_obj_f1', ...
    'jac', 'loocv'};
fn4 = {'red'};
for i1 = 1:nSubsets
    res_gs{1, i1} = res_aio{1,i1}; %place holder
    for i3 = 1:length(fn1)
        res_gs_ps{1, i1}.(fn1{i3}) = 1; 
    end
    for i3 = 1:length(fn3)
        res_gs_ps{1, i1}.(fn3{i3}) = 1; 
    end
    for i3 = 1:length(fn4)
        res_gs_ps{1, i1}.(fn4{i3}) = 1; 
    end
    for i2 = 1:nParam
        for i3 = 1:length(fn1)
            if (isfield(res_aio{i2, i1}, fn1{i3}) && isfield(res_gs{1, i1},fn1{i3}) ) && (res_aio{i2, i1}.(fn1{i3}) > res_gs{1, i1}.(fn1{i3}))
                res_gs{1, i1}.(fn1{i3}) = res_aio{i2, i1}.(fn1{i3});
                res_gs{1, i1}.(fn2{i3}) = res_aio{i2, i1}.(fn2{i3});
                res_gs_ps{1, i1}.(fn1{i3}) = i2; 
            end
        end
        for i3 = 1:length(fn3)
            if (isfield(res_aio{i2, i1}, fn3{i3}) && isfield(res_gs{1, i1}, fn3{i3}) ) && (res_aio{i2, i1}.(fn3{i3}) > res_gs{1, i1}.(fn3{i3}))
                res_gs{1, i1}.(fn3{i3}) = res_aio{i2, i1}.(fn3{i3});
                res_gs_ps{1, i1}.(fn3{i3}) = i2; 
            end
        end
        for i3 = 1:length(fn4)
            if (isfield(res_aio{i2, i1}, fn4{i3}) && isfield(res_gs{1, i1}, fn4{i3}) ) && (res_aio{i2, i1}.(fn4{i3}) > res_gs{1, i1}.(fn4{i3}))
                res_gs{1, i1}.(fn4{i3}) = res_aio{i2, i1}.(fn4{i3});
                res_gs_ps{1, i1}.(fn4{i3}) = i2;
            end
        end
    end
end

res_gs2 = res_gs{1,1};
res_gs_ps2 = res_gs_ps{1,1};
for i1 = 2:nSubsets
	for i3 = 1:length(fn1)
		res_gs2.(fn1{i3}) = [res_gs2.(fn1{i3}), res_gs{1, i1}.(fn1{i3})];
		res_gs2.(fn2{i3}) = [res_gs2.(fn2{i3}), res_gs{1, i1}.(fn2{i3})];
        res_gs_ps2.(fn1{i3}) = [res_gs_ps2.(fn1{i3}), res_gs_ps{1, i1}.(fn1{i3})];
	end
	for i3 = 1:length(fn3)
		res_gs2.(fn3{i3}) = [res_gs2.(fn3{i3}), res_gs{1, i1}.(fn3{i3})];
        res_gs_ps2.(fn3{i3}) = [res_gs_ps2.(fn3{i3}), res_gs_ps{1, i1}.(fn3{i3})];
	end
	for i3 = 1:length(fn4)
		res_gs2.(fn4{i3}) = [res_gs2.(fn4{i3}), res_gs{1, i1}.(fn4{i3})];
        res_gs_ps2.(fn4{i3}) = [res_gs_ps2.(fn4{i3}), res_gs_ps{1, i1}.(fn4{i3})];
	end
end
end