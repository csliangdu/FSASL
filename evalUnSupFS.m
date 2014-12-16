function fs_res = evalUnSupFS(X, Y, feaIdx, options)
% Evaluate the selected features
%
% [1] redundancy
% [2] jac, k=5
% [3] acc, mean and std, best_obj
% [4] nmi(max_version), mean, std, best_obj
% [5] purity, mean, std, best_obj
% [6] precision, mean, std, best_obj
% [7] recall, mean, std, best_obj
% [8] f1, mean, std, best_obj
% [9] loocv, knn, k=1
%
%
% Reference
% [1] On Similarity Preserving Feature Selection, TKDE, 2011
%
% Liang Du (csliangdu@gmail.com)
%

if ~exist('options', 'var')
    options = [];
end

if ~isfield(options, 'jac_k')
    options.jac_k = 5;
end

if ~isfield(options, 'nKm')
    options.nKm = 10;
end

if ~isfield(options, 'knn_k')
    options.knn_k = 1;
end

[nSmp, nDim] = size(X);
Xsub = X(:, feaIdx);

fs_red = compute_RED(Xsub);
fs_jac = compute_JAC(X, Xsub, options.jac_k);
fs_loocv = compute_loocv(Xsub, Y, options.knn_k);

fs_cluster = compute_Clustering(Xsub, Y, options.nKm);

fs_res = struct('red', fs_red, 'jac', fs_jac, 'loocv', fs_loocv);
fs_res = cell2struct([struct2cell(fs_res);struct2cell(fs_cluster)],[fieldnames(fs_res);fieldnames(fs_cluster)]);
end

function fs_red = compute_RED(Xsub)
[nSmp, nDim] = size(Xsub);

if ~isempty('corr') && nDim < 2000
    C1 = corr(Xsub);
    sum_corr = sum(sum(tril(C1, -1)));
else
    mX = mean(Xsub, 1);
    stdX = std(Xsub, 0, 1);
    Xsub = bsxfun(@minus, Xsub, mX);
    sum_corr = 0;
    for i1 = 1:nDim
        for i2 = 1:i1-1
            sum_corr = sum_corr + (Xsub(:,i1)' * Xsub(:, i2)) / (stdX(i1) * stdX(i2) + eps);
        end
    end
end
fs_red = sum_corr / (nDim * (nDim -1) + eps);
end

function fs_jac = compute_JAC(X, Xsub, k)
D1 = EuDist2(X, X, 0);
[~, Idx1] = sort(D1, 2, 'ascend');
Idx1 = Idx1(:, 2:k+1);
Idx1 = mat2cell(Idx1, ones(size(X,1), 1), k);
D2 = EuDist2(Xsub, Xsub, 0);
[~, Idx2] = sort(D2, 2, 'ascend');
Idx2 = Idx2(:, 2:k+1);
Idx2 = mat2cell(Idx2, ones(size(X,1), 1), k);
s1 = cellfun(@union, Idx1, Idx2, 'UniformOutput', 0);
s2 = cellfun(@intersect, Idx1, Idx2, 'UniformOutput', 0);
n1 = cellfun(@length, s1);
n2 = cellfun(@length, s2);
fs_jac = mean(n2 ./ n1);
end

function fs_cluster = compute_Clustering(Xsub, Y, nKm)
if ~exist('nKm', 'var')
    nKm = 20;
end
nClass = length(unique(Y));
acc_list = zeros(nKm, 1);
nmi_max_list = zeros(nKm, 1);
nmi_sqrt_list = zeros(nKm, 1);
purity_list = zeros(nKm, 1);
obj_list = zeros(nKm, 1);
prec_list = zeros(nKm, 1);
recall_list = zeros(nKm, 1);
f1_list = zeros(nKm, 1);
rand('twister',5489); %#ok
for iKm = 1:nKm
    [label, ~, ~, sumD] = litekmeans(Xsub, nClass,'Replicates',1);
    tmp_res = evalClustering(Y, label);
    acc_list(iKm) = tmp_res.acc;
    nmi_max_list(iKm) = tmp_res.nmi_max;
    nmi_sqrt_list(iKm) = tmp_res.nmi_sqrt;
    purity_list(iKm) = tmp_res.purity;
    obj_list(iKm) = sum(sumD);
    prec_list(iKm) = mean(tmp_res.precision);
    recall_list(iKm) = mean(tmp_res.recall);
    f1_list(iKm) = mean(tmp_res.f1);
end
[~, idx] = min(obj_list);
fs_cluster = struct('mean_acc', mean(acc_list), 'std_acc', std(acc_list), ...
    'mean_nmi_max', mean(nmi_max_list), 'std_nmi_max', std(nmi_max_list), ...
    'mean_nmi_sqrt', mean(nmi_sqrt_list), 'std_nmi_sqrt', std(nmi_sqrt_list), ...
    'mean_purity', mean(purity_list), 'std_purity', std(purity_list), ...
    'mean_prec', mean(prec_list), 'std_prec', std(prec_list), ...
    'mean_recall', mean(recall_list), 'std_recall', std(recall_list), ...
    'mean_f1', mean(f1_list), 'std_f1', std(prec_list), ...
    'best_obj_acc', acc_list(idx(1)), 'best_obj_nmi_max', nmi_max_list(idx(1)),...
    'best_obj_nmi_sqrt', nmi_sqrt_list(idx(1)), 'best_obj_purity', purity_list(idx(1)), ...
    'best_obj_prec', prec_list(idx(1)), 'best_obj_recall', recall_list(idx(1)),...
    'best_obj_f1', f1_list(idx(1)), ...
    'aio_acc', acc_list, 'aio_nmi_max', nmi_max_list, 'aio_nmi_sqrt', nmi_sqrt_list, 'aio_purity', purity_list,...
    'aio_prec', prec_list, 'aio_recall', recall_list, 'aio_f1', f1_list);
end

function fs_loocv = compute_loocv(Xsub, Y, k)
if ~exist('k', 'var')
    k = 1;
end
Dist = EuDist2(Xsub,Xsub,0);
[~, Idx] = sort(Dist, 2, 'ascend');
idx = Idx(:, 2);
label = Y(idx);
fs_loocv = mean(label == Y);
end
