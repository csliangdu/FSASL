function [feature_idx, feature_score, subset_score] = fs_unsup_traceratio(Sb, Sw, feature_num)
% Sb: a matrix to reflects the between-class or global affinity
%     relationship encoded on Graph, Sb = X*Lb*X'
% Sw: a matrix to reflects the within-class or local affinity relationship
%     encoded on Graph, Sw = X*Lw*X'
% feature_idx: the ranked feature index based on subset-level score
% feature_score: the feature-level score
% subset_score: the subset-level score


sb = abs(diag(Sb));
sw = abs(diag(Sw));
sw(find(sw == 0)) = 0.000000000000001; %#ok

% preprocessing.
t_fnum = length(sb);
[fs, fs_idx] = sort(sb./sw,'descend');

para = 0.9;

u_fnum = floor(para*t_fnum);
sb = sb(fs_idx(1:u_fnum));
sw = sw(fs_idx(1:u_fnum));


ind = 1:feature_num;
k = sum(sb(ind))/sum(sw(ind));
for i = 1: 20
    [score, I] = sort(sb - k*sw, 'descend');
    ind = I(1:feature_num);
    old_k = k;
    k = sum(sb(ind))/sum(sw(ind));
    if abs(k - old_k) < 0.000000000001
        break;
    end;
end
I = fs_idx(I);

feature_idx = I;
feature_score = score;
subset_score = k;