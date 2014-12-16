function [ fList ] = fs_unsup_spfs_sfs(X, K, numF)
% function [ fList ] = spfs_sfs(X, K, numF)
%   X - data, each row is an instance
%   K - the similarity matrix of instances
%   numF - the number of features to be selected

nF = size(X,2);
fList = zeros(numF,1);
R = K;
count = 1;
while count <= numF && count <= nF
    %     fprintf('%i,',count);
    %     if mod(count,10)==0
    %         fprintf('\n');
    %     end
    
    [R, selF] = find_best_match(X, fList, R);
    if selF == -1
        return;
    else
        fList(count) = selF;
    end
    count = count + 1;
end
end

function [ newR, selF ] = find_best_match(X, fList, R)
nF = size(X,2);
newR = R;
selF = -1;
%     smallestErr = norm(newR,'fro'); modified
smallestErr = inf;
for i = 1:nF
    if sum(fList == i)>0
        continue;
    end
    curF = X(:,i);
    curErr = (curF'*curF)^2-2*curF'*R*curF;
    if smallestErr >= curErr
        newR = R - curF*curF';
        smallestErr = curErr;
        selF = i;
    end
end

end