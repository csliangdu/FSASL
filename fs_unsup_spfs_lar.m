function [ fList, W ] = fs_unsup_spfs_lar( X, Y, numF )
% function [ fList W ] = spfs_lar( X, K, numF )
%   X - the data, each row is an instance
%   Y - the response of nY column
%   numF - the number of features we want to selected

[nD, nF] = size(X);
nY = size(Y,2);

W = zeros(nF, nY);
fList = zeros(numF, 1);
k = 1; R = Y;

% find the most correlated one
bestCor = -1; bestNor = 0;
for i = 1:nF
    curF = X(:,i);
    curNorm = norm(curF'*R,2);
    if curNorm > bestNor
        bestCor = i; bestNor = curNorm;
    end
end
fList(k) = bestCor; XA = X(:, bestCor);

while k < numF && k < nF && k < nD
    k = k + 1;
%     fprintf('%i,',k);
    
    % obtain the proceed direction
    GA = XA\R;
    
    % compute how far can we go for every f
    a = X(:,fList(1))'*R;
    b = X(:,fList(1))'*XA*GA;
    bestCor = -1; bestNor = inf;
    for i = 1:nF
        if sum(fList==i) > 0
            continue;
        end
        c = X(:,i)'*R;
        d = X(:,i)'*XA*GA;
        p1=b*b'-d*d'; p2 = a*b'-c*d'; p3 = a*a'-c*c';
        s1 = (p2+abs(sqrt(p2^2-p1*p3)))/p1;
        s2 = (p2-abs(sqrt(p2^2-p1*p3)))/p1;
        if (s1<=0 || s1>1)
            s1 = 100;
        end
        if (s2<=0 || s2>1)
            s2 = 100;
        end
        if s1==100 && s2==100
            continue;
        else
            s = min(s1,s2);
        end
        if s < bestNor
            bestNor = s;
            bestCor = i;
        end
    end
    if bestCor == -1;
        return
    else
        fList(k) = bestCor;
        XA = X(:, fList(1:k));
        W(fList(1:k-1),:) = W(fList(1:k-1),:) + bestNor*GA;
        R = Y - X*W;
        % fprintf(' R: %f, W: %f, l: %f\n',norm(R), norm(W), bestNor);
    end
end

GA = pinv(full(XA'*XA))*XA'*R;
W(fList(1:k),:) = W(fList(1:k),:) + GA;
R = Y - X*W;
% fprintf(' R: %f, W: %f\n',norm(R), norm(W));