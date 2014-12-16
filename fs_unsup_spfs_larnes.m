function [ fList W ] = fs_unsup_spfs_larnes( X, Y, numF )
% function [ fList W ] = spfs_lar( X, K, numF )
%   X - the data, each row is an instance
%   Y - the response of nY column
%   numF - the number of features we want to selected

[nD, nF] = size(X);
nY = size(Y,2);

W = zeros(nF, nY);

R = Y;

% find the most correlated one
nor = X'*R;
nor = sqrt(sum((nor.*nor),2));
[bestNor, bestCor] = max(nor);

fList = bestCor; k = length(fList);
cnt = 0;

while k < numF && k < nF && k < nD
    cnt = cnt + 1;
    
    % obtain the proceed direction
    XA = X(:, fList);
    GA = XA\R;
    
    % compute how far can we go for every f to reduce lambda
    a = X(:,fList(1))'*R;
    bestCor = -1; bestNor = inf;
    for i = 1:nF
        if sum(fList==i) > 0
            continue;
        end
        c = X(:,i)'*R;
        d = X(:,i)'*XA*GA;
        p1=a*a'-d*d'; p2 = a*a'-c*d'; p3 = a*a'-c*c';
        bb = p2^2-p1*p3;
        if bb < 0
            continue;
        end
        bb = sqrt(bb);
        s1 = (p2+abs(bb))/p1;
        s2 = (p2-abs(bb))/p1;
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
        % reduce the size of labmda and update W with nes-L2,1        
        W(fList,:) = W(fList,:) + bestNor*GA;
        R = Y - X*W;
        lam = norm(X(:,fList(1))'*R,2);
                
        % find the nes-L2,1 solution
        [ fList, WA ] = nes(X, Y, W, [fList bestCor], lam*0.995);
        W(fList,:) = WA;
        R = Y - X*W; k = length(fList);
        disp(' ');
%         fprintf('step: %5i, feature: %5i, Lambda:%f\n',cnt+1, k, lam);
%         fprintf('----------------------------------\n');
    end
end

% R = Y - X*W;
% lam = norm(X(:,fList(1))'*R,2);
% opts.q=2;
% opts.tol=1e-6;
% opts.maxIter = 10000;
% opts.x0=W;
% W = mcLeastR(X, Y, lam, opts);
% fList = find(sum(abs(W),2));

    function [newfList WAA] = nes(X, Y, W, fList, lam)
        trd = 10e-5;

        WAA = W(fList,:); XAA = X(:,fList); newfList = fList;
        
        opts.q=2;
        opts.tol=1e-7;
        opts.maxIter = 10000;
        
        stop = 0;
        maxC = 1000; counterr = 1;
        
        % obtain a solution on XAA
        while stop == 0 && counterr <= maxC
            LC = setdiff(1:nF,newfList);
            
            opts.x0=WAA;
            WAA = mcLeastR(XAA, Y, lam, opts);
            
            keepIDX = find(sum(abs(WAA),2));
            
            newfList = newfList(keepIDX); 
            WAA = WAA(keepIDX,:);
            XAA = XAA(:,keepIDX);
            
            RR = Y - XAA*WAA;
            pp = X(:,LC)'*RR; pp = sqrt(sum(pp.*pp,2)); [maxr sel] = max(pp);
            if maxr - lam >= trd
%                 fprintf('find %i voilations\n', length(find( (pp-lam) > trd )));
            end
            
            if maxr - lam < trd
                if length(keepIDX) < length(fList)
                    lam = lam*0.995;
                else
                    stop = 1;
                end
            else
                aaa = [newfList,LC((pp-lam) > trd )]; newfList = aaa;
                aaa = [WAA;zeros(length(find( (pp-lam) > trd )),size(WAA,2))]; WAA = aaa;
                XAA = X(:,newfList);                
                counterr = counterr + 1;
            end
        end
    end % end function nes
end