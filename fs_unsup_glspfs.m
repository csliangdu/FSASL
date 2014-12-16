function [feaIndx,W,obj] = fs_unsup_glspfs(X, Kmatrix, L, r1, r2, numFea) %% 
[num, dim] = size(X);
d = ones(dim,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L = computeM(X,Kmatrix,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[UY,VY] = eig(Kmatrix);
diagVal = diag(VY);
indxPos = find(diagVal>eps);
UYpos = UY(:,indxPos);
VYpos = diag(sqrt(diagVal(indxPos)));
Ypos = UYpos*VYpos;

NIter = 20;
flag =1;
objold = inf;
iter = 0;
if num<dim    
    while flag
        iter = iter +1;
        D = spdiags(d,0,dim,dim);
        DX = D*X';
        %%% Notive trick!!!
        W = DX*(((eye(num)+ r2*L)*X*DX + r1*eye(num))\Ypos);
        Xi = sqrt(sum(W.*W,2));
        d = 2*Xi;     
        XW = X*W -Ypos;
        obj(iter) = trace(XW*XW') + r2*trace(L*((X*W)*(X*W)')) + r1*sum(Xi);        
        if abs((objold-obj(iter))/obj(iter)) <1e-4 || iter>NIter
            flag = 0;
        end
        objold = obj(iter);
    end
else
    while flag
        iter = iter +1;
        D = spdiags(d,0,dim,dim);
        DX = D*X';
        %%% Notive trick!!!
        W = (DX*(eye(num)+ r2*L)*X + r1*eye(dim) )\(DX*Ypos);     
        Xi = sqrt(sum(W.*W,2));
        d = 2*Xi;    
        XW = X*W -Ypos;
        obj(iter) = trace(XW*XW') + r2*trace(L*((X*W)*(X*W)')) + r1*sum(Xi);        
        if abs((objold-obj(iter))/obj(iter)) <1e-4 || iter>NIter
            flag = 0;
        end
        objold = obj(iter);
    end  
end
Xi = sqrt(sum(W.*W,2));
[val0,indx0] = sort(Xi,'descend');
feaIndx = indx0(1:numFea);