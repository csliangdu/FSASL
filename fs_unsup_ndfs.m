function idx = fs_unsup_ndfs(X, nClass, param)

[L, F_init] = NDFS_init(X, nClass, param);

W_init = ones(size(X,2),nClass);   %W: the feature selection matrix
warning off;
[F,W,obj]=NDFS_iter(X', L, F_init, W_init, param.maxiter, param.alpha, param.beta, param.gamma);
warning on;
[~, idx] = sort(sum(W.*W,2),'descend');
end

function [L, F_init] = NDFS_init(X, nClass, param)
%construct the affinity matrix
S = constructW(X, struct('k', param.k, 'WeightMode', param.weightMode, 't', param.t));
diag_ele_arr = sum(S);
diag_ele_arr_t = diag_ele_arr.^(-1/2);
L = eye(size(X,1)) - diag(diag_ele_arr_t)* S *diag(diag_ele_arr_t);
L = (L + L')/2;
[eigvec, eigval] = eig(L);
[~, t1] = sort(diag(eigval), 'ascend');
eigvec = eigvec(:, t1(1:nClass));
eigvec = bsxfun(@rdivide, eigvec, sqrt(sum(eigvec.^2,2) + eps));

%init F and W
rand('twister',5489); %#ok
label = litekmeans(eigvec,nClass,'Replicates',10); % significantly!
F_init = rand(size(X,1),nClass);
for i = 1:size(X,1)
    F_init(i,label(i)) = 1;
end
F_init = F_init + 0.2;
end

function [F,W,obj]=NDFS_iter(X,L,F,W,maxIter,alpha,beta,gamma)
%	X: Rows of vectors of data points
%	L: The laplacian matrix.
%   F: the cluster result
%   W: the feature selection matrix

if nargin == 0
    return; 
end

[nFeat,nSamp] = size(X);

if size(L,1) ~= nSamp
    error('L is error');
end
XX=X*X';

Wi = sqrt(sum(W.*W,2)+eps);
d = 0.5./Wi;
D = diag(d);

% G=inv(XX+beta*D);
% W=G*X*F;
% Wi = sqrt(sum(W.*W,2)+eps);
% d = 0.5./Wi;
% D = diag(d);
% clear Wi
% M=L+alpha*(eye(nSamp)-X'*G*X);
% clear G
% M=(M+M')/2;
% F = F.*(gamma*F + eps)./(M*F + gamma*F*F'*F + eps);
% F = F*diag(sqrt(1./(diag(F'*F)+eps)));

iter=1;
while iter<=maxIter %|| (iter>2&& obj(end-1)-obj(end)>10^(-3)*obj(end))
    G=inv(XX+beta*D);
    W=G*X*F;
    Wi = sqrt(sum(W.*W,2)+eps);
    d = 0.5./Wi;
    D = diag(d);
    clear Wi
    M=L+alpha*(eye(nSamp)-X'*G*X);
    clear G
    M=(M+M')/2;

    F = F.*(gamma*F + eps)./(M*F + gamma*F*F'*F + eps);
    F = F*diag(sqrt(1./(diag(F'*F)+eps)));
    clear Wnew   
    
    obj(iter)=trace(F'*M*F)+gamma/4*norm(F'*F-eye(size(F,2)),'fro')^2;
    iter=iter+1;
    
end
end