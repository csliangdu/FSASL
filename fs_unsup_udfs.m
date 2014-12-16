function [X, obj]=fs_unsup_udfs(A, k, r, X0)
% quadratic loss with 21-norm regularization
%
%  min_{X'*X=I}  Tr(X'*A*X) + r * ||X||_21
% 

NIter = 20;
[m, n] = size(A); %#ok
if nargin < 4
    d = ones(n,1);
else
    Xi = sqrt(sum(X0.*X0,2)+eps);
    d = 0.5./(Xi);
end;

for iter = 1:NIter
    D = diag(d);
    M = A+r*D;
    M = max(M,M');
    [evec, eval] = eig(M);
    eval = diag(eval);
    [~, idx] = sort(eval);
    X = evec(:,idx(1:k));
    
    Xi = sqrt(sum(X.*X,2)+eps);
    d = 0.5./(Xi);
    
    obj(iter) = trace(X'*A*X) + r*sum(Xi); %#ok
end;