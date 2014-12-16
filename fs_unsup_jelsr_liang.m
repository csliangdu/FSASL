function [model_jelsr] = fs_unsup_jelsr_liang(X, param)

[nDim, nSmp] = size(X);

if ~exist('param', 'var');
    param = [];
end

if ~isfield(param, 'nClusters');
    error('The number of Clusters should be specified');
else
    nClusters = param.nClusters;
end

if isfield(param, 'k')
    k = param.k;
else
    k = 5;
end

if isfield(param, 'beta')
    beta = param.beta;
else
    beta = 1;
end

if isfield(param, 'alpha')
    alpha = param.alpha;
else
    alpha = 1;
end

t1 = cputime;
L = computeLocalStructure(X', param.weightMode, param.k, param.t);
[W, Y, obj] = JELSR_AlterOptimizer(X, L, nClusters, alpha, beta);

model_jelsr.z = sqrt(sum(W.^2,2));
model_jelsr.Y = Y;
model_jelsr.runTime = cputime - t1;
model_jelsr.obj = obj;
end

function [W, Y, obj] = JELSR_AlterOptimizer(X, L, nDimEmb, alpha, beta)
% Input
%         X: nDim * nSmp
%         L: nSmp * nSmp; Local reconstruction kernel
%         nDimEmb: low embedding dimension
%         alpha: regularization parameter
%         beta: regularization parameter
% Output
%         W: nDim * nEmb
%         Y: nEmb * nSmp
%         obj: obj history
% Optimization objective
%         min{W, U, Y} = tr(Y L Y') + beta*||W' X - Y ||^2 + beta*alpha* tr(W' U W)
%
% [1]. Feature Selection via Joint Embedding Learning and Sparse Regression.
% Chenping Hou, etc. IJCAI, 2011.
%

[nDim, nSmp] = size(X);

if nDim < nSmp
    A = X*X';
end
U = ones(nDim, 1);

nIter = 20;
obj = [];
epsilon = 1e-2;

for iter = 1:nIter
    % Step1: Fix U, update Y by solving the problem in Eq. (16);
    
    if nDim < nSmp
        % AiX = inv(A + alpha*U)*X;
        AiX = (A + alpha*diag(U))\X;
    else
        % AiX = alpha * U^-1 X [ I - (alpha I + X' U^-1 X)^-1 X' U^-1 X]
        UX = bsxfun(@times, 1./U, X);
        KiU = X' * UX;
        AiX = UX * (eye(nSmp) - (alpha * eye(nSmp) + KiU) \ KiU);
        AiX = AiX/alpha;
    end
    
    K = L + beta*eye(nSmp) - beta*X'*AiX;
    K = (K + K') / 2;
    [eigvec, eigval] = eig(K);
    [eigval, idx] = sort(diag(eigval));
    Y = eigvec(:, idx(1:nDimEmb));
    
    % Step2: Fix U, update W by using Eq. (13);
    W = AiX*Y;
    
    % Step3: FixW, update U by Eq. (9);
    U = full(0.5./(sqrt(sum(W.^2,2)) + eps));
    
    % obj(end+1) = trace(Y'*L*Y) + beta*sum(sum( (X'*W - Y).^2)) + beta*alpha*sum(sqrt(sum(W.^2,2)));
    %
    %     if iter > 1 && abs(obj(end) - obj(end-1))/abs(obj(end)) < epsilon;
    %         break;
    %     end
end
end