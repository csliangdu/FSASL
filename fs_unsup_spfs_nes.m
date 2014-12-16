function [ W, lam ]= fs_unsup_spfs_nes( X, Y, k, err, starting )
%unsupervised feature selection by 2-1 norm regression
% X - the training data, each row is an instance
% Y - the class label

if nargin < 5
    starting = 0.5;
end
% L2-1 norm
opts.q=2;

% lambda = lambda * lambda_{max}
opts.rFlag=1;

% norm( x_i - x_{i-1}, 2) <= .tol
% opts.tFlag = 3;

% Tolerance parameter.
% opts.tol=1e-4;

% opts.init=2;

% .x0= zeros(n,1), .c0=0
% opts.init=2;
opts.verbose = 0;
opts.maxIter = 500;

upL = 1; downL = 0;
lam = starting; % the initial search point
nZ = k + 2*err;
count = 1;
need = -1;

while abs(nZ - k) > err && count <= 10
    oldNZ = nZ;
    oldNeed = need;
    
%     fprintf('need %i, iteration: %2i, lam: %f\n', k, count, lam);
    W = mcLeastR(X, Y, lam, opts);
    opts.x0=W;
    nZ = sum(sum(W.^2,2)>0);
    if nZ - k > err
        need = -1;
        downL = lam; lam = (downL + upL) / 2;
    elseif nZ - k < -err
        need = 1;
        upL = lam; lam = (downL + upL) / 2;
    end
    if nZ < oldNZ && oldNeed == 1
        opts = rmfield(opts, 'x0');
        W = mcLeastR(X, Y, lam, opts);
        nZ = sum(sum(W.^2,2)>0);
%         fprintf('restart, %f, sel feat: %i\n-----\n', lam, nZ);
    end
%     fprintf('sel feat: %i\n-----\n', nZ);
    count = count + 1;
end