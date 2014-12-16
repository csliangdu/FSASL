function [W, S, A, objHistory] = FSASL(X, nClass, options)
if ~exist('options', 'var')
    options = [];
end

% Optios for global structure learning
if ~isfield(options, 'lambda1')
    options.lambda1 = 1; % [need to search]
end

if ~isfield(options, 'LassoType')
    options.LassoType = 'SLEP';
end

if ~isfield(options, 'SLEPrFlag')
    options.SLEPrFlag = 1; % the input parameter 'ReguAlpha' is a ratio in (0, 1)
end

if ~isfield(options, 'SLEPreg')
    options.SLEPreg = 0.01; % [need to search, and fix it]
end

if ~isfield(options, 'LARSk')
    options.LARSk = 5; % [need to search, and fix it]
end

if ~isfield(options, 'LARSratio')
    options.LARSratio = 2;
end

% Optios for local structure learning
if ~isfield(options, 'lambda2')
    options.lambda2 = 1; % [need to search] aim to show local structure is helpful
end

if ~isfield(options, 'Localk')
    options.Localk = 5; % [need to search, and fix it]
end

if ~isfield(options, 'LocalReg')
    options.LocalReg = estimateReg(X, options.Localk); % aim to avoid search
end

% Optios for subspace learning
if ~isfield(options, 'GroupLassoType')
    options.GroupLassoType = 'LS21';
end

if ~isfield(options, 'lambda3')
    options.lambda3 = 1; % [need to search
end

if ~isfield(options, 'maxiter')
    options.maxiter = 1; % [need to search
end
% options.lambda1 = 1 - options.lambda2;
[~, nSmp] = size(X);
X2 = X;
objHistory = [];
for iter = 1:options.maxiter
    
    S = zeros(nSmp);
    if options.lambda1 > 0 && ( options.maxiter < 5 || iter > 1)
        % update global structure LG
        for iSmp = 1:nSmp
            candIdx = ones(nSmp, 1);
            candIdx(iSmp) = 0;
            candIdx = candIdx > 0;
            switch lower(options.LassoType)
                case lower('SLEP')
                    S(candIdx, iSmp) = LeastR(X2(:, candIdx), X2(:, iSmp), options.SLEPreg, struct('rFlag', options.SLEPrFlag, 'rsL2', 0));
                case lower('LARS')
                    S(candIdx, iSmp) = LassoLARS(X2(:, candIdx), X2(:, iSmp), options.LARSk * options.LARSratio, 'verbose', 0);
                case lower('lars2')
                    Gram = X2(:, candIdx) * X2(:, candIdx)';
                    Gram = max(Gram,Gram');
                    S(candIdx, iSmp) = lars(X2(:, candIdx), X2(:, iSmp),'lasso', -(max(options.LARSk)+5),1,Gram,options.LARSk);
                case lower('lars3')
                    S(candIdx, iSmp) = lars(X2(:, candIdx), X2(:, iSmp),'lasso', -(max(options.LARSk)+5),0,[],options.LARSk);
                otherwise
                    error('method does not exist!');
            end
        end
        LG = (eye(nSmp) - S);
        LG = LG * LG';
        LG = (LG + LG') / 2;
    else
        LG = 0;
    end
    
    A = zeros(nSmp);
    if options.lambda2 > 0
        if iter > 1
            % update local structure LL
            distx = L2_distance_1(X2, X2);
            if iter>0
                [~, idx] = sort(distx,2);
            end;
            
            for iSmp = 1 : nSmp
                if options.Localk < nSmp
                    idxa0 = idx(iSmp, 2: options.Localk + 1);
                else
                    idxa0 = 1 : nSmp;
                end;
                dxi = distx(iSmp, idxa0);
                ad = - (dxi) / (2 * options.LocalReg);
                A(iSmp, idxa0) = EProjSimplex_new(ad);
            end;
        else
            A = constructW(X2', struct('k', options.Localk));
        end
        A = (A+A')/2;
        LL = diag(sum(A)) - A;
        LL = (LL + LL') / 2;
    else
        LL = 0;
    end
    
    L = options.lambda1 * LG + options.lambda2 * LL;
    
    % update embedding
    
    switch lower(options.GroupLassoType)
        case lower('JFSSL')
            Y = eig1(L, nClass, 0);
            W = FSSL_subspace(X, Y, options.lambda3);
        case lower('LS21')
            Y = eig1(L, nClass, 0);
            W = LS21(X', Y, options.lambda3);
        case lower('NDFS') % d^3
            tmp = X * L * X';
            tmp = (tmp + tmp') / 2;
            if exist('W', 'var')
                W = LquadR21_reg(tmp, nClass, options.lambda3, W);
            else
                W = LquadR21_reg(tmp, nClass, options.lambda3);
            end
        case lower('UDFS')
            tmp = X * L * X';
            tmp = (tmp + tmp') / 2;
            W = LquadR21_reg(tmp, nClass, options.lambda3);
        case lower('MCLEASTR')
            Y = eig1(L, nClass, 0);
            W = mcLeastR(X', Y, options.lambda3, struct('rFlag', 1, 'rsL2', 0));
        otherwise
            error('method does not exist!');
    end
    X2 = W' * X;
    obj = trace(X2 * L * X2');
    if options.lambda1 > 0 && strcmpi(options.LassoType, 'SLEP')
        obj = obj + options.lambda1 * options.SLEPreg * sum(sum(abs(S)));
    end
    
    if options.lambda2 > 0
        obj = obj + options.lambda2 * options.LocalReg * sum(sum(A.^2));
    end
    
    obj = obj + options.lambda3 * sum(sqrt(sum(W.^2, 2)));
    objHistory = [objHistory; obj]; %#ok
    %
end
end

function A = FSSL_subspace(X, Y, regu)
[d, ~] = size(X);
[n, nClass] = size(Y);
% Check the solutions
nSolutionCheck = 0;
r1 = rank(X');
r2 = rank([X', Y]);
if r1 == r2 && r1 < d
    % X'*A = Y has many solution == rank(X') == rank([X', Y]) < d
    nSolutionCheck = 1;
end

A = zeros(d, nClass);
% Step 2: Find A satisfies the linear system
nIter = 20;
if nSolutionCheck
    % Situation 1, Infinitely many solutions
    G = eye(d);
    for iter = 1:nIter
        Gi = inv(G);
        A = Gi*X*inv(X'*Gi*X)*Y; %#ok
        normG = sqrt(sum(A.^2,2));
        nzIdx = (normG ~= 0);
        dd = zeros(d, 1);
        dd(nzIdx) = 1./normG;
        G = diag(dd);
    end
else
    % Situation 1, Single or No solution
    G = eye(d);
    for iter = 1:nIter
        Gi = inv(G);
        A = Gi*X*inv(X'*Gi*X + 0.5/regu*eye(n))*Y; %#ok
        normG = sqrt(sum(A.^2,2));
        nzIdx = (normG ~= 0);
        dd = zeros(d, 1);
        dd(nzIdx) = 1./normG;
        G = diag(dd);
    end
end
end

function [x, ft] = EProjSimplex_new(v, k)
%
% Problem
%
%  min  1/2 || x - v||^2
%  s.t. x>=0, 1'x=1
%

if nargin < 2
    k = 1;
end;

ft=1;
n = length(v);

v0 = v-mean(v) + k/n;
%vmax = max(v0);
vmin = min(v0);
if vmin < 0
    f = 1;
    lambda_m = 0;
    while abs(f) > 10^-10
        v1 = v0 - lambda_m;
        posidx = v1>0;
        npos = sum(posidx);
        g = -npos;
        f = sum(v1(posidx)) - k;
        lambda_m = lambda_m - f/g;
        ft=ft+1;
        if ft > 100
            x = max(v1,0); %#ok
            break;
        end;
    end;
    x = max(v1,0);
    
else
    x = v0;
end;
end

% compute squared Euclidean distance
% ||A-B||^2 = ||A||^2 + ||B||^2 - 2*A'*B
function d = L2_distance_1(a,b)
% a,b: two matrices. each column is a data
% d:   distance matrix of a and b



if (size(a,1) == 1)
    a = [a; zeros(1,size(a,2))];
    b = [b; zeros(1,size(b,2))];
end

aa=sum(a.*a); bb=sum(b.*b); ab=a'*b;
d = repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab;

d = real(d);
d = max(d,0);

% % force 0 on the diagonal?
% if (df==1)
%   d = d.*(1-eye(size(d)));
% end

end

function r = estimateReg(X, k)
[d, nSmp] = size(X);
distX = L2_distance_1(X,X);
%distX = sqrt(distX);
[distX1, idx] = sort(distX,2);
A = zeros(nSmp);
rr = zeros(nSmp,1);
for i = 1:nSmp
    di = distX1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    A(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end;
r = mean(rr);
end


function [X, obj]=LquadR21_reg(A, k, r, X0)
% quadratic loss with 21-norm regularization
%  min_{X'*X=I}  Tr(X'*A*X) + r * ||X||_21


NIter = 36;
[m n] = size(A);
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
end

function [eigvec, eigval, eigval_full] = eig1(A, c, isMax, isSym)

if nargin < 2
    c = size(A,1);
    isMax = 1;
    isSym = 1;
elseif c > size(A,1)
    c = size(A,1);
end;

if nargin < 3
    isMax = 1;
    isSym = 1;
end;

if nargin < 4
    isSym = 1;
end;

if isSym == 1
    A = max(A,A');
end;
try
    [v, d] = eig(A);
    d = diag(d);
    %d = real(d);
catch
    if isMax == 0 
        [v, d] = eigs(sparse(A), c, 'sa', struct('tol', 1e-5'));
    else
        [v, d] = eigs(sparse(A), c, 'la', struct('tol', 1e-5'));
    end
end

if isMax == 0
    [d1, idx] = sort(d);
else
    [d1, idx] = sort(d,'descend');
end;
idx1 = idx(1:c);
eigval = d(idx1);
eigvec = v(:,idx1);

eigval_full = d(idx);
end


function W = LS21(X, Y, r, W0)
[n, m] = size(X);
if nargin < 4
    d = ones(m,1);
else
    Wi = sqrt(sum(W0.^2,2)+eps);
    d = 0.5./(Wi);
end;

maxiter = 10;
if n < d % n^3
    XY = X' * Y;
    for iter= 1:maxiter
        rd = 1 ./ (r * d);
        Xrd = bsxfun(@times, X, rd');
        XrdX = Xrd * X';
        A = diag(rd) - Xrd' / (eye(n) + XrdX) * Xrd;
        W = A * XY;
        Wi = sqrt(sum(W.^2,2)+eps);
        d = 0.5./(Wi);
    end
else % d^3
    XX = X' * X;
    XY = X' * Y;
    for iter= 1:maxiter
        W = (XX + r * diag(d)) \ XY;
        Wi = sqrt(sum(W.^2,2)+eps);
        d = 0.5./(Wi);
    end
end
end
