function [ W, F, G ] = fs_unsup_rufs(X, L, G0, options)

% RUFS: Robust Unsupervised Feature Selection which solves the following
% problem:
%
% min || X - G * F ||_{2,1} + nu * tr(G' * L * G) +
%     alpha * || X * W - G ||_{2,1} + beta * || W ||_{2,1}
% s.t. F >= 0, G >= 0, G'G = I_c
% 
% We update both F and G by projected limited-memory BFGS method and W by
% limited-memory BFGS method.
%
% Input:
%
%    X: NSample x NFea data matrix.
%    L: NSample x NSample local learning regularization matrix.
%    G0: Initial NSample x NClus encoding matrix.
%    options: 1-by-1 structure for parameters.
%        options.nu: graph regularization parameter.
%        options.alpha: parameter on the regression term.
%        options.beta: parameter to control the row-sparsity of the 
%                      projection matrix.
%        options.MaxIter: maximal number of iterations, e.g., 10.
%        options.epsilon: convergence precision, e.g., 1e-4.
%        options.verbose: if verbose, 1 or 0.
%
% Output:
%
%    W: NFea x NClus projection matrix. 
%    F: NClus x NFea basis matrix.
%    G: NSample x NClus encoding matrix (relaxed). The i-th row corresponds
%       to the encoding vector for the i-th data point.
%
% Author: Mingjie Qian
% *************************************************************************
% Ref:
% Mingjie Qian, Chengxiang Zhai. 
% Robust Unsupervised Feature Selection. 
% The 23rd International Joint Conference on Artificial Intelligence (IJCAI), 2013.
% *************************************************************************

% Paramters setup

MaxIter = options.MaxIter;
epsilon = options.epsilon;
nu = options.nu;
a = options.alpha;
b = options.beta;
zeta = 1e7;
verbose = options.verbose;

L = nu * L;
if issparse(L)
    L = full(L);
end

% Algorithm

G = G0;
NClus = size(G0, 2);
NFea = size(X, 2);
I = eye(NClus);
F = ((G' * G) \ I) * G' * X;
F_pos = subplus(F);
F = F_pos + 0.2 * sum(sum(abs(F_pos))) / length(find(F_pos));
W = (X' * X + a \ b * eye(NFea)) \ (X' * G);
I_c = eye(NClus);
gamma = 1e-4;
mu = 1e-1;
ofv = sum(l2NormByRows(G * F - X)) + ...
      gamma * sum(l2NormByRows(F)) + ...
      mu * sum(G(:)) + ...
      trace(G' * L * G) + ...
      a * sum(l2NormByRows(X * W - G)) + ...
      b * sum(l2NormByRows(W)) + ...
      4 \ zeta * norm(G' * G - I_c)^2;
if verbose
    fprintf('\nInitial ofv: %g', ofv);
end

R = G * F - X;
r_f = l2NormByRows(R);
Grad_f1_F = G' * (R ./ repmat(r_f, [1, NFea]));
phi_F = l2Norm(F);
Grad_phi_F = F ./ repmat(phi_F + eps, [NClus, 1]);    
Grad_F = Grad_f1_F + gamma * Grad_phi_F;

Grad_f1_G = R * F' ./ repmat(r_f, [1, NClus]);
Grad_phi_G = mu + 2 * L * G;

H_h = G - X * W;
r_h = l2NormByRows(H_h);
psi_G = H_h ./ repmat(r_h, [1, NClus]);
Grad_G = Grad_f1_G + Grad_phi_G + a * psi_G + zeta * G * (G' * G - I_c);

Grad_W = a * X' * (H_h ./ repmat(r_h, [1, NClus])) + ...
         b * W ./ repmat(l2NormByRows(W) + eps, [1, NClus]);
initgrad = norm([Grad_F'; Grad_G; Grad_W], 'fro');

tol = epsilon * initgrad;
if verbose
    fprintf('\nInit gradient norm: %f\n', initgrad);
end
ind = 0;
while true
        
    % Fixing G, updating W
    [W, Grad_W] = UpdateW(X, G, a, b, W);
    
    % Fixing F and W, updating G
    [G, PGrad_G] = UpdateG(X, F, W, L, mu, a, G);
    
    % Fixing G, updating F
    [F, PGrad_F] = UpdateF(X, G, gamma, F);
    
    ind = ind + 1;
    if ind > MaxIter
        disp( 'Maximal iterations' );
        break;
    end
    
    ofv_pre = ofv;
    
    if verbose
        ofv = sum(l2NormByRows(G * F - X)) + ...
            gamma * sum(l2NormByRows(F)) + ...
            mu * sum(G(:)) + ...
            trace(G' * L * G) + ...
            a * sum(l2NormByRows(X * W - G)) + ...
            b * sum(l2NormByRows(W)) + ...
            4 \ zeta * norm(G' * G - I_c)^2;
        fprintf('Iter = %d ofv: %g', ind, ofv);
        if abs(ofv - ofv_pre) <= options.epsilon
            fprintf( '\nObjective value doesn''t decrease.\n' );
            break;
        end 
    end
    
    d = norm([PGrad_F'; PGrad_G; Grad_W]);
    if d <= tol
        if verbose
            fprintf('\nConverge successfully!');
            fprintf('\nIter = %d Final proj-grad norm: %f\n', ind, d);
        end
        break;
    end
    
    if mod(ind, 1) == 0 && verbose
        fprintf('\nIter = %d proj-grad norm: %f\n', ind, d);
    end
    
end

end

function [ G, PGrad ] = UpdateG(X, F, W, L, mu, s, G0)

NClus = size(F, 1);

% Projected L-BFGS method

% Parameter setup

alpha = 0.2;
beta = 0.5;
zeta = 1e7;
H0 = 1;
m = 10;
epsilon = 1e-1;
I_c = eye(NClus);
G = G0;

R = G * F - X;
r_f = l2NormByRows(R);
Grad_f1_G = R * F' ./ repmat(r_f, [1, NClus]);
Grad_phi_G = mu + 2 * L * G;
H_h = G - X * W;
r_h = l2NormByRows(H_h);
psi_G = H_h ./ repmat(r_h, [1, NClus]);
Grad_G = Grad_f1_G + Grad_phi_G + s * psi_G + zeta * G * (G' * G - I_c);

f1_G = sum(r_f);
phi_G = mu * sum(G(:)) + trace(G' * L * G) + s * sum(r_h) + 4 \ zeta * norm(G' * G - I_c, 'fro')^2;

G_x = Grad_G;

tol = epsilon * norm(Grad_G);

PG_x = zeros(size(G_x));
PHG_x = zeros(size(G_x));

fval_x = f1_G + phi_G;

s_ks = cell(m, 1);
y_ks = cell(m, 1);
rou_ks = cell(m, 1);
a = zeros(m, 1);

k = 0;
while true
    
    I_k = G_x < 0 | G > 0;
    I_k_com = not(I_k);
    PG_x(I_k) = G_x(I_k);
    PG_x(I_k_com) = 0;
    
    if norm(PG_x(:)) <= tol
        break;
    end
    
    if k == 0
        H = H0;
    else
        H = ((s_k(:)' * y_k(:)) / (y_k(:)' * y_k(:)));
    end
    
    q = G_x;
    for i = min(k, length(s_ks)):-1:1
        a(i) = rou_ks{i} * (s_ks{i}(:)' * q(:));
        q = q - a(i) * y_ks{i};
    end
    r = H * q;
    for i = 1:min(length(s_ks), k)
        b = rou_ks{i} * (y_ks{i}(:)' * r(:));
        r = r + (a(i) - b) * s_ks{i};
    end
    
    HG_x = r;
    I_k = HG_x < 0 | G > 0;
    I_k_com = not(I_k);
    PHG_x(I_k) = HG_x(I_k);
    PHG_x(I_k_com) = 0;
    
    if (PHG_x(:)' * G_x(:) <= 0)
        p = -PG_x;
    else
        p = -PHG_x;
    end
    
    t = 1;
    
    while true
        
        G_t = subplus(G + t * p);
        
        R_t = G_t * F - X;
        r_f_t = l2NormByRows(R_t);

        H_h_t = G_t - X * W;
        r_h_t = l2NormByRows(H_h_t);
        
        phi_G_t = mu * sum(G_t(:)) + trace(G_t' * L * G_t) + s * sum(r_h_t) + 4 \ zeta * norm(G_t' * G_t - I_c, 'fro')^2;
        fval_x_t = sum(r_f_t) + phi_G_t;
        
        if fval_x_t <= fval_x + alpha * (G_x(:)' * (G_t(:) - G(:)))
            break;
        else
            t = beta * t;
        end
        
    end
    
    G_pre = G;    
    G_x_pre = G_x;
    
    if abs(fval_x_t - fval_x) < epsilon
        break;
    end
    
    fval_x = fval_x_t;
    
    G = G_t;
    R = R_t;
    r_f = l2NormByRows(R);
    Grad_f1_G = R * F' ./ repmat(r_f, [1, NClus]);
    Grad_phi_G = mu + 2 * L * G;
    H_h = G - X * W;
    r_h = l2NormByRows(H_h);
    psi_G = H_h ./ repmat(r_h, [1, NClus]);
    Grad_G = Grad_f1_G + Grad_phi_G + s * psi_G + zeta * G * (G' * G - I_c);
    G_x = Grad_G;

    s_k = G - G_pre;
    y_k = G_x - G_x_pre;
    rou_k = 1 / (y_k(:)' * s_k(:));

    if k < m
        s_ks{k + 1} = s_k;
        y_ks{k + 1} = y_k;
        rou_ks{k + 1} = rou_k;
    else
        s_ks(1) = [];
        y_ks(1) = [];
        rou_ks(1) = [];
        s_ks{m} = s_k;
        y_ks{m} = y_k;
        rou_ks{m} = rou_k;
    end
    
    k = k + 1;
    
end

PGrad = PG_x;

end

function [ F, PGrad ] = UpdateF(X, G, gamma, F0)

NFea = size(X, 2);
NClus = size(G, 2);

% Projected L-BFGS method

% Parameter setup

alpha = 0.2;
beta = 0.5;
H0 = 1;
m = 10;
epsilon = 1e-1;

F = F0;
R = G * F - X;
r_f = l2NormByRows(R);
phi_F = l2Norm(F);

Grad_f1_F = G' * (R ./ repmat(r_f, [1, NFea]));
Grad_phi_F = F ./ repmat(phi_F + eps, [NClus, 1]);
G_x = Grad_f1_F + gamma * Grad_phi_F;

PG_x = zeros(size(G_x));
PHG_x = zeros(size(G_x));

fval_x = sum(r_f) + gamma * sum(phi_F);

s_ks = cell(m, 1);
y_ks = cell(m, 1);
rou_ks = cell(m, 1);
a = zeros(m, 1);

k = 0;
while true
    
    I_k = G_x < 0 | F > 0;
    I_k_com = not(I_k);
    PG_x(I_k) = G_x(I_k);
    PG_x(I_k_com) = 0;
        
    if k == 0
        H = H0;
    else
        H = ((s_k(:)' * y_k(:)) / (y_k(:)' * y_k(:)));
    end
    
    q = G_x;
    for i = min(k, length(s_ks)):-1:1
        a(i) = rou_ks{i} * (s_ks{i}(:)' * q(:));
        q = q - a(i) * y_ks{i};
    end
    r = H * q;
    for i = 1:min(length(s_ks), k)
        b = rou_ks{i} * (y_ks{i}(:)' * r(:));
        r = r + (a(i) - b) * s_ks{i};
    end
    
    HG_x = r;
    I_k = HG_x < 0 | F > 0;
    I_k_com = not(I_k);
    PHG_x(I_k) = HG_x(I_k);
    PHG_x(I_k_com) = 0;
    
    if (PHG_x(:)' * G_x(:) <= 0)
        p = -PG_x;
    else
        p = -PHG_x;
    end
    
    t = 1;
    
    while true
        
        F_t = subplus(F + t * p);
        
        R_t = G * F_t - X;
        r_f_t = l2NormByRows(R_t);
        phi_F_t = l2Norm(F_t);
        fval_x_t = sum(r_f_t) + gamma * sum(phi_F_t);
        
        if fval_x_t <= fval_x + alpha * (G_x(:)' * (F_t(:) - F(:)))
            break;
        else
            t = beta * t;
        end
        
    end
    
    F_pre = F;    
    G_x_pre = G_x;
    
    if abs(fval_x_t - fval_x) < epsilon
        break;
    end
    
    fval_x = fval_x_t;
    
    F = F_t;
    R = R_t;
    phi_F = phi_F_t;
    
    Grad_f1_F = G' * (R ./ repmat(r_f, [1, NFea]));
    Grad_phi_F = F ./ repmat(phi_F + eps, [NClus, 1]);
    G_x = Grad_f1_F + gamma * Grad_phi_F;
    
    s_k = F - F_pre;
    y_k = G_x - G_x_pre;
    rou_k = 1 / (y_k(:)' * s_k(:));
    
    if k < m
        s_ks{k + 1} = s_k;
        y_ks{k + 1} = y_k;
        rou_ks{k + 1} = rou_k;
    else
        s_ks(1) = [];
        y_ks(1) = [];
        rou_ks(1) = [];
        s_ks{m} = s_k;
        y_ks{m} = y_k;
        rou_ks{m} = rou_k;
    end
    
    k = k + 1;
    
end

PGrad = PG_x;

end

function [ W, Grad ] = UpdateW(X, G, s, v, W0)

NFea = size(X, 2);
NClus = size(G, 2);

% Projected L-BFGS method

% Parameter setup

alpha = 0.2;
beta = 0.5;
H0 = 1;
m = 10;
epsilon = 1e-1;

if isempty(W0)
    W = (X' * X + alpha \ beta * eye(NFea)) \ (X' * G);
else
    W = W0;
end

R = X * W - G;
r_h = l2NormByRows(R);
phi_W = l2NormByRows(W);
Grad_f1_W = X' * (R ./ repmat(r_h, [1, NClus]));
Grad_phi_W = W ./ repmat(phi_W + eps, [1, NClus]);
G_x = s * Grad_f1_W + v * Grad_phi_W;

fval_x = s * sum(r_h) + v * sum(phi_W);

s_ks = cell(m, 1);
y_ks = cell(m, 1);
rou_ks = cell(m, 1);
a = zeros(m, 1);

k = 0;
while true
    
    if norm(G_x(:)) <= epsilon
        break;
    end
    
    if k == 0
        H = H0;
    else
        H = ((s_k(:)' * y_k(:)) / (y_k(:)' * y_k(:)));
    end
    
    q = G_x;
    for i = min(k, length(s_ks)):-1:1
        a(i) = rou_ks{i} * (s_ks{i}(:)' * q(:));
        q = q - a(i) * y_ks{i};
    end
    r = H * q;
    for i = 1:min(length(s_ks), k)
        b = rou_ks{i} * (y_ks{i}(:)' * r(:));
        r = r + (a(i) - b) * s_ks{i};
    end
    
    p = -r;
    
    t = 1;
    z = G_x(:)' * p(:);
    while true
        
        W_t = W + t * p;
        R_t = X * W_t - G;
        r_h_t = l2NormByRows(R_t);
        phi_W_t = l2NormByRows(W_t);
        fval_x_t = s * sum(r_h_t) + v * sum(phi_W_t);
        
        if fval_x_t <= fval_x + alpha * t * z
            break;
        else
            t = beta * t;
        end
        
    end
    
    W_pre = W;    
    G_x_pre = G_x;
    
    if abs(fval_x_t - fval_x) < 1e-1
        break;
    end
    
    fval_x = fval_x_t;
    
    W = W_t;
    R = R_t;
    r_h = r_h_t;
    phi_W = phi_W_t;
    Grad_f1_W = X' * (R ./ repmat(r_h, [1, NClus]));
    Grad_phi_W = W ./ repmat(phi_W + eps, [1, NClus]);
    G_x = s * Grad_f1_W + v * Grad_phi_W;
    
    s_k = W - W_pre;
    y_k = G_x - G_x_pre;
    rou_k = 1 / (y_k(:)' * s_k(:));
    
    if k < m
        s_ks{k + 1} = s_k;
        y_ks{k + 1} = y_k;
        rou_ks{k + 1} = rou_k;
    else
        s_ks(1) = [];
        y_ks(1) = [];
        rou_ks(1) = [];
        s_ks{m} = s_k;
        y_ks{m} = y_k;
        rou_ks{m} = rou_k;
    end
    
    k = k + 1;
    
end

Grad = G_x;

end

function [ res ] = l2NormByRows( X )
% l2NormByRows Calculate the l_2 norm of a vector or row vectors of a matrix.
%
% Input:
%
%    X: A vector or a matrix.
%
% Output:
%
%    res: The l_2 norm of a vecor or row vectors of a matrix.
%
% Author: Mingjie Qian
% Date: Sep. 23rd, 2012

if size(X, 1) == 1
    res = norm(X, 2);
else
    res = sum(X.^2, 2).^(0.5);
end

end

function [ res ] = l2Norm( X )
% L2Norm Calculate the l_2 norm of a vector or column vectors of a matrix.
%
% Input:
%
%    X: A vector or a matrix.
%
% Output:
%
%    res: The l_2 norm of a vector or column vectors of a matrix.
%
% Author: Mingjie Qian
% Date: Sep. 23rd, 2012

if size(X, 2) == 1
    res = norm(X, 2);
else
    res = sum(X.^2, 1).^(0.5);
end

end