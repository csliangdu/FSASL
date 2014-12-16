function [Sb, Sw] = compute_SbSw_unsup(X, nK)
% X: training data each row is a data;
% calculate L_b and L_w defined in Laplacian score
% Sb = X*L_b*X';
% Sw = X*L_w*X';
if nargin < 2
    nK = 5;
end
W = constructW(X, struct('k', nK));
Dw = sum(W,2);
L_w = diag(Dw) - W;
L_b = (Dw * Dw') / sum(Dw);

L_w = (L_w + L_w')/2;
L_b = (L_b + L_b')/2;

Sb = X'*L_b*X;
Sw = X'*L_w*X;

% very important!
Sb = (Sb + Sb')/2;
Sw = (Sw + Sw')/2;
