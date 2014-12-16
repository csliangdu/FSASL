function [ V, D ] = funG( G, t )
%function [ V, D ] = funG( G, t )
%   modify the eigenvalue of a matrix by t order.

G = full(G);
[V,D] = eig(G);
d = diag(D);
% it is important to calculate G first
[d, orderIDX] = sort(d);
V = V(:,orderIDX);
d = d.^t;
d(isnan(d)) = 0;
d(isinf(d)) = 0;
D = diag(d);