function T = localLearnMx_KRR( X, param)

% conpute K via rbf function, the width is computed by self-tunning

K = constructW(X, struct('WeightMode', 'HeatKernel', 'k', param.nNeighbors)); 

[nSmp, nDim] = size(X);

% locate neighbors for each data
W = 1*(K>0);

% compute the local learning matrices
if param.nNeighbors < nSmp - 1 && param.nNeighbors > 0
    % compute A by local regularized kernel ridge regression
    A = zeros( nSmp, nSmp );
    for n = 1 : nSmp
        idxV = find( W( n, : ) > 0 );
        A( n, idxV ) = K( n, idxV )*inv( K(idxV, idxV) + param.rLambda * eye( length( idxV ) ) );
    end

    % matrix T
    T = eye( nSmp ) - A;
    T = T' * T;

else  % all the data are neighboring to each other
    A = [];  % A can not be computed directly
    I = eye( nSmp );

    % deformed kernel
    T = K * inv( K + param.rLambda * I );

    T = I - T;
    T = inv( diag( diag( T )  ) ) * T;
    T = T' * T;
end