function [Lap, S] = computeLocalStructure(X, type, k, sigma, emb_dim)
% Input
%     X, n * nDim
%     type, 'LPP', 'LLE', 'LTSA'
%     k, neighborhood size, needed by all the three types, 5 by default;
%     sigma, gaussian kernel bandwidth, optSigma(X), by default, only used by LPP
%     emb_dim, embedding dimension, only used by LTSA
%

if ~exist('k', 'var')
    k = 5;
end

if ~exist('type', 'var')
    type = 'LPP';
end

if strcmp(type, 'LPP') && (~exist('sigma', 'var') || isempty(sigma))
    sigma = optSigma(X);
end

if strcmp(type, 'LTSA') && (~exist('emb_dim', 'var') || isempty(emb_dim))
    emb_dim = 2;
end

[n, d] = size(X);
switch lower(type)
    case lower('LPP')
        % Construct neighborhood graph
        % disp('Constructing neighborhood graph...');
        if size(X, 1) < 4000
            G = L2_distance(X', X');
            % Compute neighbourhood graph
            [tmp, ind] = sort(G);
            for i=1:size(G, 1)
                G(i, ind((2 + k):end, i)) = 0;
            end
            G = sparse(double(G));
            G = max(G, G');             % Make sure distance matrix is symmetric
        else
            G = find_nn(X, k);
        end
        G = G .^ 2;
        G = G ./ max(max(G));
        
        % Compute weights (W = G)
        % disp('Computing weight matrices...');
        
        % Compute Gaussian kernel (heat kernel-based weights)
        G(G ~= 0) = exp(-G(G ~= 0) / (sigma ^ 2));
        
        % Construct diagonal weight matrix
        D = diag(sum(G, 2));
        
        % Compute Laplacian
        L = D - G;
        L(isnan(L)) = 0; D(isnan(D)) = 0;
        L(isinf(L)) = 0; D(isinf(D)) = 0;
        Lap = L;
        S = G;
    case lower('LLE')
%         neighborhood = zeros(n,k);
        Dist = EuDist2(X);

%         for ii =1:n
%             index00 = setdiff(1:n,ii);
%             [sorted,index] = sort(Kmatrix(ii,index00),2,'descend');
%             neighborhood(ii,:) = index00(index(1:k));
%         end
        [~, neighborhood] = sort(Dist, 2, 'ascend');
        neighborhood = neighborhood(:,2:k+1);
        if(k > d)
            tol=1e-3; % regularlizer in case constrained fits are ill conditioned
        else
            tol=1e-12;
        end
        
        W = zeros(k,n);
        for ii=1:n
            z = X(neighborhood(ii,:),:)-repmat(X(ii,:),k,1); % shift ith pt to origin
            C = z*z';                                        % local covariance
            C = C + eye(size(C))*tol*trace(C);                   % regularlization
            W(:,ii) = C\ones(k,1);                           % solve Cw=1
            W(:,ii) = W(:,ii)/sum(W(:,ii));                  % enforce sum(w)=1
        end
        
        M = sparse(1:n,1:n,ones(1,n),n,n,4*k*n);
        for ii=1:n
            w = W(:,ii);
            jj = neighborhood(ii,:)';
            M(ii,jj) = M(ii,jj) - w'; %#ok
            M(jj,ii) = M(jj,ii) - w;%#ok
            M(jj,jj) = M(jj,jj) + w*w';%#ok
        end
        M = max(M,M');
        M = sparse(M);
        % For sparse datasets, we might end up with NaNs or Infs in M. We just set them to zero for now...
        M(isnan(M)) = 0;
        M(isinf(M)) = 0;
        Lap = M;
        S = sparse(repmat(1:n, k, 1), neighborhood(:), W(:), n, n, n*k);
    case lower('LTSA')
        % Compute neighborhood indices
        % disp('Find nearest neighbors...');
        n = size(X, 1);
        [D, ni] = find_nn(X, k);
        
        % Compute local information matrix for all datapoints
        % disp('Compute local information matrices for all datapoints...');
        Bi = cell(1, n);
        for i=1:n
            % Compute correlation matrix W
            Ii = ni(i,:);
            Ii = Ii(Ii ~= 0);
            kt = numel(Ii);
            Xi = X(Ii,:) - repmat(mean(X(Ii,:), 1), [kt 1]);
            W = Xi * Xi';
            W = (W + W') / 2;
            
            % Compute local information by computing d largest eigenvectors of W
            [Vi, Si] = schur(full(W));
            [s, Ji] = sort(-diag(Si));
            if length(Ji) < emb_dim
                emb_dim = length(Ji);
                % warning(['Target dimensionality reduced to ' num2str(emb_dim) '...']);
            end
            Vi = Vi(:,Ji(1:emb_dim));
            
            % Store eigenvectors in G (Vi is the space with the maximum variance, i.e. a good approximation of the tangent space at point Xi)
            % The constant 1/sqrt(kt) serves as a centering matrix
            Gi = double([repmat(1 / sqrt(kt), [kt 1]) Vi]);
            
            % Compute Bi = I - Gi * Gi'
            Bi{i} = eye(kt) - Gi * Gi';
        end
        
        % Construct sparse matrix B (= alignment matrix)
        % disp('Construct alignment matrix...');
        B = speye(n);
        for i=1:n
            Ii = ni(i,:);
            Ii = Ii(Ii ~= 0);
            B(Ii, Ii) = B(Ii, Ii) + Bi{i};							% sum Bi over all points
            B(i, i) = B(i, i) - 1;
        end
        B = (B + B') / 2;											% make sure B is symmetric
        
        % For sparse datasets, we might end up with NaNs in M. We just set them to zero for now...
        B(isnan(B)) = 0;
        B(isinf(B)) = 0;
        Lap = B;
        S = [];
    otherwise
        Lap = [];
        S = [];
        disp('not supported yet!');
end