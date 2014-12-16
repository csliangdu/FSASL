function [ wFeat, SF ] = fs_unsup_spec( W, X, Y, Pram )
%function [ wFeat, SF ] = fsSpectrum( X, Y, Pram )
%   Select feature using the spectrum information of the graph laplacian
%   W - the similarity matrix or a kernel matrix
%   X - the input data, each row is an instance
%   Y - the labels of the data
%   Pram - the prameter of the algorithm
%   Pram.style - 1: unsupervised feature selection 2: supervised feature
%                         selection
%   Pram.expLam - the exp order for the eigenvalue
%   Pram.function - 1:f'Lf; 2:using all eigenvalue except the first one; 3:
%                             using the first k eigenvalues. (In this case
%                             the wieght the bigger the better.

[numInst,dimDat] = size(X);
if size(Y,2) > 1
    numC =size(Y,2);
else
    numC = length(unique(Y));
end

% build the degree matrix
D = diag(sum(W,2));
% build the laplacian matrix
L = D - W;

% D1 = D^(-0.5)
d1 = (sum(W,2)).^(-0.5);
d1(isinf(d1)) = 1;

% D2 = D^(0.5)
d2 = (sum(W,2)).^0.5;
v = diag(d2)*ones(numInst,1);
v = v/norm(v);
%  build the normalized laplacian matrix hatW = diag(d1)*W*diag(d1)
hatL = repmat(d1,1,numInst).*L.*repmat(d1',numInst,1);
if Pram.style ~=2
    hatL = (hatL'+hatL)/2;
end

% calculate and construct spectral information
switch Pram.style
    case 1,
        [ V, EVA ] = funG( hatL, Pram.expLam );
    case 2.
        [ V, EVA ] = funG( hatL, 1 );
end

% begin to select features
wFeat = [];

switch Pram.function
    case 1, % using f'Lf formulation
        for i = 1:dimDat
            f = X(:,i);
            hatF = diag(d2)*f;
            l = norm(hatF);
            
            if l < 100*eps
                wFeat(i) = 1000;
            else
                if Pram.style ~=2
                    hatF = hatF/l;
                end
                wFeat(i) = hatF'*hatL*hatF;
            end
        end
    case 2, % using all eigenvalues except the first one
        for i = 1:dimDat
            f = X(:,i);
            hatF = diag(d2)*f;
            l = norm(hatF);
            
            if l < 100*eps
                wFeat(i) = 1000;
            else
                hatF = hatF/l;
                wFeat(i) = hatF'*hatL*hatF/(1-(hatF'*v)^2);
            end            
        end
    case 3, % use the first k eigenvalues and the weight is the bigger the better.
        eva = diag(EVA);
        % calculate the eigenvalues
        switch Pram.style
            case 1,
                eva = eva.^(1/Pram.expLam);
                eva = 2 - eva;
                eva = eva.^(Pram.expLam);
            case 2,
                eva = max(eva) - eva;
        end
        
        for i = 1:dimDat
            % normalize the feature
            f = X(:,i);
            hatF = diag(d2)*f;
            l = norm(hatF);

            % calculate the weight
            if l < 100*eps
                wFeat(i) = -1;
            else
                hatF = hatF/l;
                alphas = hatF'*V(:,2:numC);
                wFeat(i) = (alphas.^2)*eva(2:numC);
            end
        end
end

SF = 1:dimDat;