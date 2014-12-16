function paramCell = fs_unsup_spec_build_param(kernelParamCell, styleCandi, expLamCandi, funcCandi)
%   Pram.style - 1: unsupervised feature selection 2: supervised feature
%                         selection
%   Pram.expLam - the exp order for the eigenvalue
%   Pram.function - 1:f'Lf; 2:using all eigenvalue except the first one; 3:
%                             using the first k eigenvalues. (In this case
%                             the wieght the bigger the better.
if ~exist('styleCandi', 'var')
	kernelParamCell = {};
end
if ~exist('styleCandi', 'var') || isempty(styleCandi)
	styleCandi = [1];
end

if ~exist('expLamCandi', 'var') || isempty(expLamCandi)
	expLamCandi = [0.25, 1, 4];
end

if ~exist('funcCandi', 'var') || isempty(funcCandi)
	funcCandi = [1, 2, 3];
end

n0 = max(length(kernelParamCell), 1);
n1 = length(styleCandi);
n2 = length(expLamCandi);
n3 = length(funcCandi);
nP = n0 * n1 * n2 * n3;
paramCell = cell(nP, 1);
idx = 0;
for i0 = 1:n0
for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:n3
            param = [];
            if ~isempty(kernelParamCell)
            	param.kernelOption = kernelParamCell{i0};
            end
            param.style = styleCandi(i1);
            param.expLam = expLamCandi(i2);
            param.function = funcCandi(i3);
            idx = idx + 1;
            paramCell{idx} = param;
        end
    end
end
end