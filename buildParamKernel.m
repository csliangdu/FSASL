function paramCell = buildParamKernel(p1Candidates, p2Candidates, p3Candidates)
if ~exist('p1Candidates', 'var')
    p1Candidates = {'Linear', 'PolyPlus', 'Polynomial', 'Gaussian'};
end
if ~exist('p2Candidates', 'var')
    p2Candidates = {[], [2, 4], [2, 4], [0.01, 0.05, 0.1, 1, 10, 50, 100]};
end
if ~exist('p3Candidates', 'var')
    p3Candidates = {'Sample-Scale'};
end

n1 = max(1, length(p1Candidates));
n2 = zeros(n1, 1);
for i2 = 1:length(p1Candidates)
    n2(i2) = max(1, length(p2Candidates{i2}));
end
n3 = max(1, length(p3Candidates));

nP = max(sum(n2), 1) * n3;

paramCell = cell(nP, 1);
idx = 0;
for i1 = 1:n1
    for i2 = 1:max(n2(i1), 1)
        for i3 = 1:n3
            param = [];
            if ~isempty(p1Candidates) && iscell(p1Candidates)
                param.KernelType = p1Candidates{i1};
            end
            
            if ~isempty(p1Candidates) && ~isempty(p2Candidates{i1})
                tmp = p2Candidates{i1};
                switch param.KernelType
                    case 'Gaussian'
                        param.t = tmp(i2);
                    case {'Polynomial',  'PolyPlus'}
                        param.d = tmp(i2);
                end
            end
            
            if ~isempty(p3Candidates) && iscell(p3Candidates)
                param.normType = p3Candidates{i3};
            end
            
            idx = idx + 1;
            paramCell{idx} = param;
        end
    end
end
end