function paramCell = buildParam_LLKRR(knnCandi, rLamdaCandi)
n1 = length(knnCandi);
n2 = length(rLamdaCandi);
nP = n1 * n2;
paramCell = cell(nP, 1);
idx = 0;
for i1 = 1:n1
    for i2 = 1:n2
        param = [];
        param.nNeighbors = knnCandi(i1);
        param.rLambda = rLamdaCandi(i2);
        idx = idx + 1;
        paramCell{idx} = param;
    end
end