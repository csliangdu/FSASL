function paramCell = fs_unsup_rufs_build_param(llkrrParamCell, alphaCandi, betaCandi, nuCandi)
n1 = length(alphaCandi);
n2 = length(betaCandi);
n3 = length(nuCandi);
n4 = length(llkrrParamCell);
nP = n1 * n2 * n3 * n4;
paramCell = cell(nP, 1);
idx = 0;
for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:n3
            for i4 = 1:n4
                param = [];
                param.alpha = alphaCandi(i1);
                param.beta = betaCandi(i2);
                param.nu = nuCandi(i3);
                param.MaxIter = 20;
                if param.alpha + param.beta + param.nu > 1e4
                    param.MaxIter = 5; % large parameter is costly for convergence
                end
                param.epsilon = 1e-2;
                param.verbose = 0;
                
                param.llkrrParam = llkrrParamCell{i4};
                idx = idx + 1;
                paramCell{idx} = param;
            end
        end
    end
end