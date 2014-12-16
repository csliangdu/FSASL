function paramCell = fs_unsup_mcfs_build_param(knnCandi, weightCandi, weight_param_Candi)
n1 = length(knnCandi);
n2 = length(weightCandi);
n3 = zeros(n2, 1);
for i1 = 1:length(weightCandi)
    n3(i1) = max(1, length(weight_param_Candi{i1}));
end
nP = n1 * max(sum(n3), 1) ;
paramCell = cell(nP, 1);
idx = 0;
for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:max(n3(i2), 1)
            param = [];
            param.k = knnCandi(i1);
            param.weightMode = weightCandi{i2};
            if ~isempty(weightCandi) && ~isempty(weight_param_Candi{i2})
                tmp = weight_param_Candi{i2};
                param.t = tmp(i3);
            else
                param.t = 1; % place holder
            end
            idx = idx + 1;
            paramCell{idx} = param;
        end
    end
end