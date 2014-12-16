function paramCell = fs_unsup_jelsr_liang_build_param(r1Candi, r2Candi, knnCandi)
n1 = length(r1Candi);
n2 = length(r2Candi);
n3 = length(knnCandi);
nP = n1 * n2 * n3;
paramCell = cell(nP, 1);
idx = 0;
for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:n3
            param = [];
            param.r1 = r1Candi(i1);
            param.r2 = r2Candi(i2);
            param.r3 = knnCandi(i3);
            idx = idx + 1;
            paramCell{idx} = param;
        end
    end
end
end