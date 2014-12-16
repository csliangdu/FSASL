function paramCell = fs_unsup_udfs_build_param(knnCandi, gammaCandi, lamdaCandi)
n1 = length(knnCandi);
n2 = length(gammaCandi);
n3 = length(lamdaCandi);
nP = n1 * n2 * n3;
paramCell = cell(nP, 1);
idx = 0;
for i1 = 1:n1
    for i2 = 1:n2
		for i3 = 1:n3
			param = [];
			param.k = knnCandi(i1);
			param.gamma = gammaCandi(i2);
			param.lamda = lamdaCandi(i3);
			idx = idx + 1;
			paramCell{idx} = param;
		end
    end
end
