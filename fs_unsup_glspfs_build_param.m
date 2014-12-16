function paramCell = fs_unsup_glspfs_build_param(local_type_candi, local_type_param_candi, knn_size_candi, ...
    lambda1_candi, lambda2_candi, global_kernel_cell_candi)
n1 = length( local_type_candi );
n2 = zeros(n1, 1);
for i1 = 1:length( local_type_candi )
    n2(i1) = max(1, length(local_type_param_candi{i1}));
end
n3 = length( knn_size_candi );
n4 = length( lambda1_candi );
n5 = length( lambda2_candi );
n6 = length( global_kernel_cell_candi );

nP = max(sum(n2), 1) * n3 * n4 * n5 * n6;

paramCell = cell(nP, 1);
idx = 0;
for i1 = 1:n1
    for i2 = 1:max(n2(i1), 1)
        for i3 = 1:n3
            for i4 = 1:n4
                for i5 = 1:n5
                    for i6 = 1:n6
                        param = [];
                        param.local_type = local_type_candi{i1};
                        if ~isempty(local_type_candi) && ~isempty(local_type_param_candi{i1})
                            tmp = local_type_param_candi{i1};
                            param.local_lpp_sigma = tmp(i2);
                            param.local_ltsa_embedded_dim = tmp(i2);
                        else
                            param.local_lpp_sigma = []; %place holder
                            param.local_ltsa_embedded_dim = [];%place holder
                        end
                        param.local_k = knn_size_candi(i3);
                        param.lambda1 = lambda1_candi(i4);
                        param.lambda2 = lambda2_candi(i5);
                        param.global_kernel_option = global_kernel_cell_candi{i6};
                        idx = idx + 1;
                        paramCell{idx} = param;                        
                        
                    end
                end
            end
        end
    end
end
end