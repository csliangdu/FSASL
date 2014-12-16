function paramCell = fs_unsup_fsasl_build_param(sr_solver_candi, sr_solver_param_candi, knn_size_candi, ...
    lambda2_candi, lambda3_candi, fs_solver_candi, iter_candi)
n1 = length( sr_solver_candi );
n2 = zeros(n1, 1);
for i1 = 1:length( sr_solver_candi )
    n2(i1) = max(1, length(sr_solver_param_candi{i1}));
end
n3 = length( knn_size_candi );
n4 = length( lambda2_candi );
n5 = length( lambda3_candi );
n6 = length( fs_solver_candi );
n7 = length( iter_candi );

nP = max(sum(n2), 1) * n3 * n4 * n5 * n6 * n7;

paramCell = cell(nP, 1);
idx = 0;
% for i0 = 1: n0
for i1 = 1:n1
    for i2 = 1:max(n2(i1), 1)
        for i3 = 1:n3
            for i4 = 1:n4
                for i5 = 1:n5
                    for i6 = 1:n6
                        for i7 = 1:n7
                            param = [];
                            param.LassoType = sr_solver_candi{i1};
                            if ~isempty(sr_solver_candi) && ~isempty(sr_solver_param_candi{i1})
                                tmp = sr_solver_param_candi{i1};
                                param.SLEPreg = tmp(i2);
                                param.LARSk = tmp(i2);
                            end
                            param.Localk = knn_size_candi(i3);
                            param.lambda2 = lambda2_candi(i4);
                            param.lambda1 = 1;
                            param.lambda3 = lambda3_candi(i5);
                            param.GroupLassoType = fs_solver_candi{i6};
                            param.maxiter = iter_candi(i7);
                            idx = idx + 1;
                            paramCell{idx} = param;
                            
                        end
                    end
                end
            end
        end
    end
end
end