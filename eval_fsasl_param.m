function r = eval_fsasl_param(p_name, param_candi, ids, fns, paramCell, res_aio)
% 
% r1 = eval_fsasl_param('lambda3',10.^[-5:5], [11:25], {'mean_acc', 'mean_nmi_sqrt', 'loocv'}, paramCell, res_aio);
% r2 = eval_fsasl_param('SLEPreg',[10.^-3, 0.005, 10.^-2, 0.05, 0.01], [11:25], {'mean_acc', 'mean_nmi_sqrt', 'loocv'}, paramCell, res_aio);
% r3 = eval_fsasl_param('lambda1',[0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99], [11:25], {'mean_acc', 'mean_nmi_sqrt', 'loocv'}, paramCell, res_aio);
% r1
% r2
% r3

if isvector(param_candi)
    param_candi = num2cell(param_candi);
end


r = zeros(length(param_candi), length(fns));
for i1 = 1:length(param_candi)
    tmp = []; % nP_a * 3
    for i2 = 1:size(res_aio, 1);
        if isfield(paramCell{i2, 1}, p_name) && strcmp(num2str(paramCell{i2,1}.(p_name)), num2str(param_candi{i1}))
            
            tmp2 = zeros(length(fns), length(ids));
            for i3 = 1:length(ids)
                tmp3 = zeros(length(fns),1);
                for i4 = 1:length(fns)
                    tmp3(i4) = res_aio{i2, ids(i3)}.(fns{i4});
                end
                tmp2(:,i3) = tmp3;
            end
            tmp = [tmp; mean(tmp2, 2)'];
        end
        
    end
    r(i1,:) = max(tmp, [], 1);
end

