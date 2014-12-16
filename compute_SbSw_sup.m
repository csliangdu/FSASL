function [Sb, Sw, L_b, L_w] = calculate_L(X,Y)
% X: training data each row is a data;
% Y: labels
% calculate L_b and L_w defined in traditional LDA
% Sb = X*L_b*X';
% Sw = X*L_w*X';


% ====== Initialization
data_n = size(X, 1);
class_set = unique(Y);
class_n = length(class_set);

W = zeros(data_n);
% U 的每一列为一个样本的标签向量，设该样本为第i类，则标签向量的第i个元素为1，其他为0
for i = 1:class_n
    U = (Y == class_set(i)); 
    count = sum(U);	% Cardinality of each class
    index = find(U==1);
    W(index,index) = 1/count; 
end;     
L_w = eye(data_n) - W;
L_t = eye(data_n) - 1/data_n*ones(data_n,1)*ones(1,data_n);
L_b = L_t - L_w;

L_w = (L_w + L_w')/2;
L_b = (L_b + L_b')/2;

Sb = X'*L_b*X;
Sw = X'*L_w*X;

% very important!
Sb = (Sb + Sb')/2;
Sw = (Sw + Sw')/2;
