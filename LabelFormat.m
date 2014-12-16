function Y = LabelFormat(y)
[~,~,y] = unique(y);
nClass = max(y);
Y = zeros(length(y), nClass);
Y(sub2ind(size(Y), 1:length(y), y')) = 1;