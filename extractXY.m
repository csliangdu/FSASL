function [X, Y] = extractXY(dataset)
% For single view data : 'X', 'y'
load(dataset);
if ~exist('X', 'var') && exist('fea', 'var')
    X = fea;
end

if exist('Y', 'var') && size(Y,1) < size(Y,2)
    Y = Y';
end

if exist('Y', 'var') && min(size(Y)) > 1
    Y = LabelFormat(Y);
end

if ~exist('Y', 'var') && exist('gnd', 'var')
    Y = gnd(:);
end

if ~exist('Y', 'var') && exist('y', 'var')
    Y = y(:);
end
end