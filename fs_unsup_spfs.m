function fList = fs_unsup_spfs(X, K, Y, numF, options)
% A wrapper function for different solvers of SPFS
%
% [ fList ] = spfs_sfs(X, K, numF);
%
% [ W, lam ]= spfs_nes( X, Y, k, err, starting );
%
% [ fList, W ] = spfs_larnes( X, Y, numF );
%
% [ fList, W ] = spfs_lar( X, Y, numF )
%
% each solver is downloaded from the author Zheng Zhao
% https://sites.google.com/site/alanzhao/
%
% [1] Efficient Spectral Feature Selection with Minimum Redundancy, AAAI 2010
% [2] On Similarity Preserving Feature Selection, TKDE, 2013

if ~exist('options', 'var') || ~isfield(options, 'spfs_type')
    options.spfs_type = 'SFS';
end

switch lower(options.spfs_type)
    case lower('SFS')
        [ fList ] = fs_unsup_spfs_sfs(X, K, numF);
    case lower('LAR')
        error('not supported yet!');
        % the following code with LAR did not return enough features
        [eigvec, eigval] = eigs(K, options.nClass, 'LA');
        Y = eigvec * diag(sqrt(max(diag(eigval), eps)));
        [ fList, W ] = fs_unsup_spfs_lar( X, Y, numF );
    case lower('LARNES')
        error('not supported yet!');
        % the following code with LARNES did not return enough features
        [eigvec, eigval] = eigs(K, options.nClass, 'LA');
        Y = eigvec * diag(sqrt(max(diag(eigval), eps)));
        [ fList, W ] = fs_unsup_spfs_larnes( X, Y, numF );
    case lower('NES')
        [eigvec, eigval] = eigs(K, options.nClass, 'LA');
        Y = eigvec * diag(sqrt(max(diag(eigval), eps)));
        [ W, lam ]= fs_unsup_spfs_nes( X, Y, numF, 0.1*numF);
        fList = sum(W.^2,2);
        [~, fList] = sort(fList, 'descend');
        fList = fList(1:numF);
    otherwise
        error('not supported yet!');
end