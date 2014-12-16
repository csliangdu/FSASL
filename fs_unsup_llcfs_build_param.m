function paramCell = fs_unsup_llcfs_build_param(nClusters, kCandidates, betaCandidates, kTypeCandidates, maxiterCandidates, epsilonCandidates )
if ~exist('kTypeCandidates', 'var') || isempty(kTypeCandidates)
	kTypeCandidates = [1];
end

if ~exist('maxiterCandidates', 'var') || isempty(maxiterCandidates)
	maxiterCandidates = [20];
end

if ~exist('epsilonCandidates', 'var') || isempty(epsilonCandidates)
	epsilonCandidates = [1e-4];
end


n1 = length( kCandidates );
n2 = length( betaCandidates );
n3 = length( kTypeCandidates );
n4 = length( maxiterCandidates );
n5 = length( epsilonCandidates );

% number of parameter sets
nP = n1 * n2 * n3 * n4 * n5;
paramCell = cell( 1, nP );

idx = 0;
for id1 = 1 : n1
    for id2 = 1 : n2
        for id3 = 1 : n3
            for id4 = 1 : n4
                for id5 = 1 : n5
                    param = [];
                    
                    param.nClusters = nClusters;
                    param.k = kCandidates( id1 );
                    param.beta = betaCandidates( id2 );
                    param.kType = kTypeCandidates( id3 );
                    param.maxiter = maxiterCandidates( id4 );
                    param.epsilon = epsilonCandidates( id5 );
                    
                    idx = idx + 1;
                    paramCell{idx} = param;
                end
            end
        end
    end
end