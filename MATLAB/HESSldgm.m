function h2 = HESSldgm(P, whichIndices, mergedSumstats, sampleSize)
%HESSldgm is an implementation of the heritability estimation with summary
%statistics (HESS) method of H. Shi et al. 2016 AJHG. 
% Input arguments:
% P: LDGM precision matrices, as a cell array with each cell containing a 
% precision matrix
% 
% whichIndices: output from mergesnplists, encoding which indices
% (rows/columns of the LDGM precision matrices) have corresponding SNPs in
% each summary statistic file respectively. Same size as P.
% 
% mergedSumstats: merged summary statistics tables for each LD block, 
% the output from mergesnplists. Each table should have height
% equal to the length of corresponding cell of whichIndices. Tables are
% expected to report a Z score with name Z_deriv_allele. Same size as P.
% Alternatively, specify the Z scores directly as a cell array of column 
% vectors. 
% 
% sampleSize (optional): GWAS sample size for each population, as either a
% scalar or a matrix of size equal to size(P). If not
% specified, BLUPxldgm will look for a column of the summary statistics
% file called 'N' and use that if it is found.
% 
% Output arguments:
% h2: heritability estimate for each LD block, same size as P.

noBlocks = numel(mergedSumstats);

% Look up sample size if needed
if nargin < 4
    sampleSize = zeros(size(P));
    for ii = 1:noBlocks
        assert(istable(mergedSumstats{ii}),'Please specify sample size')
        assert(any(strcmpi(mergedSumstats{ii}.Properties.VariableNames,'N')),...
            'Please specify sample size, or alternatively include a column named N in the merged sumstats');
        sampleSize(ii) = mean(mergedSumstats{ii}.N);
    end
elseif isscalar(sampleSize)
    sampleSize = sampleSize * ones(size(P));
else
    assert(all(size(sampleSize) == size(P)));
end

% Look up Z scores if needed
Z = cell(size(P));
for ii = 1:noBlocks
    if istable(mergedSumstats{ii})
        assert(any(strcmpi(mergedSumstats{ii}.Properties.VariableNames,'N')),...
            'Please include a column named Z_deriv_allele in the merged sumstats');
        Z{ii} = mergedSumstats{ii}.Z_deriv_allele;
    else
        assert(iscolumn(mergedSumstats{ii}),'Please specify Z scores as a column vector');
        Z{ii} = mergedSumstats{ii};
    end
end


h2 = zeros(size(P));
for ii = 1:noBlocks
    beta = precisionMultiply(P{ii},Z{ii},whichIndices{ii});
    noSNPs = length(beta);
    h2(ii) = (beta' * Z{ii} - noSNPs) / ...
        (sampleSize(ii) - noSNPs);
    if beta'*Z{ii} > Z{ii}'*Z{ii}
        warning('Suspiciously large h2 estimate suggests something is wrong in block %d', ii)
    end
end
end

