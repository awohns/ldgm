function h2 = HESSldgm(P, whichIndices, mergedSumstats, sampleSize)
%HESSldgm is an implementation of the heritability estimation with summary
%statistics (HESS) method of H. Shi et al. 2016 AJHG. 

noBlocks = length(mergedSumstats);
if nargin < 4
    for ii = 1:noBlocks
        assert(any(strcmpi(mergedSumstats{ii}.Properties.VariableNames,'N')), 'Please specify sample size, or alternatively include a column named N in the merged sumstats');
        sampleSize(ii) = mean(mergedSumstats{ii}.N);
    end
elseif isscalar(sampleSize)
    sampleSize = sampleSize * ones(noBlocks,1);
else
    assert(numel(sampleSize) == noBlocks);
end

h2 = zeros(noBlocks,1);
for ii = 1:noBlocks
    beta = precisionMultiply(P{ii},mergedSumstats{ii}.Z_deriv_allele,whichIndices{ii});
    noSNPs = length(beta);
    h2(ii) = (beta' * mergedSumstats{ii}.Z_deriv_allele - noSNPs) / ...
        (sampleSize(ii) - noSNPs);
end
end

