function [r2, pval] = r2SumstatsPGS(mergedSumstatsValidation, PGSWeightsPerSD, P, whichIndicesSumstats, whichIndicesPGS, sampleSizeValidation, projectOutMissingSNPs, useLdlChol)
%r2SumstatsPGS estimates the PGS accuracy, i.e. the squared correlation between
%Y and X*estimated_beta_perSD, using the precision matrix
%and summary statistics from an external cohort (not used to compute the
%PGS weights).
% 
% Input arguments:
% sumstats: summary statistics, merged with precision matrices using
% mergesumstats
% 
% PGSWeightsPerSD: estimated causal effect sizes, as a cell array of
% the same size, or as a number-of-blocks by 1 cell array (in which case
% the same effect size estimates are used for every population). Should be
% specified in per-SD units.
% 
% P: LDGM precision matrices, as a cell array with each cell containing a 
% precision matrix, for the target population (as opposed to the training
% population). If useLdlChol == true, then P should instead be a cell array
% of Cholesky factors computed using ldlchol.
% 
% whichIndices: output from mergesnplists, encoding which indices
% (rows/columns of the LDGM precision matrices) have corresponding SNPs in
% each summary statistic file respectively. Same size as P.
% 
% AF: allele frequencies, as a number-of-blocks by number-of-populations
% cell array. Can also be empty. If nonempty, trueBeta and estimatedBeta
% should be specified in per-allele units; if empty, they should be in
% per-SD units.
% 
% useLdlChol: whether Cholesky factors have been pre-computed using
% ldlchol(). This option leads to significant speedups. To precompute the
% cholesky factors, apply the following code to each LD block:
%   nz = any(P); % indices not missing from precision matrix
%   A = ldlchol(P(nz,nz)); % Cholesky factor
%   whichIndicesNz = lift(whichIndices, find(nz)); % indices into A
%   trueBeta = trueBeta(nz);
%   estimatedBeta = estimatedBeta(nz);
%   r2 = r2PGS(trueBeta, estimatedBeta, A, whichIndicesNz, AF, true);
% 
% Output arguments:
% r2: the squared correlation between the predictions (this is usually
% larger than the squared correlation between the true and predicted
% weights themselves)
% 
% pval: the p-value for H1: corr(X*PGSWeightsPerSD, y) > 0

[noBlocks,noPopns] = size(P);

assert(noPopns == 1, 'Only supported for a single ancestry group')

if nargin < 6
    sampleSizeValidation = 1;
end
if nargin < 7
    projectOutMissingSNPs = true;
end
if nargin < 8
    useLdlChol = false;
end

if projectOutMissingSNPs
    [PGSweightsProjected, whichIndicesSubset, ~, subsetOfWhichIndicesTarget] = ...
        projectPGSweights(P, whichIndicesPGS,...
        whichIndicesSumstats, PGSWeightsPerSD);
else
    [whichIndicesSubset, whichIndicesSubset, PGSweightsProjected] = deal(cell(size(P)));
    for block = 1:noBlocks
        [whichIndicesSubset{block}, subsetOfWhichIndicesPGS, subsetOfWhichIndicesTarget{block}] = ...
            intersect(whichIndicesPGS{block},whichIndicesSumstats{block},'stable');
        PGSweightsProjected{block} = PGSWeightsPerSD{block}(subsetOfWhichIndicesPGS);
    end
end

% quadratic function
qf = @(beta,P,whichIndices)beta' * precisionDivide(P, beta, whichIndices, useLdlChol);

xtx = zeros(noBlocks,1);
xty = zeros(noBlocks,1);

for block = 1:noBlocks
    alphaHatGWAS = mergedSumstatsValidation{block}.Z_deriv_allele(subsetOfWhichIndicesTarget{block}) / sqrt(sampleSizeValidation);

    xtx(block) = qf(PGSweightsProjected{block},P{block},whichIndicesSubset{block});
    
    xty(block) = alphaHatGWAS' * PGSweightsProjected{block};
  
end

r2 = sum(xty).^2 ./ sum(xtx) ;

pval = normcdf(sqrt(noBlocks) * mean(xty) / std(xty), 'upper');

end

