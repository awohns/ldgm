function [PGSweightsProjected, whichIndicesSubset, subsetOfWhichIndicesInitial, subsetOfWhichIndicesTarget] = projectPGSweights(P, whichIndicesInitial,...
    whichIndicesTarget, PGSweightsInitial_perSD)
%PGSweightsProjected calculates PGS weights for a subset of the initial
%SNPs, which is useful when some of the SNPs in the training data are
%missing from the target dataset.
% 
% Input arguments:
% P: LDGM precision matrices, as a number-of-LD-blocks by 1 cell
% array with each cell containing a precision matrix.
%
% whichIndicesInitial: output from mergesnplists, encoding which indices
% (rows/columns of the precision matrices) have corresponding SNPs in
% the initial PGS weights.
% 
% whichIndicesSubset: output from mergesnplists, encoding which indices
% (rows/columns of the precision matrices) have corresponding SNPs in
% the target dataset. If not a subset of whichIndicesInitial, extra SNPs
% are given a PGS weight of 0.
% 
% PGSweightsInitial: initial PGS weights.

noBlocks = numel(P);
[whichIndicesSubset, subsetOfWhichIndicesInitial, subsetOfWhichIndicesTarget]  = deal(cell(size(P)));
for block = 1:noBlocks
    [whichIndicesSubset{block}, subsetOfWhichIndicesInitial{block}, subsetOfWhichIndicesTarget{block}] = ...
        intersect(whichIndicesInitial{block}, whichIndicesTarget{block}, 'stable');
end

x = precisionDivide(P, PGSweightsInitial_perSD,  whichIndicesInitial);

x = cellfun(@(x,j){x(j)},x,subsetOfWhichIndicesInitial);
PGSweightsProjected = precisionMultiply(P, x, whichIndicesSubset);

end

