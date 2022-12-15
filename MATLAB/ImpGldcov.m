function [Z_merged, Z_imputed, whichIndices_imputed] = ImpGldcov(R, whichIndices, mergedSumstats, lambda)
%ImpGindiv performs "Bogdan-style" linear imputation of missing summary
% statistics using the method of Pasaniuc et al. 2014 Bioinformatics
%
% Let "0" index the SNPs being imputed, and "1" index the non-missing SNPs.
% The imputation formula is: Z_0 =  R_01 * inv(R_11) * Z_1, where R is the
% LD correlation matrix.
%
% Input arguments:
% R: LDGM correlation matrices, as a number-of-LD-blocks by 1 cell
% array with each cell containing a precision matrix, or as a single
% precision matrix.
%
% whichIndices: output from mergesnplists, encoding which indices
% (rows/columns of the LDGM precision matrices) have corresponding SNPs in
% the summary statistics. Should be a cell array of indices/logicals, or a
% single vector of indices/logicals
%
% mergedSumstats: merged summary statistics tables for each LD block,
% output from mergesnplists. Each table should have height
% equal to the length of corresponding cell of whichIndices. Tables are
% expected to report a Z score with column name Z_deriv_allele. Should be a
% cell array of the same size as P and whichIndices, or a single table.
% 
% parallelize_blocks: optionally, parallelize computation across LD blocks
% (if P, whichIndices and mergedSumstats are cell arrays) (default: 0)
%
% Output arguments:
% Z_imputed: Bogdan-style imputed Z scores.

if iscell(R)
    if nargin < 4
        lambda = 1e-9;
    end
    assert(iscell(whichIndices) && iscell(mergedSumstats))
    Z_imputed = cell(size(R)); whichIndices_imputed = cell(size(R));
    
    
    for ii = 1:numel(R)
        [Z_imputed{ii}, whichIndices_imputed{ii}] = ...
            ImpGldcov(R{ii}, whichIndices{ii}, mergedSumstats{ii}, lambda);
    end

else
    if nargin < 4
        lambda = 1e-9;
    end
    
    % Indices of typed + imputed SNPs; discard SNPs that are missing from
    % the LDGM
    if islogical(whichIndices)
        whichIndices = find(whichIndices);
    end
    noSNPs = length(R);
    allIndices = 1:noSNPs;
    otherIndices = setdiff(allIndices,whichIndices);
    noTyped = length(whichIndices);
    noImputed = length(otherIndices);

    
    % Imputed Z scores
    R = (1-lambda)*R + lambda*eye(noSNPs);
    Z_imputed = R(otherIndices,whichIndices) * (R(whichIndices,whichIndices) \ mergedSumstats.Z_deriv_allele);
    whichIndices_imputed = otherIndices;
    
    Z_merged = zeros(noSNPs,1);
    Z_merged(otherIndices) = Z_imputed;
    Z_merged(whichIndices) = mergedSumstats.Z_deriv_allele;


end