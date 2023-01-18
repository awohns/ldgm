function [Z_merged, Z_imputed, whichIndices_imputed] = ImpGldgm(P, ...
    whichIndices, mergedSumstats, precomputed_ldlChol, parallelize_blocks)
%ImpGldgm performs "Bogdan-style" linear imputation of missing summary
% statistics using the method of Pasaniuc et al. 2014 Bioinformatics
%
% Let "0" index the SNPs being imputed, and "1" index the non-missing SNPs.
% The imputation formula is: Z_0 =  R_01 * inv(R_11) * Z_1, where R is the
% LD correlation matrix.
%
% This is implemented with LDGMs instead using the formula
% Z_0 =  v_0 * (P/P_11) \ (v_1' * (P/P_00) * Z_1)
% where v_0, v_1 are the projection matrices onto the coordinates
% of the missing and non-mising SNPs, respectively.
%
% Input arguments:
% P: LDGM precision matrices, as a number-of-LD-blocks by 1 cell
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

if nargin < 4
    precomputed_ldlChol = false;
end
if precomputed_ldlChol
    error('Not yet implemented');
end

if iscell(P)
    if nargin < 5
        parallelize_blocks = false;
    end
    assert(iscell(whichIndices) && iscell(mergedSumstats))
    Z_imputed = cell(size(P)); 
    whichIndices_imputed = cell(size(P));
    Z_merged = cell(size(P));
    if parallelize_blocks
        parfor ii = 1:numel(P)
            [Z_merged{ii}, Z_imputed{ii}, whichIndices_imputed{ii}] = ...
                ImpGldgm(P{ii}, whichIndices{ii}, mergedSumstats{ii}, precomputed_ldlChol);
        end
    else
        for ii = 1:numel(P)
            [Z_merged{ii}, Z_imputed{ii}, whichIndices_imputed{ii}] = ...
                ImpGldgm(P{ii}, whichIndices{ii}, mergedSumstats{ii}, precomputed_ldlChol);
        end
    end
else
    if istable(mergedSumstats)
        z = mergedSumstats.Z_deriv_allele;
    else
        z = mergedSumstats;
    end

    % Indices of typed + imputed SNPs; discard SNPs that are missing from
    % the LDGM
    if islogical(whichIndices)
        whichIndices = find(whichIndices);
    end
    Pnz = any(P);
    assert(all(Pnz(whichIndices)), 'P should have nonzero diagonal entry for every SNP in summary statistics')
    allIndices = find(Pnz);
    otherIndices = setdiff(allIndices,whichIndices);
    noSNPs = length(P);
    noTyped = length(whichIndices);
    noImputed = length(otherIndices);

    % projection matrices
    V1 = sparse(whichIndices,1:noTyped,...
        ones(noTyped,1),noSNPs,noTyped);
    V0 = sparse(otherIndices,1:noImputed,...
        ones(noImputed,1),noSNPs,noImputed);

    % Imputed Z scores
    temp = V1 * (precisionMultiply(P, z, ...
        whichIndices));
    Z_imputed = V0(allIndices,:)' * precisionDivide(P,temp(allIndices),allIndices);
    whichIndices_imputed = otherIndices;
    
    Z_merged = zeros(noSNPs,1);
    Z_merged(otherIndices) = Z_imputed;
    Z_merged(whichIndices) = z;


end