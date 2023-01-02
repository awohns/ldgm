function [P_dentist,T_dentist,Z_tilde] = DENTISTldcov(R, whichIndices, mergedSumstats, normalizeCorrelation, partition )
%DENTISTldcov computes p-values for LD mismatch or bad summary statistics QC
%using LDGM precision matrices and GWAS summary statistics
%
% Re-implements the DENTIST method of Chen et al. 2021 Nat Commun,
% "Improved analyses of GWAS summary statistics by reducing data
% heterogeneity and errors"
% 
% However, this implementation is simplified: it only implements the
% formula in equation (1), without performing the iterative SNP-removal
% procedure. This means that for loci containing bad SNPs, it might inflate
% 
% 
% Input arguments:
% R: LDGM correlation matrices, as a number-of-LD-blocks by 1 cell
% array with each cell containing a correlation matrix, or as a single
% precision matrix.
%
% whichIndices: output from mergesnplists, encoding which indices
% (rows/columns of the correlation matrices) have corresponding SNPs in
% the summary statistics. Should be a cell array of indices/logicals, or a
% single vector of indices/logicals
%
% mergedSumstats: merged summary statistics tables for each LD block, 
% output from mergesnplists. Each table should have height
% equal to the length of corresponding cell of whichIndices. Tables are
% expected to report a Z score with column name Z_deriv_allele.
% 
% Output arguments:
% P_dentist: DENTIST p-values for LD mismatch for each SNP
% 
% T_dentist: DENTIST \chi^2 statistics
%
% S: cell array containing SNPs that were assigned to each partition

if nargin < 4
    normalizeCorrelation = false;
end
if iscell(R)
    assert(iscell(whichIndices) && iscell(mergedSumstats))
    for ii = 1:numel(R)
        [P_dentist{ii},T_dentist{ii},Z_tilde{ii}] = DENTISTldcov(R{ii}, whichIndices{ii}, mergedSumstats{ii}, normalizeCorrelation);
    end
else

    % diagonal of R -> 1
    if normalizeCorrelation
        R = corrcov(R);
    end

    if islogical(whichIndices)
        whichIndices = find(whichIndices);
    end
    noSNPs = length(whichIndices);

    % Randomly partition variants
    if nargin < 5
        partition{1} = sort(randsample(1:noSNPs,floor(noSNPs/2),false));
        partition{2} = setdiff(1:noSNPs,partition{1});
    else
        assert(iscell(partition) & length(partition) == 2)
    end

    % submatrix of identity
    for pp = 1:2
        v{pp} = sparse(partition{pp},1:length(partition{pp}),ones(length(partition{pp}),1),noSNPs,length(partition{pp}) );
    end
    
    Z = mergedSumstats.Z_deriv_allele;

    T_dentist = zeros(noSNPs,1);
    Z_tilde = T_dentist;

    for pp = 1:2
        qq = 3-pp;

        % Imputed Z scores from the other half of the data
        temp = R(whichIndices(partition{pp}),whichIndices(partition{pp})) \ ...
            Z(partition{pp});
        
        Z_tilde(partition{qq}) = R(whichIndices(partition{qq}),whichIndices(partition{pp})) * temp;

        % Compute denominators
        denominators{pp} = 1 - ...
            sum(R(whichIndices(partition{pp}),whichIndices(partition{qq}))' .* ...
            (R(whichIndices(partition{pp}),whichIndices(partition{pp})) \ ...
            R(whichIndices(partition{pp}),whichIndices(partition{qq})))',2);

        % Compute dentist statistics + p-values
        T_dentist(partition{qq}) = (Z(partition{qq}) - Z_tilde(partition{qq})).^2 ./ denominators{pp};

    end

    % P-values
    P_dentist = chi2cdf(T_dentist,1,'upper');



end