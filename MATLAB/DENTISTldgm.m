function [P_dentist,T_dentist,Z_tilde] = DENTISTldgm(P, whichIndices, mergedSumstats, normalizePrecision)
%DENTISTldgm computes p-values for LD mismatch or bad summary statistics QC
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
% expected to report a Z score with column name Z_deriv_allele.
% 
% Output arguments:
% P_dentist: DENTIST p-values for LD mismatch for each SNP
% 
% T_dentist: DENTIST \chi^2 statistics
%
% S: cell array containing SNPs that were assigned to each partition


if iscell(P)
    if nargin < 4
        normalizePrecision = false;
    end
    assert(iscell(whichIndices) && iscell(mergedSumstats))
    for ii = 1:numel(P)
        [P_dentist{ii},T_dentist{ii},Z_tilde{ii}] = DENTISTldgm(P{ii}, whichIndices{ii}, mergedSumstats{ii}, normalizePrecision);
    end
else

    % diagonal of R
    Rdiag = ones(length(P),1);
    Rdiag(any(P)) = sqrt(diag(sparseinv(P(any(P),any(P)))));
    
    P = Rdiag .* P .* Rdiag';

    if islogical(whichIndices)
        whichIndices = find(whichIndices);
    end
    noSNPs = length(whichIndices);

    % Randomly partition variants
%     S{1} = sort(randsample(1:noSNPs,floor(noSNPs/2),false));
    S{1} = sort(randsample(1:noSNPs,floor(noSNPs/2),false));
    S{2} = setdiff(1:noSNPs,S{1});

    % submatrix of identity
    for pp = 1:2
        v{pp} = sparse(S{pp},1:length(S{pp}),ones(length(S{pp}),1),noSNPs,length(S{pp}) );
    end
    
    Z = mergedSumstats.Z_deriv_allele;

    % submatrix of inv(P)
    Rit{1} = precisionDivide(P,v{2},whichIndices)' * v{1};
    Rit{2} = Rit{1}';

    T_dentist = zeros(noSNPs,1);
    Z_tilde = T_dentist;

    for pp = 1:2
        qq = 3-pp;

        % Imputed Z scores from the other half of the data
        temp = v{pp} * (precisionMultiply(P,Z(S{pp}),whichIndices(S{pp})));
        Z_tilde(S{qq}) = v{qq} * precisionDivide(P,temp,whichIndices);

        % Compute denominators
        denominators{pp} = 1 - ...
            sum(Rit{pp} .* precisionMultiply(P,Rit{pp}',whichIndices(S{pp}))',2);

        % Compute dentist statistics + p-values
        T_dentist(S{qq}) = (Z(S{qq}) - Z_tilde(S{qq})).^2 ./ denominators{pp};

    end

    % P-values
    P_dentist = chi2cdf(T_dentist,1,'upper');



end