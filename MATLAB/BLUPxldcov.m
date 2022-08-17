function [betaExpectationPerAllele, betaExpectationPerSD] =...
    BLUPxldcov(R, whichIndices, mergedSumstats, betaCov, sampleSize, alpha_param)
% BLUPx computes the cross-popn best linear unbiased predictor, E(beta|GWAS, gaussian
% prior).
%
% Input arguments:
% R: LD correlation matrices, as a number-of-LD-blocks by number-of-popns cell
% array with each cell containing a correlation matrix; should be the same
% size across populations. For SNPs that are missing from the summary
% statistics in some of the populations, R can be padded with zeros in the
% appropriate entries
%
% whichIndices: output from mergesnplists, encoding which indices
% (rows/columns of the correlation matrices) have corresponding SNPs in
% each population. Should be a number-of-blocks by number-of-popns cell array.
%
% mergedSumstats: merged summary statistics tables for each LD block -
% population pair, output from mergesnplists. Each table should have height
% equal to the length of corresponding cell of whichIndices. Tables are
% expected to report either a Z score (column name Z) or an effect size and
% standard error (beta, se). Optionally, if they report an allele frequency
% (AF) and alpha_param is specified, an AF-dependent prior will be used.
%
% sampleSize: GWAS sample size for each population, as a vector.
%
% betaCov: covariance matrix for the per-allele effect size of a SNP across
% popns, which is assumed to be i.i.d.
%
% Output arguments:
% betaExpectationPerSD: expected value of beta for each popn, in per-s.d.
% (standardized) units,  as a cell array of the same size as input
% arguments. Values are reported for every index in the precision matrices,
% but when some of these indices are missing summary statistics or
% precision matrix entries, beta will be zero in those positions.
% betaExpectationPerAllele: same as betaExpectationPerSD, but in per-allele
% instead of per-s.d. units. This can only be reported if AF is specified
% in the sumstats tables.

assert(all(size(R) == size(whichIndices)), 'Input cell arrays should agree in size')
assert(all(size(R) == size(mergedSumstats)), 'Input cell arrays should agree in size')
has_Zscore = cellfun(@(T)any(strcmp(T.Properties.VariableNames,'Z_deriv_allele')),mergedSumstats);
assert(all(has_Zscore),...
    'mergedSumstats should have a field called Z_deriv_allele. Please use the mergesnplists function.')

[noBlocks, noPopns] = size(R);
assert(all(size(betaCov) == [noPopns, noPopns]), 'betaCov should be a square matrix of size noPopns')
assert(all(size(sampleSize) == [1 noPopns]), 'sampleSize should be a row vector of size noPopns')

% Get rid of empty rows/columns of precision matrices
nonempty = cell(size(R));
for block = 1:noBlocks
    for popn = 1:noPopns
        nonempty{block,popn} = find(diag(R{block,popn}...
            (whichIndices{block,popn},whichIndices{block,popn})));
        mergedSumstats{block,popn} = mergedSumstats{block,popn}(nonempty{block,popn},:);
        whichIndices{block,popn} = whichIndices{block,popn}(nonempty{block,popn});
    end
end
noSNPs = cellfun(@length, whichIndices);

% Column names of sumstats files (assumes first one is representative)
column_names = mergedSumstats{1}.Properties.VariableNames;

% Marginal effect-size estimates (i.e., sample correlations)
alphaHat = cellfun(@(T){T.Z_deriv_allele},mergedSumstats);
nn_cell = repmat(num2cell(sampleSize),noBlocks,1);
alphaHat = cellfun(@(v,n){v/sqrt(n)},alphaHat,nn_cell);

% effect-size s.d. for each SNP in each population, in standardized
% (per-sd-of-genotype) units. alpha_param is the AF-dependent architecture parameter
% of e.g. Schoech et al. 2019. For a SNP that is missing in a population,
% its effect-size s.d. is zero.
AF_col = strcmpi(column_names, 'AF');
if ~any(AF_col)
    SD = cellfun(@double,whichIndices,'uniformoutput',0);
else
    if nargin < 6
        % default no AF-dependent architecture (constant per-allele effect
        % size variance)
        alpha_param = 0;
    end
    % assign sqrt(2pq) (if alpha_param==0) to nonmissing SNPs
    SD = cell(size(R));
    for block = 1:noBlocks
        for popn = 1:noPopns
            AF = table2array(mergedSumstats{block,popn}(:,AF_col));
            assert(all(min(AF,1-AF) < 1));
            SD{block,popn} = assignto((2*AF.*(1-AF)) .^ ((alpha_param+1)/2),...
                whichIndices{block,popn});
        end
    end
end

% Concatenated effect-size estimates and correlation matrices across popns
nnAlphaHatCat = cell(noBlocks,1);
nnRCat = cell(noBlocks,1);
whichIndicesCat = cell(noBlocks,1);
for block = 1:noBlocks
    nnAlphaHat = cellfun(@times, num2cell(sampleSize), alphaHat(block,:),'UniformOutput',false);
    nnAlphaHatCat{block} = vertcat(nnAlphaHat{:});
    nnRblocks = cellfun(@(x,y)y * x, num2cell(sampleSize), R(block,:),'UniformOutput',false);
    nnRCat{block} = blkdiag(nnRblocks{:});
    temp = cellfun(@(ii,X){unfind(ii,length(X))},whichIndices(block,:),R);
    whichIndicesCat{block} = find(vertcat(temp{:}));
end

% Precision matrix of effect sizes for each concatenated block
betaPrecision = inv(betaCov); % per-allele units
SigmaInv = cell(noBlocks,1);
for block = 1:noBlocks
    S = cell(noPopns); % blocks of covariance matrix
    for ii = 1:noPopns
        for jj = 1:noPopns
            [~, i1, i2] = intersect(whichIndices{block,ii}, whichIndices{block,jj});
            S{ii,jj} = sparse(i1,i2,...
                betaPrecision(ii,jj)./(SD{block,ii}(i1).*SD{block,jj}(i2)),...
                noSNPs(block,ii),noSNPs(block,jj));
        end
    end
    SigmaInv{block} = zeros(size(nnRCat{block}));
    SigmaInv{block}(whichIndicesCat{block},whichIndicesCat{block}) = cell2mat(S);
end

% E(beta|data)
betaExpectationCat = cell(noBlocks,1);
for block = 1:noBlocks
    betaExpectationCat{block} = (nnRCat{block}(whichIndicesCat{block},whichIndicesCat{block})...
        + SigmaInv{block}(whichIndicesCat{block},whichIndicesCat{block}))...
        \ nnAlphaHatCat{block};
end

% Un-concatenate populations
betaExpectationPerSD = cell(size(alphaHat));
for block = 1:noBlocks
    betaExpectationPerSD(block,:) = mat2cell(betaExpectationCat{block}, noSNPs(block,:));
end

% per-allele effect sizes
betaExpectationPerAllele = cellfun( @(beta,sd)beta./sd, ...
    betaExpectationPerSD, SD,'UniformOutput',false);



end

