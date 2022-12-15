function [betaExpectationPerAllele, betaExpectationPerSD] =...
    BLUPxldgm(P, whichIndices, mergedSumstats, betaCov, sampleSize, alpha_param)
% BLUPx computes the best linear unbiased predictor, E(beta|GWAS), under
% a Gaussian model. It inputs data from one or more ancestry groups.
%
% Input arguments:
% P: LDGM precision matrices, as a number-of-LD-blocks by number-of-popns cell
% array with each cell containing a precision matrix; should be the same
% size across populations. For SNPs that are missing from the summary
% statistics in some of the populations, P can be padded with zeros in the
% appropriate entries
%
% whichIndices: output from mergesnplists, encoding which indices
% (rows/columns of the LDGM precision matrices) have corresponding SNPs in
% each population. Should be a number-of-blocks by number-of-popns cell array.
% If no SNPs are missinig (e.g., in simulations), specify cells as
% whichIndices{jj} = true(length(P{jj}),1);
%
% mergedSumstats: merged summary statistics tables for each LD block -
% population pair, output from mergesnplists. Each table should have height
% equal to the length of corresponding cell of whichIndices. Tables are
% expected to report a Z score with column name Z_deriv_allele.
% Optionally, if they report an allele frequency with column name AF,
% then an AF-dependent prior can be used
%
% betaCov: covariance matrix for the per-allele effect size of a SNP across
% popns, which is assumed to be i.i.d.
%
% sampleSize (optional): GWAS sample size for each population, as a vector.
% If not specified, BLUPxldgm will look for a column of the summary statistics
% file called 'N' and use that if it is found.
%
% alpha_param (optional): frequency-dependent architecture parameter;
% models per-allele effect size variance of each SNP as
% var(beta) == betaCov * heterozygosity ^ alpha. Default: 0.
%
% Output arguments:
% betaExpectationPerSD: expected value of beta for each popn, in per-s.d.
% (standardized) units,  as a cell array of the same size as input
% arguments. Values are reported for every index in the precision matrices,
% but when some of these indices are missing summary statistics or
% precision matrix entries, beta will be zero in those positions.
%
% betaExpectationPerAllele: same as betaExpectationPerSD, but in per-allele
% instead of per-s.d. units. This can only be reported if AF is specified
% in the sumstats tables.

assert(all(size(P) == size(whichIndices)), 'Input cell arrays should agree in size')
assert(all(size(P) == size(mergedSumstats)), 'Input cell arrays should agree in size')
has_Zscore = cellfun(@(T)any(strcmp(T.Properties.VariableNames,'Z_deriv_allele')),mergedSumstats);
assert(all(has_Zscore,'all'),...
    'mergedSumstats should have a field called Z_deriv_allele. Please use the mergesnplists function.')

[noBlocks, noPopns] = size(P);
assert(all(size(betaCov) == [noPopns, noPopns]), 'betaCov should be a square matrix of size noPopns')

% Get sample size from sumstats if necessary
if nargin < 5
    for ii = 1:noBlocks
        for jj = 1:noPopns
            assert(any(strcmpi(mergedSumstats{ii,jj}.Properties.VariableNames,'N')), 'Please specify sample size, or alternatively include a column named N in the merged sumstats');
            sampleSize(ii,jj) = mean(mergedSumstats{ii,jj}.N);
        end
    end
    sampleSize = mean(sampleSize);
end
assert(all(size(sampleSize) == [1 noPopns]), 'sampleSize should be a row vector of size noPopns')

% Convert whichIndices to indices if needed
for block = 1:noBlocks
    for popn = 1:noPopns
        if islogical(whichIndices{block,popn})
            whichIndices{block,popn} = find(whichIndices{block,popn});
        end
    end
end
noSNPs = cellfun(@length, whichIndices);

% Column names of sumstats files (assumes first one is representative)
column_names = mergedSumstats{1}.Properties.VariableNames;

% Marginal effect-size estimates (i.e., sample correlations)
alphaHat = cellfun(@(T){T.Z_deriv_allele},mergedSumstats);
nn_cell = repmat(num2cell(sampleSize),noBlocks,1);
alphaHat = cellfun(@(v,n){v/sqrt(n)},alphaHat,nn_cell);

% Causal effect-size estimates
betaHat = precisionMultiply(P,alphaHat,whichIndices);

% effect-size s.d. for each SNP in each population, in standardized
% (per-sd-of-genotype) units. alpha_param is the AF-dependent architecture parameter
% of e.g. Schoech et al. 2019. For a SNP that is missing in a population,
% its effect-size s.d. is zero.
AF_col = contains(column_names, 'AF', 'IgnoreCase',true);

if ~any(AF_col)
    SD = arrayfun(@(m)ones(m,1),noSNPs,'uniformoutput',0);
else
    if sum(AF_col) > 1
        warning('Multiple columns found with name containing AF; choosing %s',...
            column_names{find(AF_col,1)})
        AF_col = find(AF_col,1);
    end

    if nargin < 6
        % default no AF-dependent architecture (constant per-allele effect
        % size variance)
        alpha_param = 0;
    end
    
    % assign sqrt(2pq) (if alpha_param==0) to nonmissing SNPs
    SD = cell(size(P));
    for block = 1:noBlocks
        for popn = 1:noPopns
            AF = table2array(mergedSumstats{block,popn}(:,AF_col));
            assert(all(min(AF,1-AF) < 1));
            SD{block,popn} = (2*AF.*(1-AF)) .^ ((alpha_param+1)/2);
        end
    end
end

% Concatenated effect-size estimates and precision matrices across popns
betaHatCat = cell(noBlocks,1);
PCat = cell(noBlocks,1);
whichIndicesCat = cell(noBlocks,1);
for block = 1:noBlocks
    betaHatCat{block} = vertcat(betaHat{block, :});
    Pblocks = cellfun(@(x,y)y/x, num2cell(sampleSize), P(block,:),'UniformOutput',false);
    PCat{block} = blkdiag(Pblocks{:});
    temp = cellfun(@(ii,X){unfind(ii,length(X))},whichIndices(block,:),P(block,:));
    whichIndicesCat{block} = find(vertcat(temp{:}));
end

% Covariance matrix of effect sizes for each concatenated block
Sigma = cell(noBlocks,1);
for block = 1:noBlocks
    S = cell(noPopns); % blocks of covariance matrix
    for ii = 1:noPopns
        for jj = 1:noPopns
            [~, i1, i2] = intersect(whichIndices{block,ii}, whichIndices{block,jj});
            S{ii,jj} = sparse(i1,i2,...
                betaCov(ii,jj).*SD{block,ii}(i1).*SD{block,jj}(i2),...
                noSNPs(block,ii),noSNPs(block,jj));
        end
    end
    Sigma{block} = sparse(length(PCat{block}),length(PCat{block}));
    Sigma{block}(whichIndicesCat{block},whichIndicesCat{block}) = cell2mat(S);
end

% E(beta|data)
x = precisionDivide(cellfun(@(p,s){p+s},PCat,Sigma),...
    betaHatCat, whichIndicesCat);
betaExpectationCat = cellfun(@(s,v,idx)full(s(idx,idx)*v),Sigma,x,whichIndicesCat,...
    'UniformOutput',false);%precisionMultiply(Sigma, x, whichSNPsCat);

betaExpectationPerSD = cell(size(alphaHat));
betaExpectationPerAllele = betaExpectationPerSD;
for block = 1:noBlocks
    % Un-concatenate populations
    betaExpectationPerSD(block,:) = mat2cell(betaExpectationCat{block}, noSNPs(block,:));

    % Assign to whichIndices
    for popn = 1:noPopns
        % per-allele effect sizes
        betaExpectationPerAllele{block,popn} = assignto(...
            betaExpectationPerSD{block,popn} ./ SD{block,popn}, whichIndices{block,popn},...
            length(P{block,popn}));

        betaExpectationPerSD{block,popn} = assignto(...
            betaExpectationPerSD{block,popn}, whichIndices{block,popn},...
            length(P{block,popn}));
    end
end




end

