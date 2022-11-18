function [sumstats, whichIndices, true_beta_perallele, true_beta_perSD] = ...
    simulateSumstats(sampleSize, varargin)
% Simulates summary statistics from specified prior distribution for
% one or more populations.
%
% Required input arguments:
% sampleSize: sample size for each population as a vector
%
%
% % Optional input arguments as name-value pairs:
%
% alleleFrequency: allele frequency for each LD block and each population as
% a number-of-LD blocks by number-of-populations cell array. If specified,
% the effect sizes will be in per-allele rather than per-SD units
%
% precisionMatrices: precision matrix for each LD block and each population
% as a number-of-LD blocks by number-of-populations cell array. Specify
% either this or correlationMatrices.
%
% correlationMatrices: correlation matrix for each LD block and each population
% as a number-of-LD blocks by number-of-populations cell array. Specify
% either this or precisionMatrices
% 
% annotations: annotation matrix for each LD block as a number-of-blocks by
% one cell array. Currently, each element should be a binary matrix with number
% of rows equal to the number of SNPs in each LD block, number of
% columns equal to the number of annotations, and row sums all equal to
% one.
%
% heritability: total heritability for each population, either as a scalar, a vector, or
% a square matrix. If a matrix, it specifies both the heritability for each
% population (along its diagonal) and the genetic correlations (off the
% diagonal, i.e. with corrcof(heritability) == r_pop). If a vector, the
% cross-population genetic correlations will be (almost) 1. If a scalar,
% the cross-population genetic correlations will be (almost) 1, and the
% heritability will be the same for each population.
%
% componentVariance: per-allele effect-size covariance matrix for each
% mixture component, as a number-of-populations by number-of-populations by
% number-of-components array. These will be rescaled to match the total
% heritability.
%
% componentWeight: mixture weight for each heritability component
%
% missingness: fraction of SNPs that are missing, in addition to those for
% which precision matrix already has a zero diagonal element
%
% Output arguments:
% sumstats: cell array of tables, one per LD block, with the following
% columns:
%   Z_deriv_allele: Z scores
%   AF: allele frequencies
%   N: sample size
%
% whichIndices: number-of-blocks by number-of-populations cell array of
% indices for which sumstats are reported. If missingness input is zero
% (default), this is simply the indices for which the precision matrix is
% nonmissing.
%
% componentVariance: per-allele effect-size covariance matrix across
% populations for each variance component. This is normalized to the
% desired heritability for each population.
%
% true_beta_perallele: true per-allele effect size of each variant,
% including those that are missing in the summary statistics
%
% true_beta_perSD: same as true_beta_perallele, but in per-SD units

p=inputParser;

% sample size for each population as a vector
addRequired(p, 'sampleSize', @isnumeric);

% allele frequency for each LD block and each population as a number-of-LD
% blocks by number-of-populations cell array
addParameter(p, 'alleleFrequency', {}, @iscell);

% precision matrix for each LD block and each population as a number-of-LD
% blocks by number-of-populations cell array
addParameter(p, 'precisionMatrices', {}, @iscell);

% correlation matrix for each LD block and each population as a number-of-LD
% blocks by number-of-populations cell array
addParameter(p, 'correlationMatrices', {}, @iscell);

% annotation matrix for each LD block as a number-of-LD
% blocks by 1 cell array. Each SNP should be assigned to exactly one
% annotation, and the annotation matrix should be 0-1 valued with row sums
% equal to one
addParameter(p, 'annotations', {}, @iscell);

% total heritability for each population, either as a scalar, a vector, or
% a square matrix. If a matrix, it specifies both the heritability for each
% population (along its diagonal) and the genetic correlations (off the
% diagonal, i.e. with corrcof(heritability) == r_pop). If a vector, the
% cross-population genetic correlations will be (almost) 1. If a scalar,
% the cross-population genetic correlations will be (almost) 1, and the
% heritability will be the same for each population.
addParameter(p, 'heritability', [], @isnumeric);

% per-allele effect-size covariance matrix for each mixture component,
% as a number-of-populations by number-of-populations by
% number-of-components array. These will be rescaled to match the total
% heritability
addParameter(p, 'componentVariance', [], @isnumeric);

% mixture weight for each heritability component, as a number-of-components
% by number-of-annotations array. If row sums are smaller than one, a null
% component is added with componentVariance equal to zero.
addParameter(p, 'componentWeight', [], @isvector);

% fraction of missing SNPs (missing at random)
addParameter(p, 'missingness', 0, @isscalar);

% turns p.Results.x into just x
parse(p, sampleSize, varargin{:});
fields = fieldnames(p.Results);
for k=1:numel(fields)
    line = sprintf('%s = p.Results.%s;', fields{k}, fields{k});
    eval(line);
end

% input handling
if isempty(precisionMatrices) && isempty(correlationMatrices)
    error('Specify either precision matrices or correlation matrices')
elseif ~isempty(precisionMatrices)
    [noBlocks, noPops] = size(precisionMatrices);
    noSNPs = cellfun(@length,precisionMatrices(:,1));

else
    [noBlocks, noPops] = size(correlationMatrices);
    noSNPs = cellfun(@length,correlationMatrices(:,1));
end

% Sample size for each population
if isscalar(sampleSize)
    sampleSize = sampleSize * ones(1,noPops);
else
    assert(numel(sampleSize) == noPops);
end

% Construct all-ones annotation matrix if needed
if isempty(annotations)
    annotations = arrayfun(@(n){ones(n,1)},noSNPs);
end
noAnnot = size(annotations{1},2);

% Verify that annotation matrix looks OK
for block = 1:noBlocks
    assert(all(size(annotations{block}) == [noSNPs(block), noAnnot]))
    assert(all(sum(annotations{block},2) == 1))
    assert(all(annotations{block}.^2 == annotations{block},'all')) % 0 or 1
end

if noAnnot > 1 && noPops > 1
    error('Multiple populations and multiple annotations not supported')
end

% Heritability for each population
if isscalar(heritability)
    if noPops > 1
        heritability = heritability * ones(1,noPops);
    elseif noAnnot > 1
        error('heritability should be a vector of length noAnnot')
    end
end

% Cross-population heritability matrix
if isvector(heritability) && noPops > 1
    h = sqrt(heritability);
    M = (1-1e-6) * ones(noPops) + 1e-6 * eye(noPops);
    heritability = h' .* M .* h;
end

if isempty(componentVariance)
    assert(isempty(componentWeight),'Specify both componentWeight and componentVariance or neither')
    assert(~isempty(heritability), 'Specify either componentVariance (and componentWeight), or heritability')
    if noAnnot == 1
        componentVariance = heritability;
    else
        componentVariance = 1;
    end
    componentWeight = 1;
end

if noPops == 1
    if isvector(componentVariance)
        componentVariance = reshape(componentVariance,1,1,length(componentVariance));
    end

end

% componentWeight should specify a convex combination of components
assert(all(componentWeight<=1) & all(componentWeight>=0))
assert(sum(componentWeight)<=1)
noCpts = length(componentWeight);

% if componentWeight sums to <1, null component is added
if sum(componentWeight) < 1
    componentWeight(end+1) = 1-sum(componentWeight);
    componentVariance(:,:,end+1) = zeros(noPops);
end

% Simulate summary statistics for each LD block
for block = 1:noBlocks
    whichCpt = randsample(1:length(componentWeight),noSNPs(block),true,componentWeight);
    beta = zeros(noSNPs(block),noPops);
    for cpt = 1:noCpts
        beta(whichCpt == cpt,:) = mvnrnd(zeros(1,noPops), componentVariance(:,:,cpt), sum(whichCpt == cpt));
    end

    for pop = 1:noPops
        true_beta_perallele{block,pop} = beta(:,pop);

        % standardized effect sizes (sd of Y per sd of X)
        if ~isempty(alleleFrequency)
            true_beta_perSD{block,pop} = beta(:,pop) .* ...
                sqrt(2 * alleleFrequency{block,pop} .* (1 - alleleFrequency{block,pop}));
        else
            true_beta_perSD{block,pop} = beta(:,pop);
        end
    end

end


% Normalize effect sizes so they add up to h2 within each annotation
if ~isempty(heritability)
    if noAnnot > 1

        % annotEffectSum: sum(beta.^2) within each annotation
        for aa = 1:noAnnot
            for block = 1:noBlocks
                for block = 1:noBlocks
                    incl = annotations{block}(:,aa)==1;
                    annotEffectSum(block,aa) = sum(true_beta_perSD{block}(incl).^2);
                end
            end
        end
        annotEffectSum = sum(annotEffectSum);
        
        % normalize effect sizes
        normalizer = sqrt(heritability./annotEffectSum);
        for aa = 1:noAnnot
            for block = 1:noBlocks
                incl = annotations{block}(:,aa)==1;
                true_beta_perallele{block}(incl) = true_beta_perallele{block}(incl) * normalizer(aa);
                true_beta_perSD{block}(incl) = true_beta_perSD{block}(incl) * normalizer(aa);
            end
        end

        % If no annotations (but possibly multiple popns), normalize effect sizes
        % so they add up to h2 within each population
    else
        for pop = 1:noPops
            normalizer(pop) = sqrt(heritability(pop,pop)/sum(cellfun(@(x)sum(x.^2),true_beta_perSD(:,pop))));
            true_beta_perallele(:,pop) = cellfun(@(b){b*normalizer(pop)}, true_beta_perallele(:,pop));
            true_beta_perSD(:,pop) = cellfun(@(b){b*normalizer(pop)}, true_beta_perSD(:,pop));
        end
    end

end

% Sample summary statistics
Z = cell(noBlocks,noPops);
whichIndices = Z;
for block = 1:noBlocks
    for pop = 1:noPops

        % Simulate sumstats using precision matrices
        if ~isempty(precisionMatrices)
            % SNPs not missing in precision matrix
            incl = diag(precisionMatrices{block,pop}) ~= 0;
            whichIndices{block,pop} = find(incl);
            Z{block,pop} = precisionMatrices{block,pop}(incl,incl)...
                \ true_beta_perSD{block,pop}(incl) * sqrt(sampleSize(pop))...
                + chol(precisionMatrices{block,pop}(incl,incl)) \...
                randn(sum(incl),1);

            % Simulate using correlation matrices
        else
            incl = diag(correlationMatrices{block,pop}) ~= 0;
            whichIndices{block,pop} = find(incl);
            Z{block,pop} = correlationMatrices{block,pop}(incl,incl) *...
                true_beta_perSD{block,pop}(incl) * sqrt(sampleSize(pop))...
                + chol(correlationMatrices{block,pop}(incl,incl))' *...
                randn(sum(incl),1);
        end

        % Additional SNPs missing at random
        if missingness > 0
            incl = rand(sum(incl),1) > missingness;
            Z{block,pop} = Z{block,pop}(incl);
            whichIndices{block,pop} = whichIndices{block,pop}(incl);
        end
    end
end

% Store output in a table
sumstats = cell(noBlocks,noPops);
for block = 1:noBlocks
    for pop = 1:noPops
        sumstats{block,pop} = table('size',[length(Z{block,pop}),0]);
        sumstats{block,pop}.Z_deriv_allele = Z{block,pop};
        sumstats{block,pop}.N = sampleSize(pop) * ones(size(Z{block,pop}));
        if ~isempty(alleleFrequency)
            sumstats{block,pop}.AF_deriv_allele = ...
                binornd( 2*sumstats{block,pop}.N, ...
                alleleFrequency{block,pop}(whichIndices{block,pop}) ) ./ ...
                (2*sumstats{block,pop}.N);
        end
    end
end


end