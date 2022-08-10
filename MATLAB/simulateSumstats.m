function [sumstats, whichIndices, componentVariance, true_beta_perallele, true_beta_perSD] = ...
    simulateSumstats(sampleSize, alleleFrequency, varargin)
% Simulates summary statistics from specified prior distribution for
% one or more populations.
% 
% Required input arguments:
% sampleSize: sample size for each population as a vector
% 
% alleleFrequency: allele frequency for each LD block and each population as
% a number-of-LD blocks by number-of-populations cell array
% 
% Optional input arguments as name-value pairs:
% 
% precisionMatrices: precision matrix for each LD block and each population
% as a number-of-LD blocks by number-of-populations cell array. Specify
% either this or correlationMatrices.
% 
% correlationMatrices: correlation matrix for each LD block and each population
% as a number-of-LD blocks by number-of-populations cell array. Specify
% either this or precisionMatrices
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
addRequired(p, 'alleleFrequency', @iscell);

% precision matrix for each LD block and each population as a number-of-LD
% blocks by number-of-populations cell array
addParameter(p, 'precisionMatrices', {}, @iscell);

% correlation matrix for each LD block and each population as a number-of-LD
% blocks by number-of-populations cell array
addParameter(p, 'correlationMatrices', {}, @iscell);

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

% mixture weight for each heritability component
addParameter(p, 'componentWeight', [], @isvector);

% fraction of missing SNPs
addParameter(p, 'missingness', 0, @isscalar);


% turns p.Results.x into just x
parse(p, sampleSize, alleleFrequency, varargin{:});
fields = fieldnames(p.Results);
for k=1:numel(fields)
    line = sprintf('%s = p.Results.%s;', fields{k}, fields{k});
    eval(line);
end

% Allele frequency
[noBlocks, noPops] = size(alleleFrequency);
if isscalar(sampleSize)
    sampleSize = sampleSize * ones(1,noPops);
else
    assert(numel(sampleSize) == noPops);
end
noSNPs = cellfun(@length,alleleFrequency);

if isscalar(heritability)
    heritability = heritability * ones(1,noPops);
end
if isvector(heritability)
    h = sqrt(heritability);
    M = (1-1e-6) * ones(noPops) + 1e-6 * eye(noPops);
    heritability = h' .* M .* h;
end

if ~isempty(precisionMatrices)
    assert(all(size(alleleFrequency) == size(precisionMatrices)),...
        'Input cell arrays must have size number of LD blocks by number of populations')
    assert(isempty(correlationMatrices),'Specify precisionMatrices or correlationMatrices but not both')
elseif ~isempty(correlationMatrices)
    assert(all(size(alleleFrequency) == size(correlationMatrices)),...
        'Input cell arrays must have size number of LD blocks by number of populations')
end
if isempty(precisionMatrices) && isempty(correlationMatrices)
    correlationMatrices = arrayfun(@speye, noSNPs, 'uniformOutput', false);
end

if isempty(componentVariance)
    assert(isempty(componentWeight),'Specify both componentWeight and componentVariance or neither')
    assert(~isempty(heritability), 'Specify either componentVariance (and componentWeight) or heritability')
    componentVariance = heritability;
    componentWeight = 1;
end

if noPops == 1
    if ndims(componentVariance) == 1
        componentVariance = reshape(componentVariance,1,1,length(componentVariance));
    end
    
end
% assert(all(size(componentVariance,1:3) == [noPops, noPops, length(componentWeight)]))
assert(all(componentWeight<=1) & all(componentWeight>=0))
assert(sum(componentWeight)<=1)
noCpts = length(componentWeight);

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
        true_beta_perSD{block,pop} = beta(:,pop) .* ...
            sqrt(2 * alleleFrequency{block,pop} .* (1 - alleleFrequency{block,pop}));
    end
end

% Normalize effect sizes so they add up to h2
if ~isempty(heritability)
    for pop = 1:noPops
        normalizer(pop) = sqrt(heritability(pop,pop)/sum(cellfun(@(x)sum(x.^2),true_beta_perSD(:,pop))));
        true_beta_perallele(:,pop) = cellfun(@(b){b*normalizer(pop)}, true_beta_perallele(:,pop));
        true_beta_perSD(:,pop) = cellfun(@(b){b*normalizer(pop)}, true_beta_perSD(:,pop));
        
    end
    % normalize variance component matrices as well
    for cpt = 1:noCpts
        componentVariance(:,:,cpt) = componentVariance(:,:,cpt) .* (normalizer .* normalizer');
    end
end

% Sample summary statistics
Z = cell(size(alleleFrequency));
whichIndices = Z;
for block = 1:noBlocks
    for pop = 1:noPops
        
        
        if ~isempty(precisionMatrices)
            % SNPs not missing in precision matrix
            incl = diag(precisionMatrices{block,pop}) ~= 0;
            whichIndices{block,pop} = find(incl);
            Z{block,pop} = precisionMatrices{block,pop}(incl,incl)...
                \ true_beta_perSD{block,pop}(incl) * sqrt(sampleSize(pop))...
                + chol(precisionMatrices{block,pop}(incl,incl)) \...
                randn(sum(incl),1);
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
        sumstats{block,pop}.AF = alleleFrequency{block,pop}(whichIndices{block,pop});
        sumstats{block,pop}.N = sampleSize(pop) * ones(size(Z{block,pop}));
    end
end


end