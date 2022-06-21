function [alphaHat, beta_perallele, beta_persd, componentVariance] = simulate_sumstats_populations(sampleSize, alleleFrequency, varargin)
% Simulates summary statistics from specified prior distribution for
% multiple populations. 
% Input arguments: see below.
% Output arguments: 
% alphaHat: sample correlation between the trait and each SNP in each
% population, as a number-of-LD-blocks by number-of-populations cell array
% beta_perallele: per-allele effect size of each SNP in each population
% beta_persd: per-s.d. effect size of each SNP in each population
% componentVariance: per-allele effect-size covariance matrix across
% populations for each variance component. This is normalized to the
% desired heritability for each population.

p=inputParser;

% sample size for each population
addRequired(p, 'sampleSize', @isnumeric);

% allele frequency for each LD block and each population
addRequired(p, 'alleleFrequency', @iscell);

% SNP identifiers (eg, RSIDs) for each LD block and each population, as a
% number-of-blocks by number-of-populations cell array. These can be
% redundant across blocks; they are joined across populations.
addParameter(p, 'SNPs', {}, @iscell);

% precision matrix for each LD block and each population
addParameter(p, 'precisionMatrices', {}, @iscell);

% correlation matrix for each LD block and each population
addParameter(p, 'correlationMatrices', {}, @iscell);

% total heritability for each population
addParameter(p, 'heritability', [], @isscalar);

% per-allele effect-size covariance matrix for each heritability component,
% as a number-of-populations by number-of-populations by
% number-of-components array. These will be rescaled to match the total
% heritability
addParameter(p, 'componentVariance', [], @isnumeric);

% mixture weight for each heritability component
addParameter(p, 'componentWeight', [], @isvector);

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

if isempty(SNPs)
    assert(all(min(noSNPs,[],2) == max(noSNPs,[],2)),...
        'If SNP identifiers not specified, the number of SNPs in each LD block must be the same across populations');
    SNPs = cellfun(@(n){1:n},noSNPs);
else
    assert(all(cellfun(@iscolumn,SNPs), 'all'), 'Specify SNPs as a cell array of column vectors')
    for ii = 1:noBlocks
        allSNPs = unique(vertcat(SNPs{ii,:}));
        for pop = 1:noPops
            [~, ~, SNPs{ii,pop}] = intersect(SNPs{ii,pop}, allSNPs, 'stable');
        end
    end
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
    componentVariance = (1-1e-6) * ones(noPops) + 1e-6 * eye(noPops);
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
    totalNoSNPs = max(vertcat(SNPs{block,:}));
    whichCpt = randsample(1:length(componentWeight),totalNoSNPs,true,componentWeight);
    beta = zeros(totalNoSNPs,noPops);
    for cpt = 1:noCpts
        beta(whichCpt == cpt,:) = mvnrnd(zeros(1,noPops), componentVariance(:,:,cpt), sum(whichCpt == cpt));
    end
    
    for pop = 1:noPops
        beta_perallele{block,pop} = beta(:,pop);
        
        % standardized effect sizes (sd of Y per sd of X)
        beta_persd{block,pop} = beta(SNPs{block,pop},pop) .* ...
            sqrt(2 * alleleFrequency{block,pop} .* (1 - alleleFrequency{block,pop}));
    end
end

% Normalize effect sizes so they add up to h2
if ~isempty(heritability)
    for pop = 1:noPops
        normalizer(pop) = sqrt(heritability(pop)/sum(cellfun(@(x)sum(x.^2),beta_persd(:,pop))));
        beta_perallele(:,pop) = cellfun(@(b){b*normalizer(pop)}, beta_perallele(:,pop));
        beta_persd(:,pop) = cellfun(@(b){b*normalizer(pop)}, beta_persd(:,pop));
        
    end
    % normalize variance component matrices as well
    for cpt = 1:noCpts
        componentVariance(:,:,cpt) = componentVariance(:,:,cpt) .* (sqrt(normalizer) .* sqrt(normalizer)');
    end
end

% Sample summary statistics
for block = 1:noBlocks
    for pop = 1:noPops
        % marginal effect size estimates (i.e., sample correlations)
        if ~isempty(precisionMatrices)
            alphaHat{block,pop} = precisionMatrices{block,pop} \ beta_persd{block,pop}...
                + chol(precisionMatrices{block,pop}) \ randn(noSNPs(block,pop),1) / sqrt(sampleSize(pop));
        else
            alphaHat{block,pop} = correlationMatrices{block,pop} * beta_persd{block,pop}...
                + chol(correlationMatrices{block,pop})' * randn(noSNPs(block,pop),1)  / sqrt(sampleSize(pop));
        end
    end
end

end