function [sumstats, whichIndices, true_beta_perallele, true_beta_perSD, mergedAnnot] = ...
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
% savePath: path and file name where summary statistics should be saved.
% One file is saved per population; if there are multiple populations, this
% should be a cell array of length noPops. Sumstats will be saved as a plain
% text file with the following columns:
%   SNP: identifier of each SNP, if snpID is specified
%   Z_deriv_allele: Z scores
%   N: sample size
%   AF_deriv_allele: if alleleFrequency is specified, the AF of each SNP
%   index: row/column of the precision matrix corresponding to each SNP;
%       indexing starts at zero
%   block: which LD block the row/col belongs to; indexing starts at
%       zero
%   beta_perSD_true: true causal effect size, in per-standard-deviation
%   units, for each SNP
%
% snplists: SNP lists containing additional information to be printed with the
% summary statistics. Must be specified in order to use 'PRScs'
% file format option
%
% fileFormat: If saving summary statistics to a file, which file format to
% use. Current options are 'ldgm', 'PRScs', and 'LDSC'
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
% annotations: annotation matrix for each LD block as a number-of-LD
% blocks by 1 cell array
%
% linkFn: link function mapping annotation vector of a SNP to its relative
% per-SNP heritability. Should be a nonnegative, scalar-valued function:
% for example, @(annot)log(1 + exp(annot * tau)), where tau is a column
% vector of length equal to the number of annotations
%
% whichIndicesAnnot: which rows/columns of the precision matrices/correlation
% matrices have a corresponding annotations vector, as a
% number-of-LD-blocks by 1 cell array. This is useful when you are
% performing simulations with real functional annotations, and they are
% missing for some of the SNPs in the LDGM. SNPs not in whichIndices will
% not be assigned an effect, and will not have summary statistics.
% 
% annotationDependentPolygenicity: if false (default), annotation-dependent
% and AF-dependent differences in per-SNP heritability will be simulated by
% scaling the effect-size magnitude of causal SNPs. If true, differences
% will be simulated by modifying the proportion of causal SNPs. When this
% is set to true, sum(componentWeight) should be small enough (e.g., 0.01)
% that it is possible to produce the deesired enrichments without exceeding
% the desired polygenicity; otherwise, a warning will be printed.
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

% save path for summary statistics
addParameter(p, 'savePath', '', @(s)ischar(s) || iscell(s));

% SNP lists containing additional information to be printed with the
% summary statistics. Must be specified in order to use 'PRScs'
% file format option
addParameter(p, 'snplists', {}, @iscell);

% If saving summary statistics to a file, which file format to use. Current
% options are 'ldgm' and 'PRScs'
addParameter(p, 'fileFormat', 'ldgm', @ischar);

% allele frequency for each LD block and each population as a number-of-LD
% blocks by number-of-populations cell array
addParameter(p, 'alleleFrequency', {}, @iscell);

% precision matrix for each LD block and each population as a number-of-LD
% blocks by number-of-populations cell array
addParameter(p, 'precisionMatrices', {}, @iscell);

% correlation matrix for each LD block and each population as a number-of-LD
% blocks by number-of-populations cell array
addParameter(p, 'correlationMatrices', {}, @iscell);

% 'alpha' parameter of e.g. Schoech et al. 2019 Nat Comm. If there are
% multiple populations, it uses mean AF to determine per-allele effect size
% variance
addParameter(p, 'alphaParam', 0, @isscalar);

% annotation matrix for each LD block as a number-of-LD
% blocks by 1 cell array
addParameter(p, 'annotations', {}, @iscell);

% link function mapping annotation vector of a SNP to its relative per-SNP
% heritability
addParameter(p, 'linkFn', @(x)ones(size(x)), @(f)isa(f,'function_handle'));

% which rows/columns of the precision matrices/correlation
% matrices have a corresponding annotations vector, as a
% number-of-LD-blocks by 1 cell array. This is useful when you are
% performing simulations with real functional annotations, and they are
% missing for some of the SNPs in the LDGM. SNPs not in whichIndices will
% not be assigned an effect, and will not have summary statistics.
addParameter(p, 'whichIndicesAnnot', {}, @iscell);

% Whether differences in per-SNP heritability due to alpha model scaling
% and annotations should be modeled by modifying the fraction of causal
% SNPs, or by scaling the causal effect sizes
addParameter(p, 'annotationDependentPolygenicity', 0, @isscalar);

% total heritability for each population, either as a scalar, a vector, or
% a square matrix. If a matrix, it specifies both the heritability for each
% population (along its diagonal) and the genetic correlations (off the
% diagonal, i.e. with corrcof(heritability) == r_pop). If a vector, the
% cross-population genetic correlations will be (almost) 1. If a scalar,
% the cross-population genetic correlations will be (almost) 1, and the
% heritability will be the same for each population.
addParameter(p, 'heritability', 1, @isnumeric);

% per-allele effect-size covariance matrix for each mixture component,
% either as a vector or as a number-of-populations by number-of-populations by
% number-of-components array. Will be scaled to match total heritability in
% each population.
addParameter(p, 'componentVariance', 1, ...
    @(x)isnumeric(x) & all(x>=0,'all') & (isvector(x) | ndims(x) == 3));

% mixture weight for each heritability component, as a number-of-components
% by 1 vector. If sum is smaller than one, a null
% component is added with componentVariance equal to zero.
addParameter(p, 'componentWeight', 1, ...
    @(x)isvector(x) & sum(x,'all')<=1 & all(x>=0,'all'));

% fraction of missing SNPs (missing at random)
addParameter(p, 'missingness', 0, @isscalar);

% whether to compute sample allele frequencies with added noise
addParameter(p, 'noisySampleAF', false, @isscalar);


% turns p.Results.x into just x
parse(p, sampleSize, varargin{:});
fields = fieldnames(p.Results);
for k=1:numel(fields)
    line = sprintf('%s = p.Results.%s;', fields{k}, fields{k});
    eval(line);
end

% input handling
assert(nargin > 0 | ~isempty(savePath), 'No output requested')

if isempty(precisionMatrices) && isempty(correlationMatrices) %#ok<*USENS>
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
if isempty(annotations) %#ok<*NODEF>
    annotations = arrayfun(@(n){ones(n,1)},noSNPs);
end
noAnnot = size(annotations{1},2);

% Construct whichIndicesAnnot if needed
if isempty(whichIndicesAnnot)
    whichIndicesAnnot = arrayfun(@(m)(1:m)',noSNPs,'UniformOutput',false);
end

% Convert whichIndicesAnnot from logicals to indices if needed
if islogical(whichIndicesAnnot{1})
    whichIndicesAnnot = cellfun(@find,whichIndicesAnnot,'UniformOutput',false);
end

% Number of SNPs in each block with an annotation vector
noSNPsAnnot = cellfun(@length,whichIndicesAnnot);

% Verify that annotation matrix looks OK
for block = 1:noBlocks
    assert(all(size(annotations{block}) == [noSNPsAnnot(block), noAnnot]))
end

if ~isempty(snplists)
    assert(numel(snplists) == noBlocks,...
        'Specify snplists as a number-of-blocks by 1 cell array')
    for block = 1:noBlocks
        [ui, representatives] = unique(snplists{block}.index);
        assert(all(ui==(0:noSNPs(block)-1)'),'Something wrong with SNP list indices')
        snplists{block} = snplists{block}(representatives,:);
    end
end

% Heritability for each population
if isscalar(heritability) && noPops > 1
    heritability = heritability * ones(1,noPops);
end

% Cross-population heritability matrix
if isvector(heritability) && noPops > 1
    h = sqrt(heritability);
    M = (1-1e-6) * ones(noPops) + 1e-6 * eye(noPops);
    heritability = h' .* M .* h;
end

% componentWeight should specify a convex combination of components
assert(all(componentWeight<=1) & all(componentWeight>=0))
assert(sum(componentWeight)<=1)
noCpts = length(componentWeight);

% If componentVariance is specified as a vector, convert it into a
% noPops x noPops x noCpts array with cross-population correlations nearly
% 1
if isvector(componentVariance)
    componentVariance = reshape(componentVariance,1,1,length(componentVariance)) .* ...
        repmat(ones(noPops)+1e-6*eye(noPops),1,1,noCpts);
else
    assert(ndims(componentVariance) == 3,...
        'Specify componentVariance as either a vector or as an array of size noPops x noPops x noCpts')

    assert(all(size(componentVariance) == [noPops, noPops, noCpts]),...
        'Specify componentVariance as either a vector or as an array of size noPops x noPops x noCpts')
end

% Get rid of null components
componentVarianceNotNull = reshape(any(componentVariance,1:2),noCpts,1);
componentVariance = componentVariance(:,:,componentVarianceNotNull);
componentWeight = componentWeight(componentVarianceNotNull);
noCpts = length(componentWeight);

% mixture component assignments, which may be annotation- and
% AF-dependent
totalProbabilityCausal = sum(componentWeight(componentVarianceNotNull));
assert(totalProbabilityCausal <= 1, 'Component weights must sum to at most 1')
if ~annotationDependentPolygenicity
    probabilityCausal = arrayfun(@(n){totalProbabilityCausal * ones(n,1)},noSNPsAnnot);
else
    % annotation-dependent probability of being assigned to a non-null
    % component
    relativeProbabilityCausal = cellfun(linkFn,annotations,'UniformOutput',false); 
    
    % allele-frequency-dependent scaling
    if alphaParam ~= -1
        assert(~isempty(alleleFrequency), ...
            'If alphaParam is specified, allele frequencies must be specified')
        for block=1:noBlocks
            % mean allele frequency across populations
            meanAF = max(1e-9, min(1-1e-9, mean([alleleFrequency{block,:}],2)));
            meanAF = meanAF(whichIndicesAnnot{block});

            % scale relativeProbabilityCausal by alpha model factor
            relativeProbabilityCausal{block} = relativeProbabilityCausal{block} .* ...
                (meanAF.*(1-meanAF)).^(1 + alphaParam);
        end
    end
    
    % scale relativeProbabilityCausal to sum to totalProbabilityCausal
    scalingFactor = totalProbabilityCausal * sum(noSNPsAnnot) / sum(cellfun(@sum,relativeProbabilityCausal));
    probabilityCausal = cellfun(@(p){p * scalingFactor}, relativeProbabilityCausal);
    
    % probabilityCausal might have entries greater than 1, which are
    % truncated with a warning
    fractionSNPsProbabilityCausalGTOne = mean(cellfun(@(p)mean(p>1),probabilityCausal));
    if any(cellfun(@max,probabilityCausal) > 1)
        warning(['Around %.3f of SNPs had causal probabilities greater than 1, ',...
            'and are being truncated. Address this by reducing polygenicity ',...
            'or changing link function'], fractionSNPsProbabilityCausalGTOne)
        probabilityCausal = cellfun(@(x)min(x,1), probabilityCausal, 'UniformOutput', false);
    end
end

% Simulate summary statistics for each LD block
true_beta_perallele = arrayfun(@(m)zeros(m,1),repmat(noSNPs,1,noPops),'UniformOutput',false);
true_beta_perSD = true_beta_perallele;
for block = 1:noBlocks

    % mixture component assignments
    whichSNPsCausal = rand(noSNPsAnnot(block),1) < probabilityCausal{block};
    whichCpt = zeros(noSNPsAnnot(block),1);
    whichCpt(whichSNPsCausal) = randsample(1:noCpts,sum(whichSNPsCausal),true,componentWeight);

    % sample beta from respective mixture components
    beta = zeros(noSNPsAnnot(block),noPops);
    for cpt = 1:noCpts
        beta(whichCpt == cpt,:) = mvnrnd(zeros(1,noPops), componentVariance(:,:,cpt), sum(whichCpt == cpt));
    end
    
    % mean allele frequency across populations
    if ~isempty(alleleFrequency)
        meanAF = max(1e-9, min(1-1e-9, mean([alleleFrequency{block,:}],2)));
        meanAF = meanAF(whichIndicesAnnot{block});
    end

    % scale causal effects to desired variance
    if annotationDependentPolygenicity && ~isempty(alleleFrequency)
        % larger per-allele effect-size magnitudes for low-frequency causal SNPs
        beta = beta ./ sqrt((meanAF.*(1-meanAF)));
    else
        % apply alpha model scaling
        if ~isempty(alleleFrequency)
            beta = beta .* sqrt((meanAF.*(1-meanAF)).^alphaParam);
        end
        % scale effect sizes using linkFn for SNPs in annotation matrix
        beta = beta .* sqrt(linkFn(annotations{block}));
    end

    assert(isreal(beta) & all(beta == beta, 'all'), 'Imaginary or NaN beta; check link function')

    % Assign normalized and per-allele betas for each population
    for pop = 1:noPops
        true_beta_perallele{block,pop}(whichIndicesAnnot{block}) = beta(:,pop);

        if ~isempty(alleleFrequency)
            true_beta_perSD{block,pop}(whichIndicesAnnot{block}) = beta(:,pop) .* ...
                sqrt(2 * alleleFrequency{block,pop}(whichIndicesAnnot{block}) .* ...
                (1 - alleleFrequency{block,pop}(whichIndicesAnnot{block})));
        
        elseif ~isempty(precisionMatrices)% zero out effects for SNPs not in precision matrix
            nzAF = any(precisionMatrices{block,pop},2); 
            true_beta_perSD{block,pop}(whichIndicesAnnot{block}) = ...
                beta(:,pop) .* nzAF(whichIndicesAnnot{block});
        else
            true_beta_perSD{block,pop}(whichIndicesAnnot{block}) = beta(:,pop);
        end
    end

end

% Normalize effect sizes so they add up to h2
if ~isempty(heritability)
    for pop = 1:noPops
        normalizer = sqrt(heritability(pop,pop)/sum(cellfun(@(x)sum(x.^2),true_beta_perSD(:,pop))));
        true_beta_perallele(:,pop) = cellfun(@(b){b*normalizer}, true_beta_perallele(:,pop));
        true_beta_perSD(:,pop) = cellfun(@(b){b*normalizer}, true_beta_perSD(:,pop));
    end
end

% Sample summary statistics
Z = cell(noBlocks,noPops);
whichIndices = Z;
mergedAnnot = cell(size(annotations));
for block = 1:noBlocks
    for pop = 1:noPops

        if ~isempty(precisionMatrices)
            % SNPs not missing in precision matrix
            pnz = diag(precisionMatrices{block,pop}) ~= 0;
            incl = whichIndicesAnnot{block}(pnz(whichIndicesAnnot{block}));
            whichIndices{block,pop} = incl;

            % Simulate sumstats using precision matrices
            Z{block,pop} = precisionDivide(precisionMatrices{block,pop},...
                true_beta_perSD{block,pop}(incl) * sqrt(sampleSize(pop)), incl);
            noise = chol(precisionMatrices{block,pop}(pnz,pnz)) \...
                randn(sum(pnz),1);
            idx = lift(whichIndices{block,pop},find(pnz));
            Z{block,pop} = Z{block,pop} + noise(idx);

        else
            % SNPs not missing in correlation matrix
            pnz = diag(correlationMatrices{block,pop}) ~= 0;
            incl = whichIndicesAnnot{block}(pnz(whichIndicesAnnot{block}));
            whichIndices{block,pop} = incl;

            % Simulate using correlation matrices
            Z{block,pop} = correlationMatrices{block,pop}(incl,incl) *...
                true_beta_perSD{block,pop}(incl) * sqrt(sampleSize(pop))...
                + chol(correlationMatrices{block,pop}(incl,incl))' *...
                randn(sum(incl),1);
        end

        % Additional SNPs missing at random
        if missingness > 0
            incl = rand(length(incl),1) > missingness;
            Z{block,pop} = Z{block,pop}(incl);
            whichIndices{block,pop} = whichIndices{block,pop}(incl);
        end

        if nargout >= 5
            idx = lift(whichIndices{block,pop},whichIndicesAnnot{block});
            mergedAnnot{block} = annotations{block}(idx,:);
        end
    end
end

noNonmissingSNPs = cellfun(@length,whichIndices);

% Store output in a table
sumstats = cell(noBlocks,noPops);
for block = 1:noBlocks
    for pop = 1:noPops
        sumstats{block,pop} = table('size',[length(Z{block,pop}),0]);
        if ~isempty(snplists)
            sumstats{block,pop}.SNP = snplists{block}.site_ids(whichIndices{block,pop});
        end
        sumstats{block,pop}.Z_deriv_allele = Z{block,pop};
        sumstats{block,pop}.N = sampleSize(pop) * ones(size(Z{block,pop}));
        if ~isempty(alleleFrequency)
            if noisySampleAF
                sumstats{block,pop}.AF_deriv_allele = ...
                    binornd( 2*sumstats{block,pop}.N, ...
                    alleleFrequency{block,pop}(whichIndices{block,pop}) ) ./ ...
                    (2*sumstats{block,pop}.N);
            else
                sumstats{block,pop}.AF_deriv_allele = ...
                    alleleFrequency{block,pop}(whichIndices{block,pop});
            end
        end
    end
end

% Save to file if requested
if ~isempty(savePath)
    if noPops > 1
        assert(iscell(savePath) & numel(savePath) == noPops)
    elseif ischar(savePath)
        savePath = {savePath};
    end
    true_beta_perSD_nonmissing = cellfun(@(x,j)x(j),true_beta_perSD,whichIndices,'UniformOutput',false);
    for pop = 1:noPops
        % ldgm sumstats file format output
        if strcmpi(fileFormat,'ldgm')
            T = vertcat(sumstats{:,pop});
            T.index(:) = vertcat(whichIndices{:,pop}) - 1; % zero-indexed
            whichBlock = arrayfun(@(n,s)n*ones(s,1),(1:noBlocks)',...
                noNonmissingSNPs,'UniformOutput',false);
            T.block(:) = vertcat(whichBlock{:}) - 1; % zero-indexed
            T.beta_perSD_true(:) = vertcat(true_beta_perSD_nonmissing{:,pop});

        % PRS-CS file format ouput
        elseif strcmpi(fileFormat,'PRScs')
            assert(~isempty(snplists),'To use PRScs file format, SNP lists must be specified')
            snplistsCat = cellfun(@(T,j)T(j,:),snplists(:,pop),whichIndices(:,pop),'UniformOutput',false);
            snplistsCat = vertcat(snplistsCat{:,pop});
            nn = height(snplistsCat);
            T = table('size',[nn,0]);
            T.SNP = snplistsCat.site_ids;
            T.A1 = snplistsCat.anc_alleles;
            T.A2 = snplistsCat.deriv_alleles;
            sumstatsCat = vertcat(sumstats{:,pop});
            T.BETA = sumstatsCat.Z_deriv_allele;
            T.P = chi2cdf(sumstatsCat.Z_deriv_allele.^2,1,'upper');
            if any(sumstatsCat.Z_deriv_allele.^2 > 300)
                warning('Very large chi^2 statistics might lead to unreliable p-values in printed summary statistics')
            end
            % Discard NA SNP IDs
            T = T(~strcmpi(T.SNP,'NA'),:);

        % LDSC file format output
        elseif strcmpi(fileFormat,'LDSC')
            assert(~isempty(snplists),'To use LDSC file format, SNP lists must be specified')
            snplistsCat = cellfun(@(T,j)T(j,:),snplists(:,pop),whichIndices(:,pop),'UniformOutput',false);
            snplistsCat = vertcat(snplistsCat{:,pop});
            nn = height(snplistsCat);
            T = table('size',[nn,0]);
            T.SNP = snplistsCat.site_ids;
            T.A1 = snplistsCat.anc_alleles;
            T.A2 = snplistsCat.deriv_alleles;
            sumstatsCat = vertcat(sumstats{:,pop});
            T.Z = sumstatsCat.Z_deriv_allele;
            T.N(:) = sampleSize;

            % Discard NA SNP IDs
            T = T(~strcmpi(T.SNP,'NA'),:);
        else
            error('Current file format options are PRScs, ldgm, and LDSC')
        end

        writetable(T,savePath{pop},'FileType','text','delimiter','\t');
    end
end

end