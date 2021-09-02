addpath(genpath('.'))
path_prefix='~/Dropbox/Pouria/data/';

% Genotype matrix
X=load([path_prefix, 'genomat'])';

% Adjacency matrix
A_weighted=importGraph([path_prefix,'adjlist'],1);

% Empty rows/columns correspond to duplicate SNPs (on same brick as
% another SNP)
SNPs = find(any(A_weighted));
X = X(:,SNPs);
A_weighted = A_weighted(SNPs,SNPs);
[numHaplotypes, numSNPs] = size(X);

% Edges to retain
desired_density = 0.2;

% Threshold weighted network to desired density, and add self-edges
max_density = nnz(A_weighted) / length(A_weighted)^2;
threshold = quantile(nonzeros(A_weighted),...
    max(0, 1 - desired_density / max_density));
A = A_weighted + speye(numSNPs) > threshold;


% set missing values to the mean genotype value
missing = X==-1; 
if any(missing(:))
    warning('Some genotypes missing')
    X(missing) = 0;
    allele_freq = repmat(sum(X)./sum(~missing),numHaplotypes,1);
    X(missing) = allele_freq(missing);
    allele_freq = allele_freq(1,:);
else
    allele_freq = mean(X);
end


if numSNPs > 1e3
    % LD matrix for edges of A
    [ii,jj] = find(A);
    X = (X - mean(X,1))./std(X);
    R = arrayfun(@(i,j)dot(X(:,i),X(:,j)),ii,jj)/numHaplotypes;
    R = sparse(ii,jj,R);
else
    R = corr(X);
end

% estimated precision matrix
maxReps = 1e3;
omegaEst = speye(size(A));
tic;[omegaEst, pval] = LDPrecision(R,A,numHaplotypes,maxReps, omegaEst, 0);toc
tic;[omegaEst] = LDPrecisionDTraceLoss(R + .1 * speye(numSNPs),A,numHaplotypes,maxReps, omegaEst, 0);toc

% Plot stuff for SNPs with frequency above threshold
AF_threshold = .05;
precisionplotscript
