addpath(genpath('../precision'))
path_prefix='~/Dropbox/Pouria/data/';

% Genotype matrix
X=load([path_prefix, 'genomat'])';

% Weighted adjacency matrix
import_weighted = 1;
A_weighted=importGraph([path_prefix,'adjlist'], import_weighted);

% Empty rows/columns correspond to duplicate SNPs (on same brick as
% another SNP)
SNPs = find(any(A_weighted));
X = X(:,SNPs);
if any(X(:)==-1)
    error('Missing genotypes not supported')
end
A_weighted = A_weighted(SNPs,SNPs);
[numHaplotypes, numSNPs] = size(X);
allele_freq = mean(X);

% How many edges to retain
desired_density = 0.2;

% Threshold weighted network to desired density, and add self-edges
max_density = nnz(A_weighted) / length(A_weighted)^2;
threshold = quantile(nonzeros(A_weighted),...
    max(0, 1 - desired_density / max_density));
A = A_weighted + speye(numSNPs) > threshold;

% LD matrix for edges of A
[ii,jj] = find(A);
X = (X - mean(X));
X = X./sqrt(mean(X.^2));
R = arrayfun(@(i,j)dot(X(:,i),X(:,j)),ii,jj)/numHaplotypes;
R = sparse(ii,jj,R);


% estimated precision matrix
maxReps = 1e3;
convergenceTol=1e-6;
precisionEstimate = speye(size(A));

[precisionEstimate, obj_val] = LDPrecision(R,'P0',precisionEstimate,'max_steps',maxReps,...
    'convergence_tol',convergenceTol,'printstuff',true);

% Plot stuff for SNPs with frequency above threshold
AF_threshold = 0.05;
plotting_script
