% rng('default')
addpath(genpath('.'))
path_prefix='~/Dropbox/Pouria/data/';

% Adjacency matrix
A=importGraph([path_prefix,'mutgraph_n238']);

% Genotype matrix
X=load([path_prefix, 'genomat_n238'])';

% indices with zero diagonal correspond to SNPs that were on the same brick
% as another SNP
SNPs = find(any(A));
X = X(:,SNPs);
A = A(SNPs,SNPs);

[numHaplotypes, numNodes] = size(X);

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


if numNodes > 1e3
    % LD matrix for edges of A
    [ii,jj] = find(A);
    X = (X - mean(X,1))./std(X);
    R = arrayfun(@(i,j)dot(X(:,i),X(:,j)),ii,jj)/numHaplotypes;
    R = sparse(ii,jj,R);
else
    R = corr(X);
end

% estimated precision matrix
tol = 1e-4;
% [omegaEst, pval] = LDPrecision(R, omega~=0, numHaplotypes, reps, speye(size(omega)));
tic;[omegaEst, pval] = LDPrecision(R,A,numHaplotypes,tol, speye(numNodes));toc

% Plot stuff for SNPs with frequency above threshold
AF_threshold = 0.01;
precisionplotscript
