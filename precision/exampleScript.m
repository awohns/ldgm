% example script to infer LD precision matrix from simulated data
numNodes = 1e3; 
density = 2e-2;
numHaplotypes = 2*numNodes; 

% For compatibility with plotting script
allele_freq = ones(numNodes,1);

% True precision matrix
omega = makeSparsePrecision(numNodes,density,numNodes*10);
A = omega~=0;

% Draw MVN samples with precision matrix omega
X = randn( numHaplotypes, numNodes) * chol(inv(omega));
X = (X - mean(X,1));
R = corr(X);

% Estimate precision matrix
tol = 1e-4;
[omegaEst, pval] = LDPrecision(R,A,numHaplotypes,tol);

mse = mean((omega(A) - omegaEst(A)).^2) / mean(omega(A).^2);
fprintf('P-value for model fit: %f\n Precision matrix percent mean squared error: %f\n',pval,mse);

precisionplotscript
