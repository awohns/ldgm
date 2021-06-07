% rng('default')

numNodes = 1e2; % very slow when this is more than a few hundred
density = 1e-1;

% seems this needs to be >= numNodes^2 in order for p-values to be 
% well-calibrated, but out-of-sample log-likelihoods should still be OK
numHaplotypes = 1e3; 

% threshold to throw out some edges as positive control
wrongGThreshold = .1;

omega = makeSparsePrecision(numNodes,density,numNodes*10);
% disp(omega)

% Draw samples with precision matrix omega
A=chol(inv(omega));
X=A'*randn(numNodes, numHaplotypes);
X = X - mean(X,2);

% estimated precision matrix
gradsteps=1e1;
[omegaEst, pval] = LDPrecision(X', omega~=0, gradsteps);

% pval for wrong graphical model
[omegaWrongG, pvalWrongG] = LDPrecision(X', abs(omega)>wrongGThreshold, gradsteps);

% pval for true LD precision matrix
ObjFn = @(Om) 0.5 * numHaplotypes * log(det(Om)) - 0.5 * sum(sum((Om*X) .* X));
omegaFull = inv(X*X'/numHaplotypes);
pvalTrueOmega = chi2cdf(2*(ObjFn(omegaFull)-ObjFn(omega)), (numNodes*(numNodes+1)/2), 'upper');

% Expected is null, non-null, null
disp('P-values for G, wrong G, true omega:')
disp([(pval) (pvalWrongG) (pvalTrueOmega)])

% out-of-sample evaluation
numOutOfSampleHaplotypes = 1000;
Xtesting=A'*randn(numNodes, numOutOfSampleHaplotypes);
ObjFnIndv = @(Om) 0.5 * log(det(Om)) - 0.5 * sum((Om*Xtesting) .* Xtesting)';

logLikelihood = ObjFnIndv(omegaEst);
logLikelihoodWrongG = ObjFnIndv(omegaWrongG);
logLikelihoodFullG = ObjFnIndv(omegaFull);
logLikelihoodTrueG = ObjFnIndv(omega);

% Expected is positive, positive, negative
disp('Out-of-sample relative log likelihood ratio for wrong G, full graphical model, true omega:')
disp(mean(logLikelihood - [logLikelihoodWrongG logLikelihoodFullG logLikelihoodTrueG]))
disp('Standard errors:')
disp(std(logLikelihood - [logLikelihoodWrongG logLikelihoodFullG logLikelihoodTrueG])/sqrt(numOutOfSampleHaplotypes))

