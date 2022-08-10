function [betaExpectation] =...
    BLUPslow(alphahat, R, nn, tau, whichSNPs)
% BLUP computes the best linear unbiased predictor, E(beta|GWAS, gaussian
% prior) where the GWAS data is alphaHat, the LD precision matrix is P, and
% beta~N(0,tau^-1).
% betaExpectation: expected value of beta given the data
% alphahat: GWAS sumstats; omega: LD precision matrix; nn: GWAS sample
% size; tau: precision of beta; whichSNPs: rows/cols of P for which we have
% data alphaHat
% h2MLE: vector of MLE values for diag(beta*alpha'). Can be negative. Add
% these up to get an MLE heritability estimate for a set of SNPs.

if nargin < 5
    whichSNPs = true(size(alphahat));
end
mm = length(alphahat);

betaExpectation = nn * ((nn*R(whichSNPs,whichSNPs) + tau.*speye(mm)) \ alphahat);


end

