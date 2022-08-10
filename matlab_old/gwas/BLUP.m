function [betaExpectation, h2MLE] =...
    BLUP(alphahat, P, nn, tau, whichSNPs)
% BLUP computes the best linear unbiased predictor, E(beta|GWAS, gaussian
% prior) where the GWAS data is alphaHat, the LD precision matrix is P, and
% beta~N(0,tau^-1).
% betaExpectation: expected value of beta given the data
% alphahat: GWAS sumstats; omega: LD precision matrix; nn: GWAS sample
% size; tau: precision of beta; whichSNPs: rows/cols of P for which we have
% data alphaHat
% h2MLE: vector of MLE values for diag(beta*alpha'). Can be negative. Add
% these up to get an MLE heritability estimate for a set of SNPs.

mm = length(P);
if isvector(tau)
    sigmasq = zeros(mm,1);
    sigmasq(whichSNPs) = 1./tau;
    sigmasq = diag(sparse(sigmasq));
end
betahat = zeros(mm,1);
betahat(whichSNPs) = P(whichSNPs,whichSNPs) * alphahat - ...
    P(whichSNPs,~whichSNPs) * (P(~whichSNPs,~whichSNPs) \ ...
    (P(~whichSNPs,whichSNPs)*alphahat));

% betaExpectation = (tauMatrix + nn*inv(P)) \ (nn * alphaHat);
betaExpectation = sigmasq * ((P + nn*sigmasq) \ (nn * betahat));

if nargout > 1
    betaExpectation(~whichSNPs) = 0;
    alphaExpectation = P \ betaExpectation;
    h2MLE = betaExpectation(whichSNPs) .* alphaExpectation(whichSNPs);
end

betaExpectation = betaExpectation(whichSNPs);

end

