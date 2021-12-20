function [betaExpectation] =...
    BLUP(alphahat, P, nn, tau, whichSNPs)
% BLUP computes the best linear unbiased predictor, E(beta|GWAS, gaussian
% prior) where the GWAS data is alphaHat, the LD precision matrix is P, and
% beta~N(0,tau^-1).
% betaExpectation: expected value of beta given the data
% alphahat: GWAS sumstats; omega: LD precision matrix; nn: GWAS sample
% size; tau: precision of beta; whichSNPs: rows/cols of P for which we have
% data alphaHat

mm = length(P);
if isvector(tau)
    tauMatrix = zeros(mm,1);
    tauMatrix(whichSNPs) = tau;
    tauMatrix = diag(sparse(tau));
end
betahat = zeros(mm,1);
betahat(whichSNPs) = P(whichSNPs,whichSNPs) * alphahat - ...
    P(whichSNPs,~whichSNPs) * (P(~whichSNPs,~whichSNPs) \ ...
    (P(~whichSNPs,whichSNPs)*alphahat));

% betaExpectation = (tau + nn*inv(omega)) \ (nn * alphaHat);
betaExpectation = (P*tauMatrix + nn*speye(mm)) \ (nn * betahat);
betaExpectation = betaExpectation(whichSNPs);

end

