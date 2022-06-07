function [betaExpectation,betaVariance,tauPosterior,tauExpectation] =...
    variationalGWASx(alphaHat, P, nn, sigmasqSupport, sigmasqPrior, reps)
%variationalGWAS computes a variational approximation to the posterior
%distribution of the effect sizes, beta, and their variance parameters,
%sigamsq: P(beta,sigmasq|data) ~= Qbeta(beta) Qsigmasq(sigmasq)
% Input arguments:
% alphahat: GWAS sumstats (m x 1); 
% P: LD precision matrix (m x m sparse); 
% nn: GWAS sample size (scalar); 
% sigmasqSupport: possible values for tau (k x 1); 
% sigmasqPrior: prior probability of
% each tau value; reps: number of steps to take
% Ouput arguments
% betaExpectation: expected value of beta given the data

mm = length(alphaHat);
if iscolumn(sigmasqPrior)
    sigmasqPrior = sigmasqPrior';
end
if iscolumn(sigmasqSupport)
    sigmasqSupport = sigmasqSupport';
end
tauSupport = cellfun(@inv,sigmasqSupport);
tauExpectation = ones(1,mm)/sum(sigmasqPrior .* sigmasqSupport);

for rep = 1:reps
    tauExpectation = speye(mm) .* tauExpectation;
    denominator = tauExpectation * P + nn * speye(mm);
    
    % cov_Qtau(beta|alphaHat,tauExpectation) == P/denominator
    betaVariance = sum(P.*sparseinv(denominator),1)';
    
    % E_Qtau(beta|alphaHat,tauExpectation) == (P/denominator) * alphaHat
    betaExpectation = nn * P * (denominator \ alphaHat);
    
    % E_Qbeta(log P(tau|beta)) == const + log prior + 1/2*log tau 
    %     - 1/2 * tau * E(beta.^2)
    tauPosterior = sigmasqPrior .* sqrt(tauSupport) .* ...
        exp(-1/2 * tauSupport .* (betaVariance + betaExpectation.^2));
    tauPosterior = tauPosterior ./ sum(tauPosterior,2);
    % E_Qbeta(tau)
    tauExpectation = sum(tauPosterior .* tauSupport,2);
    
end

% [~, bestTau] = max(tauPosterior);
% tauEM = tauSupport(bestTau);
% betaVariance = inv(diag(tauEM) + nn*inv(omega));
% betaEM = betaVariance * nn * alphaHat;
% betaExpectation = omega * alphaVariance * nn * omega * alphaHat;
end

