function [h2, sigmasq, betaExpectation, betaVariance] =...
    h2EM(alphaHat, P, nn, annot, reps, tau)
%h2EM computes MLE of h2 using expectation-maximization
% Input arguments:
% alphahat: GWAS sumstats (m x 1);
% P: LD precision matrix (m x m sparse);
% nn: GWAS sample size (scalar);
% sigmasqSupport: possible values for tau (k x 1);
% sigmasqPrior: prior probability of
% each tau value; reps: number of steps to take
% Ouput arguments
% betaExpectation: expected value of beta given the data

mm = cellfun(@length,alphaHat);
noAnnot = size(annot{1},2);
noBlocks = length(annot);
annot_cat = vertcat(annot{:});
if nargin < 6
    tau = zeros(noAnnot,1);
    for ii=1:noAnnot
        tau(ii) = nn;
    end
end
mm_cum = [0; cumsum(mm)];
betaVariance = zeros(mm_cum(end),1);
betaExpectation = betaVariance;

for rep = 1:reps
    for block = 1:noBlocks
        denominator = (annot{block} * tau) .* P{block} + nn * speye(mm(block));
        
        % cov_Qtau(beta|alphaHat,tauExpectation) == P/denominator
        betaVariance(mm_cum(block)+1:mm_cum(block+1)) = ...
            sum(P{block}.*sparseinv(denominator),1)';
        
        % E_Qtau(beta|alphaHat,tauExpectation) == (P/denominator) * alphaHat
        betaExpectation(mm_cum(block)+1:mm_cum(block+1)) = ...
            nn * P{block} * (denominator \ alphaHat{block});
    end
    % E_Qbeta(log P(tau|beta)) == const + log prior + 1/2*log tau
    %     - 1/2 * tau * E(beta.^2)
    for ii=1:noAnnot
        tau(ii) = 1 / ( mean(betaExpectation(annot_cat(:,ii)).^2) + ...
            mean(betaVariance(annot_cat(:,ii))) );
    end
    
end

sigmasq = 1./tau;
h2 = sum(annot_cat * sigmasq);

end

