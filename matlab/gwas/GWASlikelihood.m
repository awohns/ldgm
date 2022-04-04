function [logLikelihood, logLikelihoodBlocks] = GWASlikelihood(alphahat,sigmasq,P,nn,whichSNPs)
% GWASlikelihood computes likelihood of the GWAS sumstats alphahat under a
% gaussian model:
%                   beta ~ MVN(mu,diag(sigmasq))
%                   alphahat|beta ~ MVN(R*beta, R/nn)
%                   inv(R) = P.
% The GWAS SNPs in alphahat should be a subset of those in P, and in the
% same order; boolean vector whichSNPs should be true for rows/columns of P
% corresponding to one of these SNPs, false elsewhere.
% Currently, sigmasq should be the same size as alphahat, such that missing
% SNPs are modeled as having zero effect size.
%
% Optionally, arguments alphahat, sigmasq, P, whichSNPs can be specified as
% cell arrays with one cell for each LD block.
%
% Units of alphahat should be either (1) normalized effect size estimates,
% ie sample correlations, or (2) Z scores. In case (1), nn should be the
% number of individuals in the GWAS. In case (2), nn should be 1, and the
% vector of sigmasq values should be multiplied by nn.
%
% If it is desired to specify the mean of beta, ie
%                   beta ~ MVN(mu, diag(sigmasq))
% call GWASlikelihood( alphahat - P \ mu, sigmasq, P, ...)

assert(isscalar(nn) && all(nn>0),'Sample size nn should be a positive scalar')

if iscell(P)
    assert(iscell(alphahat) && iscell(whichSNPs) && ...
        (iscell(sigmasq) || isscalar(sigmasq)), ...
        'If P is a cell array then alphahat, sigmasq, whichSNPs should also cell arrays of the same size');
    if isscalar(sigmasq) && ~iscell(sigmasq)
        sigmasq = cellfun(@(b){ones(size(b))*sigmasq},alphahat);
    end
    assert(all(size(P) == size(sigmasq)) && all(size(P) == size(alphahat))...
        && all(size(P) == size(whichSNPs)), 'Input cell arrays should have same size');
    
    % Iterate over cell arrays
    logLikelihoodBlocks = cellfun(@(a,s,p,w)likelihoodFn(a,s,p,nn,w), ...
        alphahat, sigmasq, P, whichSNPs);
    
    logLikelihood = sum(logLikelihoodBlocks);
else
    logLikelihood = likelihoodFn(alphahat,sigmasq,P,nn,whichSNPs);
    logLikelihoodBlocks = logLikelihood;
end

    function ll = likelihoodFn(alphahat,sigmasq,P,nn,whichSNPs)
        
        mm = length(whichSNPs);
        mmz = sum(whichSNPs);
        assert(all(sigmasq>=0),'sigmasq should be nonnegative')
        assert(mmz == length(alphahat) && mmz == length(sigmasq),...
            'whichSNPs should be a boolean vector with sum equal to length of alphahat and sigmasq')
        assert(all(length(whichSNPs) == size(P)))
        
        
        % inv(P)(whichSNPs,whichSNPs) * alphahat
        betahat = P(whichSNPs,whichSNPs) * alphahat - P(whichSNPs,~whichSNPs) * ...
            (P(~whichSNPs,~whichSNPs) \ (P(~whichSNPs,whichSNPs) * alphahat));
        
        % diag([sigmasq,0,...]) + P/nn)
        DplusP = sparse(zeros(mm,1));
        DplusP(whichSNPs) = sigmasq;
        DplusP = diag(DplusP) + P/nn;
        A = chol(DplusP);
        
        % log|1/nn*P/P11 + diag(sigmasq)| == log|DplusP| - log|P11|
        logdetDplusP = 2*sum(log(diag(A)));
        logdetP11 = 2*sum(log(diag(chol(1/nn*P(~whichSNPs,~whichSNPs)))));
        
        % betahat' * PplusD\betahat == x'*x
        y = zeros(mm,1);
        y(whichSNPs) = betahat;
        x = A' \ y;
        x = x(whichSNPs);
        
        ll = 1/2 * (-(logdetDplusP - logdetP11) - x'*x - mmz*log(2*pi));
    end


end

