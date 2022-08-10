function [logLikelihood, logLikelihoodBlocks] = GWASlikelihood(Z_deriv_allele,sigmasq,P,nn,whichIndices)
% GWASlikelihood computes likelihood of the GWAS sumstats, Z_deriv_allele,
% under a gaussian model:
%                   beta ~ MVN(mu,diag(sigmasq))
%                   Z|beta ~ MVN(sqrt(nn)*R*beta, R)
%                   inv(R) = P.
% The GWAS SNPs in Z_deriv_allele should be a subset of those in P, and in the
% same order; whichIndices should be true for rows/columns of P
% corresponding to one of these SNPs, false elsewhere (or specify indices
% instead of logicals)
% 
% sigmasq should be the same size as Z_deriv_allele, and missing SNPs are 
% modeled as having zero effect size.
%
% Optionally, arguments Z_deriv_allele, sigmasq, P, whichSNPs can be specified as
% cell arrays with one cell for each LD block.
%
% If it is desired to specify the mean of beta, ie
%                   beta ~ MVN(mu, diag(sigmasq))
% call GWASlikelihood( Z - P \ mu * sqrt(nn), sigmasq, P, ...)

assert(isscalar(nn) && all(nn>0),'Sample size nn should be a positive scalar')

if iscell(P)
    assert(iscell(Z_deriv_allele) && iscell(whichIndices) && ...
        (iscell(sigmasq) || isscalar(sigmasq)), ...
        'If P is a cell array then alphahat, sigmasq, whichSNPs should also cell arrays of the same size');
    if isscalar(sigmasq) && ~iscell(sigmasq)
        sigmasq = cellfun(@(b){ones(size(b))*sigmasq},Z_deriv_allele);
    end
    assert(all(size(P) == size(sigmasq)) && all(size(P) == size(Z_deriv_allele))...
        && all(size(P) == size(whichIndices)), 'Input cell arrays should have same size');
    
    % Iterate over cell arrays
    logLikelihoodBlocks = cellfun(@(a,s,p,w)likelihoodFn(a,s,p,nn,w), ...
        Z_deriv_allele, sigmasq, P, whichIndices);
    
    logLikelihood = sum(logLikelihoodBlocks);
else
    logLikelihood = likelihoodFn(Z_deriv_allele,sigmasq,P,nn,whichIndices);
    logLikelihoodBlocks = logLikelihood;
end

    function ll = likelihoodFn(Z,sigmasq,P,nn,whichSNPs)
        
        if ~islogical(whichSNPs)
            whichSNPs = unfind(whichSNPs,length(P));
        end
        mm = length(whichSNPs);
        mmz = sum(whichSNPs);
        assert(all(sigmasq>=0),'sigmasq should be nonnegative')
        assert(mmz == length(Z) && mmz == length(sigmasq),...
            'whichSNPs should be a boolean vector with sum equal to length of alphahat and sigmasq')
        assert(all(length(whichSNPs) == size(P)))
        
        
        % inv(P)(whichSNPs,whichSNPs) * Z
        x = precisionMultiply(P,Z,whichSNPs);
        
        % diag([sigmasq,0,0,...]) + P)
        DplusP = sparse(zeros(mm,1));
        DplusP(whichSNPs) = sigmasq * nn;
        DplusP = diag(DplusP) + P;
        
        % handling SNPs missing from P
        incl = diag(P)~=0;
        otherSNPs = incl & (~whichSNPs);
        
        A = chol(DplusP(incl,incl));
        
        % log|1/nn*P/P11 + diag(sigmasq)| == log|DplusP| - log|P11|
        logdetDplusP = 2*sum(log(diag(A)));
        logdetP11 = 2*sum(log(diag(chol(P(otherSNPs,otherSNPs)))));
        
        % betahat' * PplusD\betahat == x'*x
        y = zeros(mm,1);
        y(whichSNPs) = x;
        x = A' \ y(incl);
        x = x(whichSNPs(incl));
        
        ll = 1/2 * (-(logdetDplusP - logdetP11) - x'*x - mmz*log(2*pi));
    end


end

