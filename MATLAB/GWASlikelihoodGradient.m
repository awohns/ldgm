function [grad] = GWASlikelihoodGradient(Z, sigmasq, P, nn, delSigmaDelA, whichSNPs, intercept, fixedIntercept)
% GWASlikelihoodGradient computes gradient of the likelihood of the GWAS 
% sumstats alphaHat under a gaussian model:
%                   beta ~ MVN(mu,diag(sigmasq))
%                   alphaHat|beta ~ MVN(R*beta, R/nn)
%                   inv(R) = P.
% 
% This function requires the sparseinv() function from suitesparse:
% https://people.engr.tamu.edu/davis/suitesparse.html
% 
% The GWAS SNPs in alphahat should be a subset of those in P, and in the
% same order; boolean vector whichSNPs should be true for rows/columns of P
% corresponding to one of these SNPs, false elsewhere.
%
% sigmasq should be the same size as alphahat, such that missing
% SNPs are modeled as having zero effect-size variance.
%
% Argument delSigmaDelA is the gradient of sigmasq w.r.t. some parameters A.
% GWASlikelihoodGradient computes the gradient of the likelihood w.r.t.
% A, i.e.:
% grad_A loglikelihood = delSigmaDelA * grad_sigmasq loglikelihood 
%
% Units of alphaHat should be either (1) normalized effect size estimates,
% ie sample correlations, or (2) Z scores. In case (1), nn should be the
% number of individuals in the GWAS. In case (2), nn should be 1, and the
% vector of sigmasq values should be multiplied by nn.
% Optionally, arguments alphaHat, sigmasq, P, whichSNPs can be specified as
% cell arrays with one cell for each LD block.
%
% Optionally, arguments alphaHat, sigmasq, P, whichSNPs can be specified as
% cell arrays with one cell for each LD block.
%
% Optionally, if it is desired to specify the mean of beta, ie
%                   beta ~ MVN(mu, diag(sigmasq))
% call GWASlikelihoodGradient( alphahat - P \ mu, ...)
%
% Optionally, if it is desired to model an unknown amount of population
% stratification/relatedness, or if the sample size is unknown, specify
% fixedIntercept = true. This modifies the model as:
% alphaHat|beta ~ MVN(R*beta, R * (1/nn + a))
% where a==0, and it computes the derivative of the log-likelihood with 
% respect to a; this will be the last element of grad.

assert(isscalar(intercept) && all(intercept>=0,'all'))

if iscell(P) % handle cell-array-valued inputs
    if nargin < 5
        whichSNPs = cellfun(@(x)true(size(x),Z),'UniformOutput', false);
    end
    assert(iscell(Z) & iscell(whichSNPs) & iscell(sigmasq))
    grad = cellfun(@(a,s,p,dS,w)GWASlikelihoodGradient(a,s,p,nn,dS,w,intercept,fixedIntercept),...
        Z,sigmasq,P,delSigmaDelA,whichSNPs, 'UniformOutput', false);
else
    
    % handle missing rows/cols of P
    if ~islogical(whichSNPs)
        whichSNPs = unfind(whichSNPs,length(P));
    end

    % handling SNPs missing from P
    incl = diag(P)~=0;
    assert(all(incl(whichSNPs)))
    mm = sum(incl);
    P = P(incl,incl);
    whichSNPs = whichSNPs(incl);
    mm0 = sum(~whichSNPs); % no. missing SNPs from sumstats
    
    % M == E(xx')
    M = sparse(find(whichSNPs), find(whichSNPs), nn*sigmasq, mm, mm);
    M = M + intercept * P;

    % sparse inverse subset of M
    Minv = sparseinv(M);
    MinvDiag = diag(Minv);
    
    % x = P/P11 * z
    x = precisionMultiply(P, Z, whichSNPs);
    b = precisionDivide(M, x, whichSNPs);
    
    % two terms of the log-likelihood (quadratic term + log-determinant
    % term)
    delwtw = -nn * sum(delSigmaDelA .* b.^2);
    delLogDetM = nn * sum(delSigmaDelA .* MinvDiag(whichSNPs));
    
    grad = - 1/2 * (delLogDetM + delwtw);
    
    % gradient of minus log-likelihood wrt 1/nn
    if ~fixedIntercept
        c = precisionMultiply(P,b,whichSNPs);
        grad(end+1) = -1/2 * (sum(nonzeros(Minv.*P)) - sum(b.*c) - mm0 / intercept);
    end
end

end

