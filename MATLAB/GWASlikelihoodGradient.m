function [grad] = GWASlikelihoodGradient(alphaHat, sigmasq, P, nn, delSigmaDelA, whichSNPs, fixedIntercept)
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

if nargin < 7
    fixedIntercept = true;
end

if iscell(P) % handle cell-array-valued inputs
    if nargin < 6
        whichSNPs = cellfun(@(x)true(size(x),alphaHat),'UniformOutput', false);
    end
    assert(iscell(alphaHat) & iscell(whichSNPs) & iscell(sigmasq))
    grad = cellfun(@(a,p,w)GWASlikelihoodGradient(a,s,p,nn,delSigmaDelA,w,fixedIntercept),...
        alphaHat,sigmasq,P,whichSNPs, 'UniformOutput', false);
else
    mm = length(P);
    
    % M = Sigma + P/nn is the covariance matrix of betaHat = P/P11 * alphaHat
    M = zeros(mm,1);
    M(whichSNPs) = sigmasq;
    M = P/nn + speye(mm).*M;
    Minv = sparseinv(M);
    MinvDiag = diag(Minv);
    
    % betahat = P/P11 * alphaHat
    betahat = precisionMultiply(P, alphaHat, whichSNPs);
    b = precisionDivide(M, betahat, whichSNPs);
    
    % two terms of the log-likelihood (quadratic term + log-determinant
    % term)
    d = -sum(delSigmaDelA .* b.^2);
    delLogDetP = sum(delSigmaDelA .* MinvDiag(whichSNPs));
    
    grad = - 1/2 * (- delLogDetP - d);
    
    % gradient of minus log-likelihood wrt 1/nn
    if ~fixedIntercept
        c = precisionMultiply(P,b,whichSNPs);
        nGrad = 1/2 * (sum(nonzeros(Minv.*P)) - sum(b.*c));
        grad(end+1) = nGrad;
    end
end
end

