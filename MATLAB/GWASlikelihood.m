function [logLikelihood, logLikelihoodBlocks] = GWASlikelihood(Z_deriv_allele,sigmasq,P,nn,whichIndices,intercept)
% GWASlikelihood computes likelihood of the GWAS sumstats, Z_deriv_allele,
% under a gaussian model:
%                   beta ~ MVN(mu,diag(sigmasq))
%                   Z|beta ~ MVN(sqrt(nn)*R*beta, R)
%                   inv(R) = P.
%
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

if nargin < 6
    intercept = 1;
end
assert(isscalar(intercept))
if intercept < 0
    logLikelihood = -inf;
    return;
end

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
    logLikelihoodBlocks = cellfun(@(a,s,p,w)likelihoodFn(a,s,p,nn,w,intercept), ...
        Z_deriv_allele, sigmasq, P, whichIndices);

    logLikelihood = sum(logLikelihoodBlocks);
else
    logLikelihood = likelihoodFn(Z_deriv_allele,sigmasq,P,nn,whichIndices,intercept);
    logLikelihoodBlocks = logLikelihood;
end

    function ll = likelihoodFn(Z,sigmasq,P,nn,whichSNPs,intercept)

        if ~islogical(whichSNPs)
            whichSNPs = unfind(whichSNPs,length(P));
        end

        assert(all(sigmasq>=0),'sigmasq should be nonnegative')

        % handling SNPs missing from P
        incl = diag(P)~=0;
        assert(all(incl(whichSNPs)))
        mm = sum(incl);
        P = P(incl,incl);
        whichSNPs = whichSNPs(incl);

        % inv(P)(whichSNPs,whichSNPs) * Z
        x = precisionMultiply(P,Z,whichSNPs);

        % M == E(xx')
        M = sparse(find(whichSNPs), find(whichSNPs), nn*sigmasq, mm, mm);
        M = M + intercept * P;
        A = chol(M);

        % log|P/P11 + nn*diag(sigmasq)| == log|M| - log|P11|
        logdetM = 2*sum(log(diag(A)));
        logdetP11 = 2*sum(log(diag(chol(P(~whichSNPs,~whichSNPs)))));

        % x'*M\x == w'*w
        y = zeros(mm,1);
        y(whichSNPs) = x;
        w = A' \ y;

        ll = 1/2 * (-(logdetM - logdetP11) - w'*w - sum(whichSNPs)*log(2*pi));
    end


end

