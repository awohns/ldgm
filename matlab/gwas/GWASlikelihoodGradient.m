function [grad, M] = GWASlikelihoodGradient(alphahat, tau, P, nn, delSigmadelA, whichSNPs, fixedIntercept)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


mm = length(P);
if nargin < 6
    whichSNPs = true(mm,1);
end
if nargin < 7
    fixedIntercept = true;
end
M = zeros(mm,1);
M(whichSNPs) = 1./tau;
M = P/nn + speye(mm).*M;

Minv = sparseinv(M);
MinvDiag = diag(Minv);

% P/P11 * alphaHat
betahat = zeros(mm,1);
betahat(whichSNPs) = P(whichSNPs,whichSNPs) * alphahat - ...
    P(whichSNPs,~whichSNPs) * (P(~whichSNPs,~whichSNPs) \ ...
    (P(~whichSNPs,whichSNPs)*alphahat));

b = M \ betahat;

d = -sum(delSigmadelA .* b(whichSNPs).^2);
delLogDetP = sum(delSigmadelA .* MinvDiag(whichSNPs));

grad = - 1/2 * (- delLogDetP - d);

% gradient of minus log-likelihood wrt 1/nn
if ~fixedIntercept
    b = b(whichSNPs);
    c = P(whichSNPs,whichSNPs) * b - ...
        P(whichSNPs,~whichSNPs) * (P(~whichSNPs,~whichSNPs) \ ...
        (P(~whichSNPs,whichSNPs)*b));
    nGrad = 1/2 * (sum(nonzeros(Minv.*P)) - sum(b.*c));
    grad(end+1) = nGrad;
end
end

