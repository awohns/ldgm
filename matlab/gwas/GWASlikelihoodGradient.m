function [grad,M] = GWASlikelihoodGradient(alphahat,tau,P,nn,delSigmadelA,whichSNPs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


mm = length(P);
if nargin < 6
    whichSNPs = true(mm,1);
end
M = zeros(mm,1);
M(whichSNPs) = 1./tau;
M = P/nn + speye(mm).*M;

MinvDiag = diag(sparseinv(M));

% P/P11 * alphaHat
betahat = zeros(mm,1);
betahat(whichSNPs) = P(whichSNPs,whichSNPs) * alphahat - ...
    P(whichSNPs,~whichSNPs) * (P(~whichSNPs,~whichSNPs) \ ...
    (P(~whichSNPs,whichSNPs)*alphahat));

b = M \ betahat;

d = -sum(delSigmadelA .* b(whichSNPs).^2);
delLogDetP = sum(delSigmadelA .* MinvDiag(whichSNPs));

grad = - 1/2 * (- delLogDetP - d);


end

