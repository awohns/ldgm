function [alt_mse] = computeAltMSE(P,R,whichSNPs)
%computeAltMSE computes alternative MSE between P and R
% P: precision matrix
% R: correlation matrix
% whichSNPs: rows/column indices of P that correspond to those of R

if islogical(whichSNPs); whichSNPs = find(whichSNPs); end
mm = length(P);
otherSNPs = setdiff(1:mm, whichSNPs);
RP = R * (P(whichSNPs,whichSNPs) - ...
    P(whichSNPs,otherSNPs) * (P(otherSNPs,otherSNPs) \ P(otherSNPs,whichSNPs)));
PR = RP';
alt_mse = mean((eye(length(whichSNPs)) - RP) .* (eye(length(whichSNPs)) - PR), 'all');
end

