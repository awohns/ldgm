function r2 = r2PGS(true_beta_perSD, estimated_beta_perSD, P, whichIndices)
%r2PGS calculates the PGS accuracy, i.e. the squared correlation between
%X*true_beta_perSD and X*estimated_beta_perSD, using the precision matrix.
% 
% Input arguments:
% true_beta_perSD: 'true' causal effect sizes in per-SD units, as a cell 
% array of size number-of-blocks by number-of-populations
% 
% estimated_beta_perSD: estimated causal effect sizes, as a cell array of
% the same size
% 
% P: LDGM precision matrices, as a cell array with each cell containing a 
% precision matrix
% 
% whichIndices: output from mergesnplists, encoding which indices
% (rows/columns of the LDGM precision matrices) have corresponding SNPs in
% each summary statistic file respectively. Same size as P.
% 
% Output arguments:
% r2: the squared correlation between the predictions (this is usually
% larger than the squared correlation between the true and predicted
% weights themselves)

% quadratic function
qf = @(beta1,beta2,P,whichIndices)beta1(whichIndices)' * precisionDivide(P, beta2(whichIndices), whichIndices);

[noBlocks,noPopns] = size(P);
xtx = zeros(1,noPopns);
xty = zeros(1,noPopns);
yty = zeros(1,noPopns);

for block = 1:noBlocks
    for popn = 1:noPopns
        xtx(popn) = xtx(popn) + qf(true_beta_perSD{block,popn},...
            true_beta_perSD{block,popn},P{block,popn},whichIndices{block,popn});
        xty(popn) = xty(popn) + qf(true_beta_perSD{block,popn},...
            estimated_beta_perSD{block,popn},P{block,popn},whichIndices{block,popn});
        yty(popn) = yty(popn) + qf(estimated_beta_perSD{block,popn},...
            estimated_beta_perSD{block,popn},P{block,popn},whichIndices{block,popn});
    
    end
    
end

r2 = xty.^2 ./ (xtx .* yty);

end

