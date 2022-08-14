function r2 = r2PGS(true_beta_perallele, estimated_beta_perallele, P, whichIndices, AF)
%r2PGS calculates the PGS accuracy, i.e. the squared correlation between
%X*true_beta_perSD and X*estimated_beta_perSD, using the precision matrix.
% 
% Input arguments:
% true_beta_perSD: 'true' causal effect sizes in per-allele units, as a cell 
% array of size number-of-blocks by number-of-populations
% 
% estimated_beta_perallele: estimated causal effect sizes, as a cell array of
% the same size, or as a number-of-blocks by 1 cell array (in which case
% the same effect size estimates are used for every population)
% 
% P: LDGM precision matrices, as a cell array with each cell containing a 
% precision matrix, for the target population (as opposed to the training
% population)
% 
% whichIndices: output from mergesnplists, encoding which indices
% (rows/columns of the LDGM precision matrices) have corresponding SNPs in
% each summary statistic file respectively. Same size as P.
% 
% AF: allele frequencies, as a number-of-blocks by number-of-populations
% cell array
% 
% Output arguments:
% r2: the squared correlation between the predictions (this is usually
% larger than the squared correlation between the true and predicted
% weights themselves)

[noBlocks,noPopns] = size(P);
if size(estimated_beta_perallele,2) == 1
    estimated_beta_perallele = repmat(estimated_beta_perallele,1,noPopns);
end

% convert per-allele to per-SD units
if nargin < 5
    warning('Allele frequencies not specified; assuming that input effect sizes are in per-SD units, which is not the recommended usage')
else
    conversion_fn = @(beta_perallele, af){beta_perallele .* sqrt(2*af.*(1-af))};
    true_beta_perSD = cellfun(conversion_fn,true_beta_perallele,AF);
    estimated_beta_perSD = cellfun(conversion_fn,estimated_beta_perallele,AF);
end

% quadratic function
qf = @(beta1,beta2,P,whichIndices)beta1(whichIndices)' * precisionDivide(P, beta2(whichIndices), whichIndices);

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

