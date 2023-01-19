function r2 = r2PGS(trueBeta, estimatedBeta, P, whichIndices, AF, useLdlChol)
%r2PGS calculates the PGS accuracy, i.e. the squared correlation between
%X*true_beta_perSD and X*estimated_beta_perSD, using the precision matrix.
% 
% Input arguments:
% trueBeta: 'true' causal effect sizes, as a cell array of size 
% number-of-blocks by number-of-populations. If AF is specified, trueBeta
% should be specified in per-allele units; this will be the most convenient
% choice when making comparisons across ancestry groups (with
% ancestry-specific LD + AF). If AF is empty, trueBeta should be specified
% in per-SD units.
% 
% estimatedBeta: estimated causal effect sizes, as a cell array of
% the same size, or as a number-of-blocks by 1 cell array (in which case
% the same effect size estimates are used for every population).
% If AF is specified, estimatedBeta
% should be specified in per-allele units; this will be the most convenient
% choice when making comparisons across ancestry groups (with
% ancestry-specific LD + AF). If AF is empty, trueBeta should be specified
% in per-SD units.
% 
% P: LDGM precision matrices, as a cell array with each cell containing a 
% precision matrix, for the target population (as opposed to the training
% population). If useLdlChol == true, then P should instead be a cell array
% of Cholesky factors computed using ldlchol.
% 
% whichIndices: output from mergesnplists, encoding which indices
% (rows/columns of the LDGM precision matrices) have corresponding SNPs in
% each summary statistic file respectively. Same size as P.
% 
% AF: allele frequencies, as a number-of-blocks by number-of-populations
% cell array. Can also be empty. If nonempty, trueBeta and estimatedBeta
% should be specified in per-allele units; if empty, they should be in
% per-SD units.
% 
% useLdlChol: whether Cholesky factors have been pre-computed using
% ldlchol(). This option leads to significant speedups. To precompute the
% cholesky factors, apply the following code to each LD block:
%   nz = any(P); % indices not missing from precision matrix
%   A = ldlchol(P(nz,nz)); % Cholesky factor
%   whichIndicesNz = lift(whichIndices, find(nz)); % indices into A
%   trueBeta = trueBeta(nz);
%   estimatedBeta = estimatedBeta(nz);
%   r2 = r2PGS(trueBeta, estimatedBeta, A, whichIndicesNz, AF, true);
% 
% Output arguments:
% r2: the squared correlation between the predictions (this is usually
% larger than the squared correlation between the true and predicted
% weights themselves)

[noBlocks,noPopns] = size(P);

if isempty(whichIndices)
    whichIndices = cellfun(@any, P, 'UniformOutput', false);
else
    error('whichIndices input argument is now required to be empty for this function to avoid misleading results')
end
if size(estimatedBeta,2) == 1
    estimatedBeta = repmat(estimatedBeta,1,noPopns);
end
if size(trueBeta,2) == 1
    trueBeta = repmat(trueBeta,1,noPopns);
end

assert(all(cellfun(@(a,b,c)length(a)==length(b) && length(b)==length(c),...
    P,trueBeta,estimatedBeta)),...
    'P, trueBeta, and estimatedBeta should be same-sized cell arrays with entries having the same length');

% convert per-allele to per-SD units
if nargin < 5 || isempty(AF)
    if noPopns > 1
        warning('Allele frequencies are not specified, and multiple populations are specified. This is unusual, since different populations have different allele frequencies.')
    end
    true_beta_perSD = trueBeta;
    estimated_beta_perSD = estimatedBeta;
else
    conversion_fn = @(beta_perallele, af){beta_perallele .* sqrt(2*af.*(1-af))};
    true_beta_perSD = cellfun(conversion_fn,trueBeta,AF);
    estimated_beta_perSD = cellfun(conversion_fn,estimatedBeta,AF);
end

if nargin < 6
    useLdlChol = false;
end



% quadratic function
qf = @(beta1,beta2,P,whichIndices)beta1(whichIndices)' * precisionDivide(P, beta2(whichIndices), whichIndices, useLdlChol);

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

