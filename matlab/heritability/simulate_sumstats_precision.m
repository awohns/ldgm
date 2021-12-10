function [alphaHat,beta,alpha,sigmasq] = simulate_sumstats_precision(P, nn, sigmasqSupport, varargin)
% Simulates summary statistics alphahat from specified prior distribution
% Input parameters:
% P: LD precision matrix
% nn: GWAS sample size
% sigmasqSupport: the possible sigma^2 values for each annotation (which will
% be added across annotations by default).
% Optional input parameters:
% sigmasqPrior: effect-size priors for each annotation, as a noAnnot x noCpts matrix. Each
% row encodes the probability that a SNP from each annotation will be
% assigned variance corresponding to each element of sigmasqSupport. To
% give every SNP an annotation-specific variance, set sigmasqPrior =
% eye(noAnnot) and sigamsqSupport to the desired values
% annot: mm x noAnnot matrix of annotation values for each SNP
% blocks: if P is block diagonal, specify blocks to speed things up, as
% indices in a cell array
% linkFn: link function to convert from annotation-specific variances into
% a per-SNP effect size variance in the presence of overlapping
% annotations.


p=inputParser;
mm = length(P);
addRequired(p, 'P', @ismatrix)
addRequired(p, 'nn', @isscalar)
addRequired(p, 'sigmasqSupport', @isnumeric)
addOptional(p, 'sigmasqPrior', 1, @(x)all(sum(x,2)==1))
addOptional(p, 'annot', ones(mm,1), @(x)size(x,1)==mm)
addOptional(p, 'blocks', {1:mm}, @iscell)
addOptional(p, 'linkFn', @(x)max(x,0), @(f)isa(f,'function_handle'))

parse(p,P, nn, sigmasqSupport,varargin{:});

sigmasqPrior = p.Results.sigmasqPrior;
annot = p.Results.annot;
blocks = p.Results.blocks;
linkFn = p.Results.linkFn;

sigmasq = zeros(mm,1);
noCpts = length(sigmasqSupport);
noAnnot = size(annot,2);

if noCpts > 1
    annotSigmasq = zeros(mm,1);
    for k = 1:noAnnot
        annotSigmasq(:) = randsample(sigmasqSupport,mm,true,sigmasqPrior(k,:));
        sigmasq = sigmasq + annot(:,k) .* annotSigmasq;
    end
else
    sigmasq = sigmasqSupport;
end

sigmasq = linkFn(sigmasq);
if any(sigmasq < 0)
    error('Negative effect-size variances produced; check link function gives nonnegative output')
end

beta = randn(mm,1) .* sqrt(sigmasq);

alpha = zeros(mm,1);
alphaHat = alpha;
for b = 1:length(blocks)
    cholP = chol(P(blocks{b},blocks{b}));
    alpha(blocks{b}) = P(blocks{b},blocks{b}) \ beta(blocks{b});
    alphaHat(blocks{b}) = alpha(blocks{b}) + ...
        cholP \ (randn(length(blocks{b}),1) / sqrt(nn));
end


end

