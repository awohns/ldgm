function [estimate, steps] = h2MLE(alphahat,P,varargin)
%h2MLE computes maximum-likelihood heritability estimates using gradient
%descent and LDGMs
%   Inputs:
%   alphahat: effect-size estimates, units of sd(Y)/sd(X), or Z scores.
%   Should be a cell array with one cell per LD block.
%   P: LD graphacal model as a cell array.
%   nn: GWAS sample size
%   whichSNPs: which SNPs in the LDGM have nonmissing summary statistics.
%   Recommended to use as many SNPs as possible. All SNPs in the summary
%   statistics vectors should be in the LDGM.
%   annot: annotation matrices, one for each LD block, same number of SNPs
%   as alphahat
%   linkFn: function mapping from annotation vectors to per-SNP h2
%   linkFnGrad: derivative of linkFn
%   params: initial point for parameters
%   convergenceTol: terminates when objective function improves less than
%   this
%   maxReps: maximum number of steps to perform
%   minReps: starts checking for convergence after this number of steps
%   method: please use 'gradient'

% initialize
p=inputParser;
if iscell(P)
    mm = sum(cellfun(@length,P));
else
    mm = length(P);
end

addRequired(p, 'alphahat', @(x)isvector(x) || iscell(x))
addRequired(p, 'P', @(x)ismatrix(x) || iscell(x))
addOptional(p, 'nn', 1, @isscalar)
addOptional(p, 'whichSNPs', cellfun(@(x){true(size(x))},alphahat), @(x)isvector(x) || iscell(x))
addOptional(p, 'annot', ones(mm,1), @(x)size(x,1)==mm || iscell(x))
addOptional(p, 'linkFn', @(a,x)max(a*x,0), @(f)isa(f,'function_handle'))
addOptional(p, 'linkFnGrad', @(a,x)a.*(a*x>=0), @(f)isa(f,'function_handle'))
addOptional(p, 'params', [], @isvector)
addOptional(p, 'convergenceTol', 1e-6, @isscalar)
addOptional(p, 'maxReps', 1e4, @isscalar)
addOptional(p, 'minReps', 2, @isscalar)
addOptional(p, 'method', 'gradient', @isstr)
addOptional(p, 'noVarianceEstimate', 0, @isscalar)

parse(p, alphahat, P, varargin{:});

% turns p.Results.x into just x
fields = fieldnames(p.Results);
for k=1:numel(fields)
    line = sprintf('%s = p.Results.%s;', fields{k}, fields{k});
    eval(line);
end

noAnnot = size(annot{1},2);
noBlocks = length(annot);
if isempty(params)
    params = zeros(noAnnot,1);
end
noParams = length(params);
smallNumber = 1e-12;

objFn = @(params)-GWASlikelihood(alphahat,...
    cellfun(@(x){linkFn(x, params)}, annot),...
    P, nn, whichSNPs);

newObjVal = objFn(params);


gradient = 1e-9*randn(noParams,1);


method = p.Results.method;
if strcmp(method,'coordinate')
    minReps = minReps * noParams;
end
if 1 %strcmp(method,'coordinate') || strcmp(method,'weighted_gradient')
    
    stepSize = [1./sum(annot{1})'.^2/nn; ones(noParams-noAnnot,1)];
    %         stepSize = ones(noParams,1)/nn/mm^2;
    
elseif strcmp(method,'gradient')
    warning('This method seems to work poorly')
    stepSize = ones(noParams,1)/nn/mm^2;
else
    error('Valid options for optimization method are coordinate, weighted_gradient, gradient')
end
% dampening = 0;

allSteps=zeros(min(maxReps,1e6),noParams);
allValues=zeros(min(maxReps,1e6),1);
allStepSizes=allSteps;
allGradients=allSteps;

for rep=1:maxReps
    
    oldGradient=gradient;
    if strcmp(method,'gradient')
        gradient = 0;
        whichParam = 1:noParams;
        for block = 1:noBlocks
            sigmasq = linkFn(annot{block}, params);
            sigmasqGrad = linkFnGrad(annot{block}, params);
            gradient = gradient + ...
                GWASlikelihoodGradient(alphahat{block},1./sigmasq,P{block},nn,sigmasqGrad,whichSNPs{block})';
        end
    else
        gradient = zeros(noParams,1);
        if strcmp(method,'weighted_gradient')
            whichParam = 1:noParams;
        elseif strcmp(method,'coordinate')
            whichParam = mod(rep-1,noParams)+1;
        end
        for k = whichParam
            eps = zeros(noParams,1);
            eps(k) = smallNumber;
            gradient(k) = (objFn(params + eps) - newObjVal) ./ eps(k);
        end
        
    end
    
    % dampening to eliminate oscillations
    %     gradient = dampening * oldGradient + (1 - dampening) * gradient;
    %     grad_corr = dot(gradient,oldGradient)/sqrt(sum(gradient.^2)*sum(oldGradient.^2));
    %     dampening = max(0, min(0.99, dampening - grad_corr * 0.05));
    
    oldObjVal = newObjVal;
    
    % dealing with threshold effects
    if all(gradient(whichParam)==0)
        gradient(whichParam) = -1;
        stepSize(whichParam) = 1;
    end
    
    stepsize_factor = 1.5 * 4.^(1/(1 + floor((rep-1)/noParams)));
    
    
    [params, stepSize(whichParam), newObjVal] = linesearch(params, ...
        oldObjVal, gradient, objFn, stepSize(whichParam), stepsize_factor);
    
    allValues(rep)=newObjVal;
    if nargout > 1
        allSteps(rep,:)=params;
        allStepSizes(rep,:) = stepSize;
        allGradients(rep,:) = gradient;
    end
    if rep > minReps
        if allValues(rep-minReps) - newObjVal < convergenceTol * minReps
            break;
        end
    end
    
    %     converged = (newObjVal-oldObjVal)/oldObjVal < convergenceTol;
end

% Convergence report
if nargout > 1
    steps.reps = rep;
    steps.params = allSteps(1:rep,:);
    steps.obj = allValues(1:rep);
    steps.stepsize = allStepSizes(1:rep,:);
    steps.gradient = allGradients(1:rep,:);
end

% From paramaters to h2 estimates
h2Est = 0;annotSize=0;totalh2=0;
for block=1:noBlocks
    perSNPh2 = linkFn(annot{block}, params);
    h2Est = h2Est + sum(perSNPh2.*annot{block});
    annotSize = annotSize + sum(annot{block});
end

estimate.params = params;
estimate.h2 = h2Est;
estimate.annotSum = annotSize;

% Enrichment only calculated if first annotation is all-ones vector
if all(cellfun(@(a)all(a(:,1)==1),annot))
    estimate.enrichment = (h2Est./annotSize) / (h2Est(1)/annotSize(1));
end

% Compute covariance of the parameter estimates
if ~noVarianceEstimate
    smallNumber = 1e-6;
    gradient = 0;
    for block = 1:noBlocks
        sigmasq = linkFn(annot{block}, params);
        sigmasqGrad = linkFnGrad(annot{block}, params);
        gradient = gradient + ...
            GWASlikelihoodGradient(alphahat{block},1./sigmasq,P{block},nn,sigmasqGrad,whichSNPs{block})';
    end
    
    FI = zeros(noParams);
    for k = 1:noParams
        newParams = params;
        newParams(k) = params(k) * (1+smallNumber);
        gradient_k = 0;
        for block = 1:noBlocks
            sigmasq = linkFn(annot{block}, newParams);
            sigmasqGrad = linkFnGrad(annot{block}, newParams);
            gradient_k = gradient_k + ...
                GWASlikelihoodGradient(alphahat{block},1./sigmasq,P{block},nn,sigmasqGrad,whichSNPs{block})';
        end
        
        FI(:,k) = (gradient_k - gradient) / (smallNumber * params(k));
    end
    
    % Variance of h2 estimates
    dh2da = 0;
    for block=1:noBlocks
        dh2da = dh2da + linkFnGrad(annot{block}, params)'*annot{block};
        
    end
    h2Var = dh2da'/FI * dh2da;
    
    
    
    estimate.paramsVar = pinv(FI);
    estimate.FI = FI;
    estimate.h2Var = h2Var;
    
    if all(cellfun(@(a)all(a(:,1)==1),annot))
        estimate.enrichmentZscore = (h2Est./annotSize - h2Est(1)/annotSize(1)) ./ ...
            sqrt(h2Var(1,:)./annotSize.^2 + h2Var(1)./annotSize(1).^2 - 2*h2Var(1,:)*h2Var(1)/annotSize(1)./annotSize);
        estimate.enrichmentZscore(1) = 0;
    end
    
    
end

end

function [thetaNew, step, newObj] = linesearch(initTheta, initObj, grad, objFn, step, stepsize_factor)

oldObj = initObj;

if ~isreal(oldObj); error('objective function value should be real at initial point for line search'); end

step = step * stepsize_factor;

newObj = objFn(initTheta - step .* grad);


while newObj > initObj -  sum(step .* grad.^2) / stepsize_factor
    step = step / stepsize_factor;
    newObj = objFn(initTheta - step .* grad);
end

thetaNew = initTheta - step .* grad;
end

