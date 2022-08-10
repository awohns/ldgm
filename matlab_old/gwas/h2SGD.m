function [params,newObjVal,allSteps,allValues] = h2SGD(alphahat,P,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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
addOptional(p, 'blocks', {1:mm}, @iscell)
addOptional(p, 'linkFn', @(a,x)max(a*x,0), @(f)isa(f,'function_handle'))
addOptional(p, 'linkFnGrad', @(a,x)a.*(a*x>=0), @(f)isa(f,'function_handle'))
addOptional(p, 'params', [], @isvector)
addOptional(p, 'convergenceTol', 1e-6, @isscalar)
addOptional(p, 'maxReps', 1e4, @isscalar)
addOptional(p, 'minReps', 2, @isscalar)

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

% Initial step size
stepSize = 1/nn/mm;% * 0.9.^(1:maxReps*noBlocks);%[1./sum(annot{1})'.^2/nn; ones(noParams-noAnnot,1)];
stepSizeFactor = 2;

% Dampening/momentum
dampening = 0;
gradient = randn(1, noParams);

allSteps=zeros(maxReps,noParams);allValues=zeros(maxReps,1);
for rep=1:maxReps*noBlocks
    
    % random order for SGD
    if mod(rep,noBlocks)==1
        blockOrder = randperm(noBlocks);
    end
    
    % current block for SGD
    whichBlock = blockOrder(mod(rep-1,noBlocks) + 1);
    
    % gradient of the current block
    oldGradient = gradient;
    sigmasq = linkFn(annot{whichBlock}, params);
    sigmasqGrad = linkFnGrad(annot{whichBlock}, params);
    gradient = GWASlikelihoodGradient(alphahat{whichBlock}, ...
        1./sigmasq, P{whichBlock}, nn, sigmasqGrad)'; %, whichSNPs{whichBlock}
    
    % dampening to eliminate oscillations
%     gradient = dampening * oldGradient + (1 - dampening) * gradient;
%     grad_corr = dot(gradient,oldGradient)/sqrt(sum(gradient.^2)*sum(oldGradient.^2));
%     dampening = max(0, min(0.99, dampening - grad_corr * 0.05));
    
    objFn = @(params)-GWASlikelihood(alphahat{whichBlock},...
            linkFn(annot{whichBlock}, params),...
            P{whichBlock}, nn, whichSNPs{whichBlock});
        
    % Every pass through the data, update the step size
    if rep <= 3 || mod(rep,noBlocks)==1
        oldObjVal = objFn(params);
        
        [params, step, newObjVal] = linesearch(params, ...
            oldObjVal, gradient, objFn, stepSize./sum(annot{whichBlock})', stepSizeFactor);
        stepSize = step .* sum(annot{whichBlock})';
    else
%         stepSize = stepSize*0.9;
        params = params - (stepSize./sum(annot{whichBlock})') .* gradient;
%         newObjVal = objFn(params);
    end
    
    allSteps(rep,:)=params;
    allValues(rep)=newObjVal;
%     allStepSizes(rep,:)=stepSize./sum(annot{whichBlock});
    
    
    if mod(rep,noBlocks) ==0 && rep > noBlocks * (minReps + 1)
        passes = rep / noBlocks;
        oldObj = sum(allValues(noBlocks * (passes - (minReps + 1)) + 1 : ...
            noBlocks * (passes - minReps)));
        newObj = sum(allValues(noBlocks * (passes - 1) + 1 : ...
            noBlocks * passes));
        if (oldObj - newObj)/abs(oldObj) < convergenceTol * minReps
            break;
        end
    end
    
    %     converged = (newObjVal-oldObjVal)/oldObjVal < convergenceTol;
end

allSteps = allSteps(1:rep,:);

end

function [thetaNew, step, newObj] = linesearch(initTheta, initObj, grad, objFn, step, stepsize_factor)

oldObj = initObj;

if ~isreal(oldObj); error('objective function value should be real at initial point for line search'); end

newObj = objFn(initTheta - step .* grad);

if newObj < initObj -  sum(step .* grad.^2) / stepsize_factor
    oldObj = initObj;
    while newObj < oldObj
        step = step * stepsize_factor;
        oldObj = newObj;
        newObj = objFn(initTheta - step .* grad);
    end
    step = step / stepsize_factor;
    newObj = objFn(initTheta - step .* grad);
else
    while newObj > initObj -  sum(step .* grad.^2) / stepsize_factor
        step = step / stepsize_factor;
        newObj = objFn(initTheta - step .* grad);
    end
end
thetaNew = initTheta - step .* grad;
end

