function [params,newObjVal,allSteps,allValues] = h2MLE(alphahat,P,varargin)
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
addOptional(p, 'whichSNPs', 1:mm, @(x)isvector(x) || iscell(x))
addOptional(p, 'annot', ones(mm,1), @(x)size(x,1)==mm || iscell(x))
addOptional(p, 'blocks', {1:mm}, @iscell)
addOptional(p, 'linkFn', @(a,x)max(a*x,0), @(f)isa(f,'function_handle'))
addOptional(p, 'params', [], @isvector)
addOptional(p, 'convergenceTol', 1e-6, @isscalar)
addOptional(p, 'maxReps', 1e4, @isscalar)
addOptional(p, 'minReps', 2, @isscalar)
addOptional(p, 'method', 'coordinate', @isstr)

parse(p, alphahat, P, varargin{:});

whichSNPs = p.Results.whichSNPs;
% if min(sum(whichSNPs),length(whichSNPs)) ~= length(alphahat)
%     error('whichSNPs should specify indices corresponding to alphahat')
% end

nn = p.Results.nn;

annot = p.Results.annot;
if ~iscell(annot); annot = {annot}; end
noAnnot = size(annot{1},2);

blocks = p.Results.blocks;

linkFn = p.Results.linkFn;

params = p.Results.params;
if isempty(params)
    params = zeros(noAnnot,1);
end
noParams = length(params);

convergenceTol = p.Results.convergenceTol;

minReps = p.Results.minReps;

maxReps = p.Results.maxReps;

smallNumber = 1e-12;

whichSNPs = p.Results.whichSNPs;

objFn = @(params)-GWASlikelihoodMissingness(alphahat,...
    cellfun(@(x){linkFn(x, params)}, annot),...
    P, nn, whichSNPs, blocks);

newObjVal = objFn(params);


gradient = 1e-9*randn(noParams,1);


method = p.Results.method;
if strcmp(method,'coordinate')
    minReps = minReps * noParams;
end
if strcmp(method,'coordinate') || strcmp(method,'weighted_gradient')
    
    stepSize = [1./sum(annot)'.^2/nn; ones(noParams-noAnnot,1)];
elseif strcmp(method,'gradient')
    warning('This method seems to work poorly')
    stepSize = ones(noParams,1)/nn/mm^2;
else
    error('Valid options for optimization method are coordinate, weighted_gradient, gradient')
end
% dampening = 0;

allSteps=zeros(maxReps,noParams);allValues=zeros(maxReps,1);
for rep=1:maxReps
    if mod(rep-1,minReps)==0
        %         disp(params')
    end
    %     sigmasq = linkFn(annot, params);
    %     sigmasqGrad = linkFnGrad(sigmasq).*annot;
    
    oldGradient=gradient;
    if 0
        warning('this isnt working yet')
        gradient = - GWASlikelihoodGradient(alphahat,1./sigmasq,P,nn,sigmasqGrad);
    else
        gradient = zeros(noParams,1);
        if strcmp(method,'gradient') || strcmp(method,'weighted_gradient')
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
    
    allSteps(rep,:)=params;
    allValues(rep)=newObjVal;
    
    
    
    if rep > minReps
        if allValues(rep-minReps) - newObjVal < convergenceTol * minReps
            break;
        end
    end
    
    %     converged = (newObjVal-oldObjVal)/oldObjVal < convergenceTol;
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

