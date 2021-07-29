function [P, pval] = LDPrecision(R,G,nn,convergenceTol,P0,printstuff)
% Calculates maximum likelihood precision matrix from genotype matrix X and graphical model
% G. X should be an N x M matrix and G a (possibly sparse) M x M
% symmetric logical-valued matrix. P is an M x M precision
% matrix with the same entries as G. pval tests whether G fits the data.
% reps is the number of gradient descend steps.
% P0 is initial guess for gradient descent.

[mm] = length(G);
if ~issymmetric(G) || ~islogical(G)
    error('G should be an undirected adjacency matrix with logical entries')
end

if any(size(R)~=size(G))
    error('R and G should agree in size')
end

% max number of steps for gradient descent
maxReps = 100;

% convergence tolerance

% initial guess for gradient descent: take full sample precision matrix and
% zero out non-edges
[ii,jj] = find(G);
if nargin < 5
    P = speye(mm);
else
    P = P0;
end
if nargin < 4
    convergenceTol = 1e-6;
end
if nargin < 3
    if nargout > 1
        error('To calculate p-values, sample size must be specified')
    else
        nn = 1;
    end
end
if nargin < 6
    printstuff = 1;%mm >=1000;
end

% objective function
obj = @objFn;
currentObj = obj(P);
tic;

stepsize = 1;
for step = 1:maxReps
    if printstuff
        fprintf('%d ', step)
    end
    % gradient of objective function
    Pinv = sparseinv(P);
    gradient = sparse(ii,jj,2 * (R(G) - Pinv(G)),mm,mm);
    
    % line search to determine step size
    oldObj = currentObj;
    [P, stepsize, currentObj] = linesearch(P, currentObj, gradient, obj, stepsize);
    
    % convergence
    if currentObj > (1 - convergenceTol) * oldObj
        break;
    end
end
if printstuff
    fprintf('\n time = %f\n', toc)
end
if currentObj <= (1 - convergenceTol) * oldObj
    warning('Gradient descent failed to converge in %d steps', maxReps);
end

% p-value for whether graph fits data
if nargout > 1
    pval = chi2cdf(2 * (obj(inv(R)) - obj(P)), mm*(mm+1)/2 - (nnz(G) + mm)/2, 'upper');
end



    function val = objFn(omega)
        [L, p] = chol(omega);
        if p == 0 % omega is positive definite
            val =  -0.5 * nn * (sum(2*log(diag(L))) - dot(nonzeros(omega), R(omega~=0)));
        else
            val = inf;
        end
    end
end

% line search to find optimal step size
function [thetaNew, step, newObj] = linesearch(initTheta, initObj, grad, objFn, step)

oldObj = initObj;
stepsize_factor = 2;
if ~isreal(oldObj); error('objective function value should be real at initial point for line search'); end

newObj = objFn(initTheta - step * grad);
if newObj < oldObj && isreal(newObj)
    while newObj < oldObj && isreal(newObj)
        
        step = step * stepsize_factor;
        oldObj = newObj;
        newObj = objFn(initTheta - step * grad);
        
        
    end
    step = step / stepsize_factor;
    newObj = oldObj;
end

while newObj > initObj || ~isreal(newObj)
    step = step / stepsize_factor;
    newObj = objFn(initTheta - step * grad);
end

thetaNew = initTheta - step * grad;
end

