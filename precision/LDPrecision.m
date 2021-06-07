function [P, pval] = LDPrecision(X,G,reps,P0)
% Calculates maximum likelihood precision matrix from genotype matrix X and graphical model
% G. X should be an N x M matrix and G a (possibly sparse) M x M
% symmetric logical-valued matrix. P is an M x M precision
% matrix with the same entries as G. pval tests whether G fits the data.
% reps is the number of gradient descend steps.
% P0 is initial guess for gradient descent.

[nn] = length(G);
if ~issymmetric(G) || ~islogical(G)
    error('G should be an undirected adjacency matrix with logical entries')
end

if size(X,2) ~= nn
    error('columns of X should correspond to rows/columns of G')
end

% number of steps for gradient descent
if nargin < 3
    reps = 100;
end

% initial guess for gradient descent: take full sample precision matrix and
% zero out non-edges
R = cov(X);
[ii,jj] = find(G);
edges = find(G); % linear indices
if nargin < 4
    P = inv(R);
    P = sparse(ii,jj,P(edges));
    
    % P should be positive definite for gradient to compute; ridge
    % regularize until it is.
    lambda = .1;
    while det(P) <= 0
        P = lambda * P + (1-lambda)*speye(nn);
    end
else
    P = P0;
end


% objective function
obj = @(Om) -0.5 * size(X,1) * (log(det(Om)) - sum(sum(Om.*R)));

stepsize = 1e-6;
for step = 1:reps
    % gradient of objective function
    Pinv = inv(P);
    gradient = sparse(ii,jj,2 * (R(edges) - Pinv(edges)),nn,nn);
    
    % line search to determine step size
    [P, stepsize] = linesearch(P, gradient, obj, stepsize);
end

% p-value for whether graph fits data
if nargout > 1
    ObjFn = @(Om) 0.5 * size(X,1) * (log(det(Om)) - sum(sum(Om.*R)));

    pval = chi2cdf(2 * (ObjFn(inv(R)) - ObjFn(P)), nn*(nn+1)/2 - (nnz(G) + nn)/2, 'upper');
end

end

% line search to find optimal step size
function [thetaNew, step, newObj] = linesearch(thetaInit, grad, objFn, stepInit)
oldObj = objFn(thetaInit);
if ~isreal(oldObj); error('objective function value should be real at initial point for line search'); end
step = stepInit;
newObj = objFn(thetaInit - step * grad);
while newObj < oldObj && isreal(newObj)
    step = step * 2;
    oldObj = newObj;
    newObj = objFn(thetaInit - step * grad);
end
step = step / 2;

oldObj = objFn(thetaInit);
while newObj > oldObj || ~isreal(newObj)
    step = step / 2;
    newObj = objFn(thetaInit - step * grad);
end

thetaNew = thetaInit - step * grad;
end

