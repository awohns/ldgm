function [P, objective_function_value, p_value] = LDPrecision(R,varargin)
% Calculates maximum likelihood precision matrix from correlation
% matrix R. R is a SNPs x SNPs matrix, possibly sparse, where edge (i,j)
% is the correlation between SNPs i and j. By default, LDPrecision
% estimates a precision matrix P whose edges correspond to the nonzero
% entries of R.  
% Optionally input arguments:
% P0: starting point for gradient descent
% graphical_model: specify a set of edges that are different from the
% nonzero entries of R. If p-value is desired, then this must be specified
% and R must be the full correlation matrix.
% convergence_tol: gradient descent convergence criterion
% max_steps: maximum no. gradient descent steps
% sample_size: affects objective_function_value and p_value but not P
% printstuff.

% Input data handling
p=inputParser;
addRequired(p, 'R', @(M)issymmetric(M) & isnumeric(M));
mm = length(R);
addParameter(p, 'convergence_tol', 1e-5, @isscalar);
addParameter(p, 'max_steps', 1e3, @isscalar);
addParameter(p, 'P0', speye(mm), @issymmetric);
addParameter(p, 'graphical_model', R~=0, @issymmetric);
addParameter(p, 'printstuff', false, @islogical);
addParameter(p, 'sample_size', 1, @isscalar);

parse(p,R,varargin{:});
P = p.Results.P0; % Precision matrix
G = p.Results.graphical_model;% Graphical model
if ~all(G(P~=0))
    error('Initial precision matrix conflicts with graphical model')
end
[ii,jj] = find(G);

% objective function
obj = @(precision)objFn(R, precision, p.Results.sample_size);
objective_function_value = obj(P);

% gradient descent
tic;
stepsize = 1;
for step = 1:p.Results.max_steps
    if p.Results.printstuff
        fprintf('%d ', step)
    end
    % gradient of objective function
    Pinv = sparseinv(P);
    gradient = sparse(ii,jj,2 * (R(G) - Pinv(G)),mm,mm);
    
    % line search to determine step size
    oldObj = objective_function_value;
    [P, stepsize, objective_function_value] = linesearch(P, objective_function_value, gradient, obj, stepsize);
    
    % convergence
    if abs(objective_function_value - oldObj) / abs(objective_function_value) < p.Results.convergence_tol
        break;
    end
end

if p.Results.printstuff
    fprintf('\n time = %f\n', toc)
end

if objective_function_value <= (1 - p.Results.convergence_tol) * oldObj
    warning('Gradient descent failed to converge in %d steps', p.Results.max_steps);
end

% p-value for whether graph fits data
if nargout > 2
    if ~issparse(R) && p.Results.sample_size ~= 1
        pval = chi2cdf(2 * (obj(inv(R)) - obj(P)), mm*(mm+1)/2 - (nnz(G) + mm)/2, 'upper');
    else
        warning('p_value output requested, but inputs might not be correct. Please specify a full covariance matrix R and the sample size that was used to compute it.')
        pval = [];
    end
end


    % objective function
    function val = objFn(R, omega, nn)
        [L, pp] = chol(omega);
        if pp == 0 % omega is positive definite
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

