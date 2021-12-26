function [P, objective_function_value, converged, data] = LDPrecision(R,varargin)
% Calculates maximum likelihood precision matrix from correlation
% matrix R. R is a SNPs x SNPs matrix, possibly sparse, where edge (i,j)
% is the correlation between SNPs i and j. By default, LDPrecision
% estimates a precision matrix P whose edges correspond to the nonzero
% entries of R.
% Optional input arguments:
% P0: starting point for gradient descent
% graphical_model: specify a set of edges that are different from the
% nonzero entries of R. If p-value is desired, then this must be specified
% and R must be the full correlation matrix.
% convergence_tol: gradient descent convergence criterion
% max_steps: maximum no. gradient descent steps
% num_steps_check_convergence: checks for convergence every this many steps
% printstuff: 0, no output; 1, some output; 2, verbose
% lambda: l1 penalty to increase sparsity of P.
% Output arguments:
% objective_function_value: final objective function value
% converged: whether it terminated at convergence criterion (1) or max
% steps (0)
% data: contains information from each step of gradient descent

% Input data handling
p=inputParser;
addRequired(p, 'R', @(M) isnumeric(M));
mm = length(R);
addParameter(p, 'convergence_tol', 1e-5, @isscalar);
addParameter(p, 'max_steps', 1e4, @isscalar);
addParameter(p, 'P0', speye(mm));
addParameter(p, 'graphical_model', R~=0);
addParameter(p, 'printstuff', 1, @isscalar);
%addParameter(p, 'sample_size', 1, @isscalar);
addParameter(p, 'num_steps_check_convergence', 100, @isscalar);
addParameter(p, 'lambda', 0, @isscalar);

%addParameter(p, 'dampening', 0, @isscalar);


parse(p,R,varargin{:});
P = p.Results.P0; % Precision matrix
G = p.Results.graphical_model;% Graphical model
convergenceTol = p.Results.convergence_tol;
numStepsCheckConvergence = p.Results.num_steps_check_convergence;
printstuff = p.Results.printstuff;
lambda = p.Results.lambda;

%dampening = p.Results.dampening;
if ~all(G(P~=0))
    error('Initial precision matrix conflicts with graphical model')
end
[ii,jj] = find(G);

% objective function
obj = @(precision)objFn(R, precision);
objective_function_value = obj(P);

if nargout > 3 || printstuff == 2
    data.stepsize = zeros(p.Results.max_steps,1);
    data.dampening = zeros(p.Results.max_steps,1);
    data.obj = zeros(p.Results.max_steps,1);
    data.gradcorr = zeros(p.Results.max_steps,1);
    data.gradnorm = zeros(p.Results.max_steps,1);
end
    
% gradient descent
tic;
stepsize = 0.01;
dampening = 0;
gradient = G;
converged = false;
lastCheckedObj = inf;
for step = 1:p.Results.max_steps
    
    Pinv = sparseinv(P);
    
    % gradient of objective function
    oldGradient = gradient;
    gradient =  sparse(ii,jj,2 * (R(G) - Pinv(G)) + lambda * sign(P(G)),mm,mm);
    
    % dampening to eliminate oscillations
    gradient = dampening * oldGradient + (1 - dampening) * gradient;
    grad_corr = corr(gradient(G),oldGradient(G));
    dampening = max(0, min(0.99, dampening - grad_corr * 0.05));
    
    
    % line search to determine step size
    oldObj = objective_function_value;
    [P, stepsize, objective_function_value] = ...
        linesearch(P, objective_function_value, ...
        gradient, ...
        obj, stepsize);
    
    if nargout > 3 || printstuff == 2
        data.stepsize(step) = stepsize;
        data.dampening(step) = dampening;
        data.obj(step) = objective_function_value / mm;
        data.gradcorr(step) = grad_corr;
        data.gradnorm(step) = sqrt(sum(nonzeros(gradient.^2)) / mm) ;
    end
    
    % sometimes bad descent direction causes step size to drop
    % dramatically, and it takes a long time to reset
    if stepsize < 1e-6 / mm
        stepsize = 1;
        if printstuff == 2 && step > 10
            warning('Bad descent direction led to small stepsize on step %d', step)
        end
    end
    
    % convergence
    if mod(step,numStepsCheckConvergence)==0
        
        if printstuff == 2
            fprintf('step %d: F=%.4f, min ||grad||=%.4f\n', step, ...
                objective_function_value / mm, min(nonzeros(data.gradnorm)))
        elseif printstuff == 1
            fprintf('%d ',step)
        end
        
        if objective_function_value - lastCheckedObj > - mm  * numStepsCheckConvergence * convergenceTol
            converged = true;
            break;
        else
            lastCheckedObj = objective_function_value;
        end
        
    end
    
end

if printstuff == 2
    fprintf('\n time = %f\n', toc)
end

if ~converged && printstuff >= 1
    warning('Gradient descent failed to converge in %d steps', p.Results.max_steps);
end


% objective function
    function val = objFn(R, omega)
        [L, pp] = chol(omega);
        if pp == 0 % omega is positive definite
            val =  -0.5 * (sum(2*log(diag(L))) - dot(nonzeros(omega), R(omega~=0))) + lambda * sum(abs(nonzeros(omega)));
        else
            val = inf;
        end
    end
end

% line search to find optimal step size
function [thetaNew, step, newObj] = linesearch(initTheta, initObj, grad, objFn, step)

oldObj = initObj;
stepsize_factor = 1.5;

if ~isreal(oldObj); error('objective function value should be real at initial point for line search'); end

step = step * stepsize_factor;

newObj = objFn(initTheta - step * grad);
if 0
    while newObj < oldObj
        step = step * stepsize_factor;
        oldObj = newObj;
        newObj = objFn(initTheta - step * grad);
    end
    
    step = step / stepsize_factor;
    newObj = oldObj;
end

while newObj > initObj %- step * sum(nonzeros(grad).^2) / stepsize_factor
    step = step / stepsize_factor;
    newObj = objFn(initTheta - step * grad);
end

thetaNew = initTheta - step * grad;
end

