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

smallNumber = 1e-12;

%dampening = p.Results.dampening;
if ~all(G(P~=0))
    error('Initial precision matrix conflicts with graphical model')
end
[ii,jj] = find(G);
GnoSelfEdges = sparse(ii,jj,ii~=jj,mm,mm);

% objective function
obj = @(precision)objFn(R,G,precision);
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
for rep = 1:p.Results.max_steps
    
    Pinv = sparseinv(P);
    
    % gradient of objective function
    oldGradient = gradient;
    penalty = lambda .* sign(P(G));
    penalty(ii==jj) = 0;
    gradient =  sparse(ii,jj,(R(G) - Pinv(G)) + penalty,mm,mm);
    
    % dampening to eliminate oscillations
%     gradient = dampening * oldGradient + (1 - dampening) * gradient;
     grad_corr = corr(gradient(G),oldGradient(G));
%     dampening = max(0, min(0.99, dampening - grad_corr * 0.05));
%     grad_corr=0;
    
    % backtracking line search to determine step size
    [P, stepsize, objective_function_value] = ...
        linesearch(P, objective_function_value, ...
        gradient, ...
        obj, stepsize);

    
    if nargout > 3 || printstuff == 2
        data.stepsize(rep) = stepsize;
        data.dampening(rep) = dampening;
        data.obj(rep) = objective_function_value / mm;
        data.gradcorr(rep) = grad_corr;
        data.gradnorm(rep) = sqrt(sum(nonzeros(gradient.^2)) / mm) ;
    end
    
    % sometimes bad descent direction causes step size to drop
    % dramatically, and it takes a long time to reset
    if stepsize < 1e-12 / mm
        stepsize = 0.01;
        if printstuff == 2 && rep > 10
            warning('Bad descent direction led to small stepsize on step %d', rep)
        end
    end
    
    % convergence
    if mod(rep,numStepsCheckConvergence)==0
        
        if printstuff == 2
            fprintf('step %d: F=%.4f, min ||grad||=%.4f\n', rep, ...
                full(objective_function_value / mm), full(min(nonzeros(data.gradnorm))))
        elseif printstuff == 1
            fprintf('%d ',rep)
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
elseif printstuff == 1
    fprintf('\n')
end

if ~converged && printstuff >= 1
    warning('Gradient descent failed to converge in %d steps', p.Results.max_steps);
end


% objective function
    function val = objFn(R,G, omega)
        [L, pp] = chol(omega);
        if pp == 0 % omega is positive definite
            val =  -1 * (sum(2*log(diag(L))) - sum(omega(G).* R(G))) + lambda * sum(abs(omega(GnoSelfEdges)));
        else
            val = inf;
        end
    end

% Step function with penalty
    function thetaNew = step(thetaInit, grad, stepsize)
        thetaInit = full(thetaInit(G));
        grad = full(grad(G));
        
        % when starting at 0, shrink gradient toward 0 by lambda
        grad(thetaInit == smallNumber) = sign(grad(thetaInit == smallNumber))...
            .* max(0, abs(grad(thetaInit == smallNumber)) - lambda);
        
        % update theta
        thetaNew = thetaInit - stepsize * grad;
        
        % if theta switches sign, set to smallNumber (not 0 for numerical
        % reasons)
        thetaNew(sign(thetaNew).*sign(thetaInit)==-1) = smallNumber;
        
        thetaNew = sparse(ii,jj,thetaNew,mm,mm);
        
    end

% line search to find optimal step size
    function [thetaNew, stepsize, newObj] = linesearch(thetaInit, initObj, grad, objFn, stepsize)
        
        oldObjVal = initObj;
        stepsize_factor = 2;
        
        if ~isreal(oldObjVal); error('objective function value should be real at initial point for line search'); end
        
        stepsize = stepsize * sqrt(stepsize_factor);
        
        thetaNew = step(thetaInit, grad, stepsize);
        newObj = objFn(thetaNew);
        if 0
            while newObj < oldObjVal
                stepsize = stepsize * stepsize_factor;
                oldObjVal = newObj;
                newObj = objFn(step(thetaInit, grad, stepsize));
            end
            
            stepsize = stepsize / stepsize_factor;
            newObj = oldObjVal;
        end
        
        while newObj > initObj %- step * sum(nonzeros(grad).^2) / stepsize_factor
            stepsize = stepsize / stepsize_factor;
            thetaNew = step(thetaInit, grad, stepsize);
            newObj = objFn(thetaNew);
        end
        
    end

end


