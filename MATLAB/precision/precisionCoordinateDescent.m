function [P] = precisionCoordinateDescent(A,R_sample,maxIters,lambda,lassoIters, P)
%precisionCoordinateDescent estimates the precision matrix P from the LDGM
%A and the sample correlation matrix R_sample, using a modified DP-GLASSO
%algorithm (Mazumder and Hastie, 2012).
% Input arguments:
% A: graphical model as a sparse symmetric matrix; precision matrix edges
% will be a subset of those of A
% R_sample: sample correlation matrix; can be sparse with entries
% corresponding to those of A
% maxIters (optional): number of coordinate descent outer loops (default
% 25)
% lambda (optional): L1 penalty term, usually something like 0.1 (default
% 0)
% lassoIters (optional): number of iterations for lasso regression inner
% loop (default 10)
% P (optional): initialization for P (default identity matrix)

nn = length(A);
% assert(all(~diag(A)))
if nargin < 6
    P = speye(nn);
end
if nargin < 5
    lassoIters = 10;
end
if nargin < 4
    lambda = 0;
end
if nargin < 3
    maxIters = 25;
end

% Lasso regression step unnecessary if there is no L1 penalty
if lambda == 0
    lassoIters = 0;
end

% initial minimum stepsize for lasso. Adaptively reduces this minimum if
% the lasso solver is failing to find good solutions
init_lasso_min_stepsize = 0.1;
lasso_increment = 2;
time0 = toc;
R_sample = R_sample - speye(nn);
fprintf('\n Beginning %d iterations by %d nodes of coordinate descent with lambda=%.2f \n', maxIters, nn, lambda)
for iter = 1:maxIters
    lasso_min_stepsize = init_lasso_min_stepsize;
    for ii = 1:nn
        neighbors = A(:,ii);
        
        % no-L1-penalty solution
        coef_init = precisionMultiply(P, R_sample(neighbors,ii), neighbors);            
        
        % proximal gradient descent Lasso solver
        if lassoIters > 0
            coef = solveLasso(P, R_sample(neighbors,ii), lambda, neighbors, lassoIters, lasso_min_stepsize, coef_init);
        else
            coef = coef_init;
        end
        
        % update precision matrix
        P_old = P;
        P(neighbors,ii) = - coef;
        P(ii,neighbors) = - coef;
        % different formula for the diagonal entry
        P(ii,ii) = 1 + R_sample(neighbors,ii)' * coef;
        
        % Check if solveLasso succeeded in finding a good solution
        if lambda > 0
            [~,p] = chol(P);
            counter = 0;
            
            while p > 0 % P is not positive definite
                counter = counter + 1;
                lasso_min_stepsize = lasso_min_stepsize / lasso_increment; % try a smaller min stepsize
                P = P_old;
                coef = solveLasso(P, R_sample(neighbors,ii), lambda, neighbors, lassoIters, lasso_min_stepsize, coef);
                if counter == 3 % give up on lasso regression for this coordinate
                    coef = precisionMultiply(P, R_sample(neighbors,ii), neighbors);
                    lasso_min_stepsize = lasso_min_stepsize * lasso_increment^3;
                    warning('Lasso solver failing to improve objective function on coordinate %d and outer loop %d; OK if this happens occasionally',ii,iter)
                elseif counter > 3
                    error('Fallback to no-L1-penalty coordinate descent solution failed to find p.s.d. precision matrix')
                end
                P(neighbors,ii) = - coef;
                P(ii,neighbors) = - coef;
                P(ii,ii) = 1 + R_sample(neighbors,ii)' * coef;
                [~,p] = chol(P);
                
            end
        end
    end
    time1 = toc;
    fprintf('Time for %d iteration(s): %.1f minutes\n', iter, (time1-time0)/60)
    
end
end

