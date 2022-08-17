function u = solveLasso(P,y,lambda,neighbors,noIter,min_stepsize,u)
%solveLasso uses proximal gradient descent (Combettes and Wajs, 2005) to
%solve a lasso regression problem:
% P_inv(neighbors,neighbors) * u - y + lambda*sign(u) == 0
% where P_inv == inv(P). This is the gradient of the lasso regression
% problem:
% argmin 0.5*sum((P_inv(neighbors,neighbors) * beta - y).^2) + lambda * sum(abs(beta))
% 
% Input arguments:
% P, y, lambda, idx: problem definition
% noIter: number of iterations
% u: initial guess for the coefficients

if nargin < 5
    noIter = 10;
end
if nargin < 6
    min_stepsize = 0.1;
end
if nargin < 7
    u = precisionMultiply(P, y, neighbors);
end

% Cholesky factor of P_inv(neighbors,neighbors)
A = precisionReorderCholesky(P, neighbors);

obj = @(u,x,y)0.5*sum(u.*x)+y'*u + lambda*sum(abs(u));
obj_val_new = inf;
gamma = 1;
counter = 0;
while counter < noIter
    counter = counter + 1;
    obj_val_old = obj_val_new;
    
    % equivalent to x = precisionDivide(P, u, neighbors), but much faster
    % with small number of neighbors
    x = precisionCholeskyDivide(A, u);
    
    % objective function value
    obj_val_new = obj(u, x, y);
    
    % if failing to improve objective, usually indicates stepsize too
    % large; reduce stepsize by half and try again
    if obj_val_new > obj_val_old
        counter = counter - 1;
        gamma = gamma / 2;
        u = u_old;
        x = precisionCholeskyDivide(A, u);
        obj_val_new = obj_val_old;
        if gamma < min_stepsize
            break;
        end
    end
    
    u_old = u;
    
    % update u
    u = prox(u - gamma  * (x - y) , gamma * lambda);
    
end
    
% proximity function
    function x = prox(x,gamma)
        x = sign(x) .* max(0,abs(x) - gamma);
    end
end

