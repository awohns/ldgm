function u = solveLasso(P,y,lambda,neighbors,noIter,min_stepsize,init)
%solveLasso solves:
% X(idx,idx) * beta - y + lambda*sign(beta) == 0
% where Xinv == inv(X). This is the gradient of the lasso regression
% problem:
% argmin 0.5*sum((X(idx,idx) * beta - y).^2) + lambda * sum(abs(beta))
% Input arguments:
% Xinv, y, lambda, idx: problem definition
% noIter: number of iterations
% init: initial guess

if nargin < 5
    noIter = 3;
end
if nargin < 6
    min_stepsize = 0.01;
end
if nargin < 7
    u =  precisionMultiply(P, y, neighbors);% zeros(length(idx),1);
else
    u = init;
end

obj = @(u,x,y)0.5*sum(u.*x)+y'*u + lambda*sum(abs(u));
obj_val_new = inf;
gamma = 1;
counter = 0;
while counter < noIter
    counter = counter + 1;
    obj_val_old = obj_val_new;
    x = precisionDivide(P, u, neighbors);
    obj_val_new = obj(u, x, y);
    if obj_val_new > obj_val_old
        counter = counter - 1;
        gamma = gamma / 2;
        u = u_old;
        x = precisionDivide(P, u, neighbors);
        obj_val_new = obj_val_old;
        if gamma < min_stepsize
            break;
%             warning('solveLasso getting stuck')
        end
    else
%         fprintf('%.3f ',gamma)
    end
    u_old = u;
    u = prox(u - gamma  * (x - y) , gamma * lambda);
    
    %u = u_next .* (sign(u).*sign(u_next) >= 0);
end
% fprintf('\n')

    function x = prox(x,gamma)
        x = sign(x) .* max(0,abs(x) - gamma);
    end
end

