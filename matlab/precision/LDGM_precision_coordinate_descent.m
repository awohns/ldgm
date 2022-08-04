function [P] = LDGM_precision_coordinate_descent(A,R_sample,maxIters,lambda,lassoIters, P)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
nn = length(A);
% assert(all(~diag(A)))
if nargin < 7
    P = speye(nn);
end
if nargin < 4
    lassoIters = 2;
end
if nargin < 4
    lambda = 0;
end

% initial minimum stepsize for lasso. Adaptively reduces this minimum if
% the lasso solver is failing to find good solutions
init_lasso_min_stepsize = 0.1;
lasso_increment = 2;
tic;
R_sample = R_sample - speye(nn);
for iter = 1:maxIters
    lasso_min_stepsize = init_lasso_min_stepsize;
    for ii = 1:nn%rp
        if mod(ii,100)==0
            disp(ii);toc
        end
        neighbors = A(:,ii);
        
        
        if lassoIters > 0
            if 1
                coef_init = precisionMultiply(P, R_sample(neighbors,ii), neighbors);
            else
                coef_init = -P(neighbors,ii);
                coef_init(find(neighbors)==ii) = coef_init(find(neighbors)==ii) + 1;
            end
            
            coef = solveLasso(P, R_sample(neighbors,ii), lambda, neighbors, lassoIters, lasso_min_stepsize, coef_init);
            
        else
            coef = precisionMultiply(P, R_sample(neighbors,ii), neighbors);
        end
        P_old = P;
        P(neighbors,ii) = - coef;
        P(ii,neighbors) = - coef;
        P(ii,ii) = 1 + R_sample(neighbors,ii)' * coef;
        if lambda > 0
            [~,p] = chol(P);
            counter = 0;
            while p > 0
                counter = counter + 1;
                lasso_min_stepsize = lasso_min_stepsize / lasso_increment;
                P = P_old;
                coef = solveLasso(P, R_sample(neighbors,ii), lambda, neighbors, lassoIters, lasso_min_stepsize, coef);
                if counter == 3
                    coef = solveLasso(P, R_sample(neighbors,ii), 0, neighbors, 0, lasso_min_stepsize);
                    lasso_min_stepsize = lasso_min_stepsize * lasso_increment^3;
                    warning('Lasso solver failing to improve objective function')
                elseif counter > 3
                    error('Failing to find p.s.d. precision matrix')
                end
                P(neighbors,ii) = - coef;
                P(ii,neighbors) = - coef;
                P(ii,ii) = 1 + R_sample(neighbors,ii)' * coef;
                [~,p] = chol(P);
                
            end
        end
    end
    %     R = sparseinv(P);
    %     d = diag(sqrt(diag(R)));
    %     P = d * P * d;
    %iter_timer(iter) = toc
end
end

