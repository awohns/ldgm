function [Omega, X] = makeSparsePrecision(nodes,density,nn)
% Generates random sparse precision matrix Omega with specified number of
% nodes, approximate density. Omega is equal to XX' where X has size nodes
% by nn. nn/nodes controls spectrum of matrix (make this larger to make
% Omega more well behaved)

if nargin < 3
    nn = nodes * 2;
end

% proportion of nonzeros in X s.t. XX' has about the right number of edges.
% proportion of edges wil be around 1-(1-Xdensity^2)^nn ~=
% 1-exp(-nn*Xdensity^2)
Xdensity = sqrt(-log(1-density)/nn);

Omega = 0;
while det(Omega)<=0
    % nonzero entries of X
    ind = rand(nodes, nn) < Xdensity;
    
    X = zeros(nodes, nn);
    X(ind) = randn(nnz(ind),1);
    Omega = X*X';
    
    D=diag(Omega);
    
    % get rid of degeneracies
    Omega = Omega + diag(D==0);
    D(D==0)=1;
    
    % normalize
    D=diag(sqrt(1./D));
    Omega = D*Omega*D;
end
Omega = sparse(Omega);

end

