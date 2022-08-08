function [A, p, s] = precisionReorderCholesky(P, idcs)
%precisionReorderCholesky computes the Cholesky factor of P(s,s) where s 
% is a permutation of the indices of P with s(end-k+1:end) == idcs.
% This is useful when computing things the form
% x = invP(idcs,idcs)*y, invP == inv(P)
% because it allows the use of the precisionCholeskyDivide function:
% instead of 
% invP = inv(P); x = invP(idcs,idcs)*y;
% or the faster
% x = precisionDivide(P,y,idcs);
% compute:
% A = precisionReorderCholesky(P, idcs);
% x = precisionCholeskyDivide(A, y);
% This is much faster if you need to repeat this for many vectors y.
% 
% With three output arguments, precisionReorderCholesky returns the entire
% matrix A, the flag p (which should be 0), and the permutation that it 
% used, s; with 1 output argument, it only returns the last part of A, 
% correponding to idcs.

if iscell(P)
    assert(iscell(idcs))
    A = cellfun(@precisionReorderCholesky, P, idcs, 'UniformOutput', false);
else
    if isrow(idcs); idcs = idcs'; end
    if islogical(idcs); idcs = find(idcs); end
    
    s = [setdiff(1:length(P),idcs)'; idcs];
    
    [A, p] = chol(P(s,s));
    
    if nargout == 1
        k = length(idcs);
        A = A(end-k+1:end, end-k+1:end);
    end
        
end


end

