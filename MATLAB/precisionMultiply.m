function x = precisionMultiply(P, y, whichIndices, ldl_factor)
%precisionMultiply computes x = (P/P00)y where P = [P11, P10; P01 P00],
% P11 = P(whichIndices,whichIndices), and P/P00 is the Schur complement
%
% Input arguments:
% P: precision matrix or cell array of precision matrices
% y: vector or cell array of vectors with same dimension as P
% whichIndices: indices or cell array of indices, or boolean vector/cell
% array of boolean vectors, indicating which entries of P correspond to the
% entries of y
% ldl_factor: Optionally, precompute the LDL decomposition of P using
% ldlchol() and use it to speed up the solving step in the presence of
% many missing values. This will only help when the number of non-mmissing
% values is small (e.g., 10), and it will actually be slower when the number
% of non-missing values is large (e.g., 100). If this option is used,
% you must subset P to the rows/columns with nonzeros on the diagonal, then
% call [LD, ~, perm] = ldlchol(P), then permute both P and
% whichIndices, then call precisionMultiply(P(perm,perm), y,
% perm(whichIndices), LD). If you compute LD = ldlchol(P) without first 
%
% Output arguments:
% x == (P/P00)*y

if iscell(P)
    assert(iscell(y) & iscell(whichIndices))
    if ~exist('ldl_factor')
        ldl_factor = cell(size(P));
    end
    x = cellfun(@precisionMultiply, P, y, whichIndices, ldl_factor, 'UniformOutput', false);
else
    if nargin < 4
        ldl_factor = [];
    end

    % handle SNPs not in the LDGM
    incl = diag(P)~=0;
    assert(all(incl(whichIndices)),'Precision matrix should have nonzero diagonal entries for all non-missing SNPs')
    
    
    % If LDL factor specified, and y is sufficiently short,     
    % speed up the solving step by taking subset of inverse of P
    if isempty(ldl_factor) || length(y) > 100
        otherIndices = incl; 
        otherIndices(whichIndices) = false;

        x = P(otherIndices,whichIndices) * y;
        x = P(otherIndices,otherIndices) \ x;
        x = P(whichIndices,whichIndices) * y - P(whichIndices,otherIndices) * x;

    else
        assert(all(incl), 'If supplying LDL factor, P should be positive definite with no zeros on its diagonal')
        
        % columns of inv(P)
        I = speye(sum(incl));
        PInvCols = ldlsolve(ldl_factor, I(:,whichIndices));
        x = PInvCols(whichIndices,:) \ y;

    end


end
end

