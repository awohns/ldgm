function x = precisionDivide(P, y, whichIndices, use_ldlchol)
%precisionDivide solves (P/P00)*x = y where P = [P11, P10; P01 P00],
% P11 = P(whichIndices,whichIndices), and P/P00 is the Schur complement
% Input arguments:
% P: precision matrix or cell array of precision matrices
% y: vector, matrix or cell array of vectors/matrices with same dimension as P
% whichIndices: indices or cell array of indices, or boolean vector/cell
% array of boolean vectors, indicating which entries of P correspond to the
% entries of y
% use_ldlchol: optionally, specify use_ldlchol == true to use a precomputed
% LDL cholesky decomposition in place of the precision matrix. This requires
% installing SuiteSparse. Then, run:
%   [A, ~, perm] = ldlchol(P); 
%   x(perm) = precisionDivide(A, y(perm), whichIndices(perm), true);
% where perm is a fill-reducing permutation and A is the LD factor of P.
% 
% Output arguments:
% x == (P/P00) \ y

if nargin < 4
    use_ldlchol = false;
end

if iscell(P)
    assert(iscell(y) & iscell(whichIndices))
    use_ldlchol = num2cell(repmat(use_ldlchol,size(P)));
    x = cellfun(@precisionDivide, P, y, whichIndices, use_ldlchol, 'UniformOutput', false);
else

    if isa(whichIndices,'logical')
        whichIndices = find(whichIndices);
    end

    % handle zero diagonal elements in P
    incl = diag(P) ~= 0;
    if use_ldlchol
        assert(all(incl), 'Cholesky factor should have all nonzero diagonal entries; please compute it using ldlchol')
        assert(istril(P), 'Cholesky factor should be lower triangular; please compute it using ldlchol')
    else
        assert(all(incl(whichIndices)),...
            'Precision matrix should have nonzero diagonal entries for all SNPs in y')
    end

    % yp is y augmented with zeros
    yp = zeros(size(P,1),size(y,2));
    yp(whichIndices,:) = y;

    % xp is x augmented with entries that can be ignored
    xp = zeros(size(yp));
    if use_ldlchol
        xp = ldlsolve(P, yp);
    else
        xp(incl,:) = P(incl,incl) \ yp(incl,:);
    end
    x = xp(whichIndices,:);
end


end

