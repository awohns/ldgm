function x = precisionDivide(P, y, whichIndices)
%precisionDivide solves (P/P00)*x = y where P = [P11, P10; P01 P00],
% P11 = P(whichIndices,whichIndices), and P/P00 is the Schur complement
% Input arguments:
% P: precision matrix or cell array of precision matrices
% y: vector, matrix or cell array of vectors/matrices with same dimension as P
% whichIndices: indices or cell array of indices, or boolean vector/cell
% array of boolean vectors, indicating which entries of P correspond to the
% entries of y
% Output arguments:
% x == (P/P00) \ y

if iscell(P)
    assert(iscell(y) & iscell(whichIndices))
    x = cellfun(@precisionDivide, P, y, whichIndices, 'UniformOutput', false);
else
    if isa(whichIndices,'logical')
        whichIndices = find(whichIndices);
    end
    
    % handle zero diagonal elements in P
    incl = diag(P) ~= 0;
    assert(all(incl(whichIndices)),'Precision matrix should have nonzero diagonal entries for all SNPs in y')
    
    % yp is y augmented with zeros
    yp = zeros(size(P,1),size(y,2));
    yp(whichIndices,:) = y;
    
    % xp is x augmented with entries that can be ignored
    xp = zeros(size(yp));
    xp(incl,:) = P(incl,incl) \ yp(incl,:);
    x = xp(whichIndices,:);
end


end

