function x = precisionMultiply(P, y, whichIndices)
%precisionMultiply computes x = (P/P00)y where P = [P11, P10; P01 P00],
% P11 = P(whichIndices,whichIndices), and P/P00 is the Schur complement
% 
% Input arguments:
% P: precision matrix or cell array of precision matrices
% y: vector or cell array of vectors with same dimension as P
% whichIndices: indices or cell array of indices, or boolean vector/cell
% array of boolean vectors, indicating which entries of P correspond to the
% entries of y
% Output arguments:
% x == (P/P00)*y

if iscell(P)
    assert(iscell(y) & iscell(whichIndices))
    x = cellfun(@precisionMultiply, P, y, whichIndices, 'UniformOutput', false);
else
    
    % handle SNPs not in the LDGM
    incl = diag(P)~=0;
    assert(all(incl(whichIndices)),'Precision matrix should have nonzero diagonal entries for all non-missing SNPs')
    otherIndices = incl; otherIndices(whichIndices) = false;
    
    x = P(whichIndices,whichIndices) * y - P(whichIndices,otherIndices) * ...
        (P(otherIndices,otherIndices) \ (P(otherIndices,whichIndices) * y));
end
end

