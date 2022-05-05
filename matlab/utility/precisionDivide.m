function x = precisionDivide(P, y, notMissing)
%precisionDivide solves (P/P00)*x = y where P = [P11, P10; P01 P00] and
% P11 = P(notMissing,notMissing)
% Input arguments:
% P: precision matrix or cell array of precision matrices
% y: vector or cell array of vectors with same dimension as P
% notMissing: indices or cell array of indices, or boolean vector/cell
% array of boolean vectors, indicating which entries of P correspond to the
% entries of y (i.e. P11 = P(notMissing,notMissing))
% Output arguments:
% x: solves (P/P00)*x == y

if iscell(P)
    assert(iscell(y) & iscell(notMissing))
    x = cellfun(@precisionDivide, P, y, notMissing, 'UniformOutput', false);
else
    % yp is y augmented with zeros
    if isa(notMissing,'logical')
        notMissing = find(notMissing);
    end
    yp = sparse(notMissing,ones(length(notMissing),1),y,length(P),1);
    
    % xp is x augmented with entries that can be ignored
    xp = P \ yp;
    x = xp(notMissing);
end


end

