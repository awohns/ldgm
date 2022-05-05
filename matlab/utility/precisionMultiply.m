function x = precisionMultiply(P, y, notMissing)
%precisionMultiply computes x = (P/P00)y where P = [P11, P10; P01 P00] and
%P11 = P(notMissing,notMissing)
% Input arguments:
% P: precision matrix or cell array of precision matrices
% y: vector or cell array of vectors with same dimension as P
% notMissing: indices or cell array of indices, or boolean vector/cell
% array of boolean vectors, indicating which entries of P correspond to the
% entries of y (i.e. P11 = P(notMissing,notMissing))
% Output arguments:
% x: solves x == (P/P00)*y

if iscell(P)
    assert(iscell(y) & iscell(notMissing))
    x = cellfun(@precisionMultiply, P, y, notMissing, 'UniformOutput', false);
else
    x = P(notMissing,notMissing) * y - P(notMissing,~notMissing) * ...
        (P(~notMissing,~notMissing) \ (P(~notMissing,notMissing) * y));
end
end

