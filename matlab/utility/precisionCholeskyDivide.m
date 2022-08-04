function x = precisionCholeskyDivide(A, y)
%precisionCholeskyDivide solves (P/P00)*x = y where 
% P = [P11, P10; P01 P00] = A'*A and
% P11 = P(end-k+1:end,end-k+1:end), k = length(y).
% A must be upper-triangular. Optionally, it can be subsetted to the
% nonmissing SNPs already.

if iscell(A)
    assert(iscell(y))
    x = cellfun(@precisionCholeskyDivide, A, y, 'UniformOutput', false);
else
    k = length(y);
    x = A(end-k+1:end,end-k+1:end) \ (A(end-k+1:end,end-k+1:end)' \ y);
end


end

