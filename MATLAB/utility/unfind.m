
function [lgc, permutation, inverse_permutation] = unfind(idx, N)
%Go from indicies into logical. Supply two output arguments to handle
%non-sorted indices; afterwards, y = x(idx) is the same as
% y = x(lgc); y = y(permutation);

if nargin < 2
    N = max(idx);
end
assert(N >= max(idx));
if nargout < 2
    assert(issorted(idx),'Only use unfind when indices are sorted')
else
    [~, inverse_permutation] = sort(idx);
    permutation(inverse_permutation) = 1:length(inverse_permutation);
end
lgc = false(N,1);
lgc(idx) = true;

end