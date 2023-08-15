
function y = assignto(x, lgc, n)
% Assign values in x to the nonzeros of lgc. lgc can be specified
% instead as a list of indices, in which case third input argument sets
% the number of rows of output y.
if isempty(lgc)
    assert(isempty(x))
    if exist('n')
        y = zeros(n,1);
    else
        y = [];
    end
else
    if ~islogical(lgc)
        lgc = unfind(lgc,n);
    end
    assert(sum(lgc) == length(x))
    y = zeros(size(lgc,1), size(x,2));
    y(lgc,:) = x;
end
end