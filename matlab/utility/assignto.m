
function y = assignto(x, lgc, n)
    if ~islogical(lgc)
        lgc = unfind(lgc,n);
    end
    % Assign values in x to the nonzeros of lgc
    assert(sum(lgc) == length(x))
    y = zeros(size(lgc));
    y(lgc) = x;
end