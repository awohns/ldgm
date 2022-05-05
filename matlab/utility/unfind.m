
function lgc = unfind(idx, N)
    %Go from indicies into logical (for vectors only)
    assert(N >= max(idx));
    lgc = false(N,1);
    lgc(idx) = true;
end