
function lgc = unfind(idx, N)
%Go from indicies into logical (for vectors only)
if ~isempty(idx)
    assert(N >= max(idx));
    assert(issorted(idx),'Only use unfind when indices are sorted')
end
lgc = false(N,1);
lgc(idx) = true;

end