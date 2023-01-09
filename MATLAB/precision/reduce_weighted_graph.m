function [A] = reduce_weighted_graph(A, missing, is_weighted)
%reduce_weighted_graph removes rows/columns "missing" from A and patches
% paths between remaining vertices. 
% A: weighted graph, with '1' corresponding to a strong edge.
% missing: which SNPs are missing from A (logical or index)
% is_weighted: true -> edge weights should be between 0 and 1 (default); 0
% -> edge weights represent distances and should be between 0 and inf

if islogical(missing)
    missing = find(missing);
end
if nargin < 3
    is_weighted = true;
end
if is_weighted
    assert(max(A(:)) <= 1, 'A should have weights between 0 and 1')
end
assert(min(A(:)) >= 0, 'A should have nonnegative weights/distances')
assert(length(A) >= max(missing))

% convert distances to weights
if ~is_weighted
    A = spfun(@(x)1./(x+1), A);
end

mm=length(A);
for jj=1:length(missing)
    ii = missing(jj);
    Ni=find(A(:,ii));
    Ni=setdiff(Ni,ii);

    % 1/weight - 1 converts weight -> distance; 1/(distance + 1) converts
    % distance back to a weight
    A(Ni,Ni)=max( A(Ni,Ni),...
        1./( (1./A(Ni,ii) - 1) + (1./A(ii,Ni) - 1) + 1) );

end

% convert weights to distances
if ~is_weighted
    A = spfun(@(x)1./x - 1, A);
end

notmissing = setdiff(1:mm,missing);
A = A(notmissing,notmissing);

end

