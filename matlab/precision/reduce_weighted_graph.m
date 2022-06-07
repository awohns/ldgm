function [A] = reduce_weighted_graph(A, missing)
%reduce_weighted_graph removes rows/columns "missing" from A and patches
%paths between remaining vertices. 
% A: weighted graph, with '1' corresponding to a strong edge.
% missing: which SNPs are missing from A (logical or index)

if islogical(missing)
    missing = find(missing);
end
assert(max(A(:)) <= 1, 'A should have weights between 0 and 1')
assert(min(A(:)) >= 0, 'A should have weights between 0 and 1')
assert(length(A) >= max(missing))
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

notmissing = setdiff(1:mm,missing);
A = A(notmissing,notmissing);

end

