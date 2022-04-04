function [A] = reduce_weighted_graph(A,missing)
%reduce_weighted_graph removes rows/columns "missing" from A and patches
%paths between remaining vertices

mm=length(A);
for j=1:length(missing)
    i = missing(j);
    Ni=find(A(:,i));
    Ni=setdiff(Ni,i);
    A(Ni,Ni)=max( A(Ni,Ni),...
        1./( (1./A(Ni,i) - 1) + (1./A(i,Ni) - 1) + 1) );
end

notmissing = setdiff(1:mm,missing);
A = A(notmissing,notmissing);

end

