function [A] = reduce_weighted_graph(A,missing)
%reduce_weighted_graph removes rows/columns "missing" from A and patches
%paths between remaining vertices

mm=length(A);
for jj=1:length(missing)
    ii = missing(jj);
    Ni=find(A(:,ii));
    Ni=setdiff(Ni,ii);
%     A(Ni,Ni)=max( A(Ni,Ni),...
%         1./( (1./A(Ni,ii) - 1) + (1./A(ii,Ni) - 1) + 1) );
    A(Ni,Ni) = min(
end

notmissing = setdiff(1:mm,missing);
A = A(notmissing,notmissing);

end

