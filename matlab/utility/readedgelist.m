function A = readedgelist(filepath)
%readedgelist inputs an edge list filename and returns a sparse matrix with
%edges listed in the file

A = readmatrix(filepath, 'FileType', 'text');
assert(size(A,2) == 3, 'edgelist file should have three columns')
A = sparse(A(:,1)+1, A(:,2)+1, A(:,3));

end

