function A = readedgelist(filepath, noVertices)
%readedgelist inputs an edge list filename and returns a sparse matrix with
%edges listed in the file. Optionally, specify noVertices and the matrix
%will have the appropriate size (otherwise, it might be smaller than
%desired if there are no edges for the last vertex or vertices)

A = readmatrix(filepath, 'FileType', 'text');
assert(size(A,2) == 3, 'edgelist file should have three columns')

% Convert to triangular matrix if necessary
if any(A(:,1) > A(:,2))
    offdiagonal = A(:,1) ~= A(:,2);
    A = [A; A(offdiagonal,[2 1 3])];
    A = A(A(:,2) >= A(:,1), :);
end

% Convert to symmetric matrix
offdiagonal = A(:,1) ~= A(:,2);
A = [A; A(offdiagonal,[2 1 3])];

if nargin < 2
    noVertices = max(max(A(:,1)))+1;
end

% Create one-indexed sparse matrix
A = sparse(A(:,1)+1, A(:,2)+1, A(:,3),noVertices,noVertices);

assert(issymmetric(A))
end

