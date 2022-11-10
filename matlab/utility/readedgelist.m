function A = readedgelist(filepath, noVertices, replace_zeros)
%readedgelist inputs an edge list filename and returns a sparse matrix with
%edges listed in the file. Optionally, specify noVertices and the matrix
%will have the appropriate size (otherwise, it might be smaller than
%desired if there are no edges for the last vertex or vertices).
%Optionally, specify replace_zeros and edges with zero entries will be
%replaced with this number.

A = readmatrix(filepath, 'FileType', 'text');
assert(size(A,2) == 3, 'edgelist file should have three columns')

% Replace zeros with replace_zeros
if nargin >= 3
    A(A(:,3)==0,3) = replace_zeros;
end

% Convert to triangular matrix if necessary
if any(A(:,1) > A(:,2))
    offdiagonal = A(:,1) ~= A(:,2);
    A = [A; A(offdiagonal,[2 1 3])];
    A = A(A(:,2) >= A(:,1), :);
end

if nargin < 2
    noVertices = max(A(:,2))+1;
end

% Get rid of redundant rows if present
[~,ui] = unique(A(:,1)+noVertices*A(:,2),'stable');
A = A(ui,:);

% Convert to symmetric matrix
offdiagonal = A(:,1) ~= A(:,2);
A = [A; A(offdiagonal,[2 1 3])];



% Create one-indexed sparse matrix
A = sparse(A(:,1)+1, A(:,2)+1, A(:,3),noVertices,noVertices);

assert(issymmetric(A))
end

