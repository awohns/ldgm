function newFile = convertEdgelist(filename, snplist, one_indexed, newFile, convertIndices)
% Converts a .adjlist file to a .edgelist file

if isstr(snplist)
    snplist = readtable(snplist,'FileType','text');
end
if nargin < 3
    one_indexed = 0;
end
if nargin < 5
    convertIndices = true;
end

noDigits = 6;
assert(strcmp(filename(end-7:end),'.adjlist'))
strformat = "%d %d {'weight': %f}";
file=fopen(filename);
data=textscan(file,strformat);

% SNP indices of entries of sparse matrix
ii = double(data{1}+1-one_indexed);
jj = double(data{2}+1-one_indexed);
weights = data{3};
assert(max(max(ii(:)),max(jj(:))) <= height(snplist));

% Convert from SNP indices to row/col indices
if convertIndices
    ii = snplist.index(ii);
    jj = snplist.index(jj);
else
    % make zero-indexed
    ii = ii - 1;
    jj = jj - 1;
end

% Upper-triangular
incl = ii <= jj;

if nargin < 4
    newFile = [filename(1:end-8),'.edgelist'];
end
writematrix([ii(incl), jj(incl), round(weights(incl), noDigits)],...
    newFile, 'filetype','text');

end