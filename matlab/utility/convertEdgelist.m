function newFilename = convertEdgelist(filename, snplist, one_indexed)
% Converts a .adjlist file to a .edgelist file

if isstr(snplist)
    snplist = readtable(snplist,'FileType','text');
end
if nargin < 3
    one_indexed = 0;
end

noDigits = 6;
assert(strcmp(filename(end-7:end),'.adjlist'))
newFilename = [filename(1:end-8),'.edgelist'];
strformat = "%d %d {'weight': %f}";
file=fopen(filename);
data=textscan(file,strformat);

% SNP indices of entries of sparse matrix
ii = double([data{1}+1-one_indexed; data{2}+1-one_indexed]);
jj = double([data{2}+1-one_indexed; data{1}+1-one_indexed]);
weights = [data{3}; data{3}];
assert(max(ii(:)) <= height(snplist));

% Convert from SNP indices to row/col indices
ii = snplist.index(ii);
jj = snplist.index(jj);

% Upper-triangular
incl = ii <= jj;

writematrix([ii(incl), jj(incl), round(weights(incl), noDigits)],...
    newFilename, 'filetype','text');

end