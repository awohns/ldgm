function convertEdgelist(filename)
% Converts a .adjlist file to a .edgelist file
noDigits = 6;
assert(strcmp(filename(end-7:end),'.adjlist'))
newFilename = [filename(1:end-8),'.edgelist'];
strformat = "%d %d {'weight': %f}";
file=fopen(filename);
data=textscan(file,strformat);
% Generate sparse matrix
ii = double([data{1}+1; data{2}+1]);
jj = double([data{2}+1; data{1}+1]);
weights = [weights; weights];
incl = ii < jj;

writematrix([jj(incl), ii(incl), round(edges(incl), noDigits)],...
    newFilename, 'filetype','text');

end