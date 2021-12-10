function A = importGraph(filename, weighted)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    weighted = 0;
end
if weighted
    strformat = "%d %d {'weight': %f}";
else
    strformat = '%d %d {}';
end
file=fopen(filename);
data=textscan(file,strformat);
if ~feof(file)
    error('failed to import file')
end
nn=double(max([data{1}; data{2}])+1);
if weighted
    weights = 1./(1+data{3});
else
    weights = true(length(data{1}),1);
end
ii = double([data{1}+1; data{2}+1]);
jj = double([data{2}+1; data{1}+1]);
weights = [weights; weights];
incl = ii < jj;

A = sparse(ii(incl),jj(incl),weights(incl),nn,nn);
A = A + A';

end

