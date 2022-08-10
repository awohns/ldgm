function A = importGraph(filename, weighted, noSNPs)
%importGraph loads an adjacency matrix exported from NetworkX to 'filename'
%   weighted: whether matrix has edge weights
%   noSNPs: optional number of nodes different from maximum node ID

% Whether to import a weighted or binary adjacency matrix
if nargin < 2
    weighted = 1;
end

if weighted
    strformat = "%d %d {'weight': %f}";
else
    strformat = '%d %d {}';
end

% Load data
file=fopen(filename);
data=textscan(file,strformat);
if ~feof(file)
    error('failed to import file')
end

if nargin < 3
    noSNPs=double(max([data{1}; data{2}])+1);
end

% Weights 
if weighted
    weights = 1./(1+data{3});
else
    weights = true(length(data{1}),1);
end

% Generate sparse matrix
ii = double([data{1}+1; data{2}+1]);
jj = double([data{2}+1; data{1}+1]);
weights = [weights; weights];
incl = ii < jj;
A = sparse(ii(incl),jj(incl),weights(incl),noSNPs,noSNPs);
A = A + A';

end

