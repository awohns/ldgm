function A = importGraph(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
file=fopen(filename);
data=textscan(file,'%d %d {}');
if ~feof(file)
    error('failed to import file')
end
nn=max([data{1}; data{2}])+1;
<<<<<<< Updated upstream
A = sparse(data{1}+1,data{2}+1,true(length(data{1}),1),nn,nn);
A = A + A' ~= 0;
=======
if weighted
    weights = 1./(1+data{3});
else
    weights = true(length(data{1}),1);
end
ii = [data{1}+1; data{2}+1];
jj = [data{2}+1; data{1}+1];
weights = [weights; weights];
incl = ii < jj;

A = sparse(ii(incl),jj(incl),weights(incl),nn,nn);
A = A + A';
>>>>>>> Stashed changes

end

