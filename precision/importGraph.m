function A = importGraph(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
file=fopen(filename);
data=textscan(file,'%d %d {}');
if ~feof(file)
    error('failed to import file')
end
nn=max([data{1}; data{2}])+1;
A = sparse(data{1}+1,data{2}+1,true(length(data{1}),1),nn,nn);
A = A + A' ~= 0;

end

