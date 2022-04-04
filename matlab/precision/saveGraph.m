function [outputArg1,outputArg2] = saveGraph(filename,A,v)
%saveGraph saves a networkX-formatted weighted adjacency matrix, A, with
%vertex labels v, with file directory + name filename

assert(isstr(filename))
assert(issparse(A))
assert(isvector(v))
if isrow(v); v = v'; end
assert(isnumeric(v))
strformat = "%d %d {'weight': %f}\n";

[i,j,e] = find(A);

file = fopen(filename,'w');
fprintf(file,strformat,[v(i)';v(j)';e']);
fclose(file);
end

