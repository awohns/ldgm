function writeedgelist(filename,A)
%writeedgelist writes an edgelist file from a symmetric matrix
% Input arguments:
% filename: name of file to be written, normally ending with .edgelist
% A: symmetric matrix, usually sparse, to be recorded

noDigits = 6;
if ~issymmetric(A)
    warning('A should be a symmetric matrix')
end
[ii,jj,entries] = find(A);

% only record upper triangle
incl = ii <= jj;

writematrix([ii(incl)-1, jj(incl)-1, round(entries(incl), noDigits)],...
    filename, 'filetype','text');
end