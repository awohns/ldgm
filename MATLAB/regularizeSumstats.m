function z_regularized = regularizeSumstats(P,z,whichIndices,lambdaParam)
%regularizeSumstats regularizes a vector of Z scores or marginal effect-size
%estimates (in per-sd units) to make them conform with the LD patterns of
%an LDGM precision matrix P. It performs a ridge-like regularization:
% 
%   z_regularized = R * ((1-lambda) * R + lambda * I)^-1 * z, where R = P^-1
% 
%   Input arguments:
%   P: cell array of LDGM precision matrices
%   z: cell array of Z scores, merged with the LDGM using mergesnplists
%   whichIndices: output of mergesnplists
%   lambdaParam: ridge regularization parameter

noIndices = cellfun(@length, P);
x = precisionMultiply(P,z,whichIndices);
x = cellfun(@assignto, x, whichIndices, num2cell(noIndices), ...
    'UniformOutput', false);
y = cellfun(@mtimes, P, x, 'UniformOutput', false);
I_lambda_P = cellfun(@(p)lambdaParam*speye(size(p)) + (1 - lambdaParam) * p, ...
    P, 'UniformOutput',false);
y = cellfun(@mldivide,I_lambda_P,y,'UniformOutput',false);
y = cellfun(@(x,j)x(j), y, whichIndices, 'UniformOutput', false);
z_regularized = precisionDivide(P,y,whichIndices);

end