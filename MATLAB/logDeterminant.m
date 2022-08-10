function [logDetM, A] = logDeterminant(P,whichIndices,d)
% logDeterminant computes the log-determinant of the matrix
% M = R(whichIndices,whichIndices) + D,
% where R == inv(P) is the correlation matrix, D == diag(d) is a diagonal
% matrix, and whichIndices is a list of indices.
% 
% Input arguments:
% P: LDGM precision matrices, as a cell array with each cell containing a 
% precision matrix or as a single precision matrix
% 
% whichIndices (optional): output from mergesnplists, encoding which indices
% (rows/columns of the LDGM precision matrices) have corresponding SNPs in
% each summary statistic file respectively. If P is a cell array, this
% should be a cell array of the same size; if P is a single matrix, this
% should be a vector of logicals or indices
% 
% d (optional): vectors to be converted into a diagonal matrix and added to
% P. Can also specify scalars. If P is a cell array, this
% should be a cell array of the same size; if P is a single matrix, this
% should be a scalar or a vector of size equal to the number of nonmissing
% SNPs
%
% Output arguments:
% logDetM: log-determinant of M for each LD block
% 
% A: cholesky factor of L, which has submatrix M; has size equal to number
% of nonzero elements on the diagonal of P

if iscell(P)
    if nargin < 2
        whichIndices = cellfun(@(P_el)diag(P_el)~=0,P,'UniformOutput',false);
    end
    if nargin < 3
        d = num2cell(zeros(size(P)));
    end
    assert(iscell(whichIndices) && ...
        (iscell(d)), ...
        'If P is a cell array, then d and whichSNPs should also cell arrays');

    assert(all(size(P) == size(d)) && all(size(P) == size(whichIndices)),...
        'Input cell arrays should have same size');
    
    % Iterate over cell arrays
    logDetM = cellfun(@(P_el,whichIndices_el,d_el)logDeterminant(P_el,whichIndices_el,d_el), ...
        P, whichIndices, d,'UniformOutput',false);
    
else
    noSNPs = length(P);
    
    % L has submatrix M==L(whichIndices,whichIndices)
    if any(d<0)
        warning('d has negative elements, which is probably not desired')
    end
    L = sparse(zeros(noSNPs,1));
    L(whichIndices) = d;
    L = diag(L) + P;
    
    % handling SNPs missing from P
    if ~islogical(whichIndices)
        whichIndices = unfind(whichIndices,noSNPs);
    end
    incl = diag(P)~=0;
    otherSNPs = incl & (~whichIndices);
    assert(all(incl(whichIndices)), 'P should have nonzero diagonal entries for every nonmissing SNP')
    
    A = chol(L(incl,incl));
    
    % log|M| == log|L| - log|P(otherSNPs,otherSNPs)|
    logDetL = 2*sum(log(diag(A)));
    logDetPOther = 2*sum(log(diag(chol(P(otherSNPs,otherSNPs)))));
    logDetM = logDetL - logDetPOther;
end