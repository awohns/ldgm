% This directory contains MATLAB code to infer LDGM precision matrices and
% to analyze GWAS data using LDGMs.
% 
% New users should start in the tutorials directory, where there are
% tutorials on how to analyze real and simulated GWAS summary statistics
% using LDGM precision matrices.
% 
% This directory contains the following functions:
% 
%   BLUPxldgm: computes polygenic scores from GWAS summary statistics and
%   LDGM precision matrices, for one or more populations
%
%   BLUPxldcov: computes polygenic scores from GWAS summary statistics and
%   LD correlation matrices, for one or more populations
% 
%   GWASlikelihood: computes log-likelihood of the GWAS summary statistics
%   under a Gaussian model
% 
%   GWASlikelihoodGradient: this function requires the sparseinv function,
%   which can be obtained by installing suitesparse: 
%   https://people.engr.tamu.edu/davis/suitesparse.html
%   It calculates the gradient of the log-likelihood w.r.t. the variance
%   parameters of the Gaussian prior distribution.
% 
%   loadLDGMs: loads LDGM precision matrices and SNP lists from a specified
%   directory
% 
%   logDeterminant: calculates the log-determinant of a matrix M=R+D, where
%   R is the inverse of the LDGM precision matrix and D is a diagonal
%   matrix.
% 
%   mergesnplists: merges LDGM precision matrix SNPs with those from a
%   table of GWAS summary statistics. It also takes care of allele phasing.
%   We highly recommend using this function, as documented in the tutorial
% 
%   precisionDivide: divides a vector y by the precision matrix, accounting
%   for missing SNPs (i.e., multiplies y by the inverse of the precision
%   matrix restricted to nonmissing SNPs)
% 
%   precisionMultiply: multiplies a vector y by the precision matrix
%   accounting for missing SNPs
% 
%   r2PGS: calculates the PGS r2 of a polygenic predictor in a population
%   with the specified LD precision matrix and allele frequencies
% 
%   simulateSumstats: simulates GWAS summary statistics with the specified
%   LD patterns and effect size distribution, for one or more populations
%   
%   ImpGldgm and ImpGldcov: perform linear imputation of GWAS summary
%   statistics using an LDGM precision matrix or an LD correlation matrix,
%   respectively
% 
% The utility subdirectory contains functions that users will probably not
% use directly, but which are called from the functions above.
