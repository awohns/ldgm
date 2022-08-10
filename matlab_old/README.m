% This directory contains MATLAB code to estimate LDGM precision matrices
% and use apply them to GWAS data. Subdirectories are:
% 
%   examples: new users should start here. There are three example
%   scripts, which should be sufficient to get up and running using LDGMs.
%       
%   utility: functions for data processing and basic matrix operations. In
%   particular, we *strongly* recommend using the mergesnplists.m function
% 
%   blup: functions related to BLUPx-ldgm. This directory contains the
%   BLUPxldgm function as well as functions computing the likelihood of the
%   GWAS summary statistics under the BLUP (infinitesimal) model
% 
%   simulation: functions to simulate GWAS summary statistics
%   
%   precision: functions related to precision matrix estimation
% 
%   sparseinv: the sparseinv() MATLAB function written by Timothy A. Davis.
%   Please refer to sparseinv/Contents.m file for details.
% 
%   manuscript_scripts: use these scripts to replicate the results of our
%   manuscript [ref]. These are documented unevenly.
