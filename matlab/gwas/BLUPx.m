function [betaExpectation] =...
    BLUPx(alphaHat, P, nn, betaCov, SNPs, whichSNPs, AF, alpha_param)
% BLUPx computes the cross-popn best linear unbiased predictor, E(beta|GWAS, gaussian
% prior) where the GWAS data is alphaHat, the LD precision matrix is P, and
% beta~N(0,tau^-1).
% 
% Input arguments
% 
% alphaHat: GWAS sumstats, as a number-of-popns by number-of-LD-blocks cell
% array with each cell containing an association vector; 
% 
% P: LD precision matrices, as a number-of-popns by number-of-LD-blocks cell
% array with each cell containing a precision matrix; 
% 
% nn: GWAS sample size for each population, as a vector;
% 
% betaCov: covariance matrix for the per-allele effect size of a SNP across 
% popns, which is assumed to be i.i.d.;
% 
% SNPs: (optional) a number-of-popns by number-of-LD-blocks cell
% array with each cell containing a list of SNPs or indices, one for each
% row/column of the corresponding precision matrix. This tells the method
% which SNPs across popn-specific LD matrices correspond to each other.
% 
% whichSNPs: (optional) a number-of-popns by number-of-LD-blocks cell
% array with each cell containing a list of indices, such that the SNPs
% corresponding to alphaHat{i,j}(:) are those corresponding to 
% P{i,j}(SNPs{i,j}, SNPs{i,j})
% 
% AF: (optional) a number-of-popns by number-of-LD-blocks cell
% array with each cell containing allele frequencies. The SNPs in each cell
% should correspond to those in alphaHat.
% 
% Output arguments:
% betaExpectation: expected value of beta for each popn. Reported for same
% SNPs as alphaHat.

[noBlocks, noPopns] = size(P);
mm = cellfun(@length, alphaHat);

% turn whichSNPs into boolean vectors if needed
if isa(whichSNPs{1},'boolean')
    assert(all(cellfun(@(x,y)length(x)==length(y),whichSNPs,P)))
else
    whichSNPs = cellfun(@(ii,X){unfind(ii,length(X))},whichSNPs,P);
end

% turn SNPs into vectors of indices if needed
if isa(SNPs{1},'boolean')
    SNPs = cellfun(@find,SNPs,'UniformOutput',false);
end
assert(all(cellfun(@(x,y)length(x)==length(y),SNPs,P)))

% effect-size s.d. for each SNP in each population, in standardized
% (per-sd-of-genotype) units. alpha_param is the AF-dependent architecture parameter
% of e.g. Schoech et al. 2019. For a SNP that is missing in a population, 
% its effect-size s.d. is zero.
if nargin < 7
    sd = cellfun(@double,whichSNPs,'uniformoutput',0);
else
    
    if nargin < 8
         % default no AF-dependent architecture (constant per-allele effect
         % size variance)
        alpha_param = 0;
    end
    % assign sqrt(2pq) (if alpha_param==0) to nonmissing SNPs
    sd = cellfun(@(p,tf){assignto((2*p.*(1-p)) .^ ((alpha_param+1)/2), tf)},...
        AF, whichSNPs);
    
end

% Causal effect-size estimates
betaHat = precisionMultiply(P,alphaHat,whichSNPs);

% Concatenated effect-size estimates and precision matrices across popns
betaHatCat = cell(noBlocks,1);
PCat = cell(noBlocks,1);
whichSNPsCat = cell(noBlocks,1);
nn_cell = num2cell(nn);
for block = 1:noBlocks
    betaHatCat{block} = vertcat(betaHat{block, :});
    Pblocks = cellfun(@(x,y)y/x, nn_cell, P(block,:),'UniformOutput',false);
    PCat{block} = blkdiag(Pblocks{:});   
    whichSNPsCat{block} = vertcat(whichSNPs{block, :});
end

% Covariance matrix of effect sizes for each concatenated block
Sigma = cell(noBlocks,1);
for block = 1:noBlocks
    S = cell(noPopns); % blocks of covariance matrix
    for ii = 1:noPopns
        for jj = 1:noPopns
            [~, i1, i2] = intersect(SNPs{block,ii}, SNPs{block,jj});
            S{ii,jj} = sparse(i1,i2,...
                betaCov(ii,jj).*sd{block,ii}(i1).*sd{block,jj}(i2),...
                length(SNPs{block,ii}),length(SNPs{block,jj}));
        end
    end
    Sigma{block} = cell2mat(S);
end

% E(beta|data)
x = precisionDivide(cellfun(@(p,s){p+s},PCat,Sigma),...
    betaHatCat, whichSNPsCat);
betaExpectationCat = cellfun(@(s,v,idx)s(idx,idx)*v,Sigma,x,whichSNPsCat,...
    'UniformOutput',false);%precisionMultiply(Sigma, x, whichSNPsCat);

% Un-concatenate populations
betaExpectation = cell(size(alphaHat));
for block = 1:noBlocks
    betaExpectation(block,:) = mat2cell(betaExpectationCat{block}, mm(block,:));
end



end

