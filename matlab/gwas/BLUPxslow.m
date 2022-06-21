function [betaExpectationPerSD, betaExpectationPerAllele] =...
    BLUPxslow(alphaHat, R, nn, betaCov, SNPs, whichSNPs, AF, alpha_param)
% BLUPxslow computes the cross-popn best linear unbiased predictor, E(beta|GWAS, gaussian
% prior) where the GWAS data is alphaHat, the LD correlation matrix is R, and
% beta~N(0,tau^-1).
% 
% Input arguments
% 
% alphaHat: GWAS sumstats, as a number-of-LD-blocks by number-of-popns cell
% array with each cell containing an association vector; 
% 
% R: LD correlation matrices, as a number-of-LD-blocks by number-of-popns cell
% array with each cell containing a correlation matrix of the same size as 
% corresponding entry of alphaHat
% 
% nn: GWAS sample size for each population, as a vector;
% 
% betaCov: covariance matrix for the per-allele effect size of a SNP across 
% popns, which is assumed to be i.i.d.;
% 
% SNPs: (optional) a number-of-LD-blocks by number-of-popns cell
% array with each cell containing a list of SNPs or indices, one for each
% row/column of the corresponding precision matrix. This tells the method
% which SNPs across popn-specific LD matrices correspond to each other.
% 
% whichSNPs: (optional) a number-of-LD-blocks by number-of-popns cell
% array with each cell containing a list of indices, such that the SNPs
% corresponding to alphaHat{i,j}(:) are those corresponding to 
% P{i,j}(SNPs{i,j}, SNPs{i,j})
% 
% AF: (optional) a number-of-LD-blocks by number-of-popns cell
% array with each cell containing allele frequencies. The SNPs in each cell
% should correspond to those in alphaHat.
% 
% Output arguments:
% betaExpectation: expected value of beta for each popn. Reported for same
% SNPs as alphaHat.

[noBlocks, noPopns] = size(R);
mm = cellfun(@length, alphaHat);

% subset correlation matrices to SNPs with sumstats
if exist('whichSNPs','var')
    assert(all(cellfun(@(x,y)length(x)==length(y),whichSNPs,R),'all'))
    for block = 1:noBlocks
        R{block} = R{block}(whichSNPs{block},whichSNPs{block});
        SNPs{block} = SNPs{block}(whichSNPs{block});
    end
end

% turn whichSNPs into boolean vectors if needed
if isa(whichSNPs{1},'boolean')
    assert(all(cellfun(@(x,y)length(x)==length(y),whichSNPs,R),'all'))
else
    whichSNPs = cellfun(@(ii,X){unfind(ii,length(X))},whichSNPs,R);
end

% turn SNPs into vectors of indices if needed
if isa(SNPs{1},'boolean')
    SNPs = cellfun(@find,SNPs,'UniformOutput',false);
end
assert(all(cellfun(@(x,y)length(x)==length(y),SNPs,R),'all'))

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

% Concatenated effect-size estimates and correlation matrices across popns,
% multiplied by popn-specific sample sizes
nnAlphaHatCat = cell(noBlocks,1);
nnRCat = cell(noBlocks,1);
nn_cell = num2cell(nn);
for block = 1:noBlocks
    nnAlphaHat = cellfun(@times, nn_cell, alphaHat(block,:),'UniformOutput',false);
    nnAlphaHatCat{block} = vertcat(nnAlphaHat{:});
    nnRblocks = cellfun(@times, nn_cell, R(block,:),'UniformOutput',false);
    nnRCat{block} = blkdiag(nnRblocks{:});   
end

% Precision matrix of beta_persd for each concatenated block
betaPrecision = inv(betaCov); % per-allele units
betaPrecision_persd = cell(noBlocks,1);
L = cell(noPopns); % population-pair blocks of precision matrix
for block = 1:noBlocks
    for ii = 1:noPopns
        for jj = 1:noPopns
            [~, i1, i2] = intersect(SNPs{block,ii}, SNPs{block,jj});
            L{ii,jj} = sparse(i1,i2,...
                betaPrecision(ii,jj)./(sd{block,ii}(i1).*sd{block,jj}(i2)),...
                length(SNPs{block,ii}),length(SNPs{block,jj}));
        end
    end
    betaPrecision_persd{block} = cell2mat(L);
end

% E(beta|data)
for block = 1:noBlocks
    betaExpectationCat{block} = (nnRCat{block} + betaPrecision_persd{block}) \ nnAlphaHatCat{block};
end

% Un-concatenate populations
betaExpectationPerSD = cell(size(alphaHat));
for block = 1:noBlocks
    betaExpectationPerSD(block,:) = mat2cell(betaExpectationCat{block}, mm(block,:));
end

% per-allele effect sizes
if nargout > 1
    betaExpectationPerAllele = cellfun( @(beta,af)beta./sqrt(2*af.*(1-af)), betaExpectationPerSD, AF,'UniformOutput',false);
end


end

