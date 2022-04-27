function [betaExpectation] =...
    BLUPx(alphaHat, P, nn, betaCov, SNPs, whichSNPs)
% BLUPx computes the cross-popn best linear unbiased predictor, E(beta|GWAS, gaussian
% prior) where the GWAS data is alphaHat, the LD precision matrix is P, and
% beta~N(0,tau^-1).
% Input arguments
% alphaHat: GWAS sumstats, as a number-of-popns by number-of-LD-blocks cell
% array with each cell containing an association vector; 
% P: LD precision matrices, as a number-of-popns by number-of-LD-blocks cell
% array with each cell containing a precision matrix; 
% nn: GWAS sample size for each population, as a vector;
% betaCov: covariance matrix for the effect size of each SNP across popns,
% which is assumed to be i.i.d.;
% SNPs: (optional) a number-of-popns by number-of-LD-blocks cell
% array with each cell containing a list of SNPs or indices, one for each
% row/column of the corresponding precision matrix. This tells the method
% which SNPs across popn-specific LD matrices correspond to each other.
% whichSNPs: (optional) a number-of-popns by number-of-LD-blocks cell
% array with each cell containing a list of indices, such that the SNPs
% corresponding to alphaHat{i,j}(:) are those corresponding to 
% P{i,j}(SNPs{i,j}, SNPs{i,j})
% Output arguments:
% betaExpectation: expected value of beta for each popn given the data
% h2MLE: vector of MLE values for diag(beta*alpha'). Can be negative. Add
% these up to get an MLE heritability estimate for a set of SNPs.

[noBlocks, noPopns] = size(P);

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

% Causal effect-size estimates
betaHat = precisionMultiply(P,alphaHat,whichSNPs);

% Concatenated effect-size estimates and precision matrices across popns
betaHatCat = cell(noBlocks,1);
PCat = cell(noBlocks,1);
whichSNPsCat = cell(noBlocks,1);
for block = 1:noBlocks
    betaHatCat{block} = vertcat(betaHat{block, :});
    PCat{block} = blkdiag(P{block,:});   
    whichSNPsCat{block} = vertcat(whichSNPs{block, :});
end

% Covariance matrix of effect sizes for each concatenated block
Tau = cell(noBlocks,1);
for block = 1:noBlocks
    S = cell(noPopns); % blocks of covariance matrix
    for ii = 1:noPopns
        for jj = 1:noPopns
            [~, i1, i2] = intersect(SNPs{block,ii}, SNPs{block,jj});
            S{ii,jj} = sparse(i1,i2,betaCov(ii,jj),...
                length(SNPs{block,ii}),length(SNPs{block,jj}));
        end
    end
    Tau{block} = cell2mat(S);
end

% for each cell, x = (P/nn + Sigma) \ betaHat
x = precisionDivide(cellfun(@(p,s,n){p/n+s},PCat,Sigma,nn),...
    betaHat,whichSNPsCat);

% double-check this
betaExpectation = precisionMultiply(Sigma,x,whichSNPs);
    

function lgc = unfind(idx, N)
    %Go from indicies into logical (for vectors only)
    lgc = false(N,1);
    lgc(idx) = true;
end
end

