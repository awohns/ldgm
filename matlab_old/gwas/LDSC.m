function  LDSCout = LDSC( chisq,l2,annot,varargin)
%LDSC (cross-trait LD score regression).
%
%   REQUIRED INPUT ARGUMENTS: chisq: Mx1 vector of chi^2 stats or squared
%   per-norm-genotype effect size estimates;
%   l2: MtotxP matrix of LD scores (LD second moments), where P is #annotations
%   and Mtot is the number of reference SNPs;
%   Optionally, prepend an all-ones vector in first column of l2,
%   which will cause an intercept to be used in the regression.
%   annot: MtotxP matrix of annotation values;

%   OPTIONAL INPUT ARGUMENTS as name-value pairs:
%   RegressionIndices: Mx1 vector of indices corresponding to the reference SNPs that
%   will be used in the regression (regression SNPs must be a subset of
%   reference SNPs);
%   l2Weights: Mx1 vector of weights for LDSC;
%   JackknifeBlocks: cell vector of custom jackknife blocks, corresponding to independent
%   loci. Primarirly for use with MixedModel option. Each cell contains a
%   list of indices between 1 and Mtot. All integers between 1 and Mtot
%   should be contained in exactly one cell.
%   NoJackknifeBlocks: number of jackknife blocks (default 100 or
%   length(JackknifeBlocks) if specified);
%   OutputSubspace: which annotations to report estimates for. This can be
%   either a list of P' indices <=P, or a PxP' projection matrix whose columns
%   correspond to linear combinations of input annotations.
%   MixedModel: pass "1" to treat large-effect SNPs as fixed effects.
%   Highly recommended to also specify JackknifeBlocks.
%   MixedModelThreshold:
%   SNPs with larger effects than threshold will be treated as fixed
%   effects.
%   NoBootstrapIter: number of bootstrap iterations (by default, block
%   jackknife is used instead)
%   SeedBootstrap: random seed for bootstrap for reproducibility.
%   default: 0->shuffle
%
%   OUTPUT ARGUMENTS:
%   h2est: "heritability", if input summary stats are in
%   appropriate units. Please be careful about units when reporting total
%   heritability estimates.
%   h2err: heritability jackknife standard error.
%   All other outputs are matrices of size NoJacknifeBlocks x
%   P', ie each row = a leave-one-block-out estimate and each column = an
%   annotation.
%   h2: "heritability"
%   tau: regression coefficients
%   intercept: LDSC intercept.


mm_regression=length(chisq);
[mm_tot,pp_l2]=size(annot);

checkvector=@(x)isreal(x) && isvector(x) && length(x)==mm_regression;
checkmatrix=@(x)isreal(x) && size(x,1)==mm_tot && size(x,2)==pp_l2;
checkmatrix_l=@(x)isreal(x) && size(x,1)==mm_regression && (size(x,2)==pp_l2 || size(x,2)==pp_l2+1);
checkIndexlist=@(x)(isvector(x) && all(mod(x,1)==0) && all(x<=pp_l2) && all(x>0));
checkProjectionmatrix=@(x)isreal(x) && (size(x,1)==pp_l2 && rank(x)==size(x,2));
is2x2CovarianceMatrix= @(x)size(x,1)==2 && size(x,2)==2 && x(1,2)==x(2,1) && x(1,2)<=sqrt(x(1,1)*x(2,2));
p=inputParser;
asColumn=@(x)reshape(x,numel(x),1);
checkJackknifeBlocks=@(x)iscell(x) && all(sort(vertcat(asColumn(x{:})))==asColumn(1:mm_tot));

addRequired(p,'chisq',@(x)true)
addRequired(p,'l2',checkmatrix_l)
addRequired(p,'annot',checkmatrix)
addOptional(p,'RegressionIndices',1:mm_tot,@(x)checkvector(x) && all(mod(x,1)==0) && all(x<=mm_tot))
addOptional(p,'l2Weights',[],checkvector)
addOptional(p,'JackknifeBlocks',[],@(x)checkJackknifeBlocks(x))
addOptional(p,'NoJackknifeBlocks',100,@(x)isscalar(x) && all(mod(x,1)==0) && all(x<=mm_regression) && all(x>1))
addOptional(p,'NoBootstrapIter',0,@(x)isscalar(x) && all(mod(x,1)==0) && all(x>=0))
addOptional(p,'SeedBootstrap',0,@(x)isscalar(x) && all(mod(x,1)==0) && all(x>=0))
addOptional(p,'OutputSubspace',[],@(x)checkProjectionmatrix(x) || checkIndexlist(x))
addOptional(p,'NoiseVar',[],@(x)(isscalar(x) | all(size(x) == [mm_regression 1])) & all(x>=0));
addOptional(p,'MixedModel',0,@(x)isscalar(x) && (x==0 || x==1));
addOptional(p,'MixedModelThreshold',100,@(x)isscalar(x) && x>0 );

parse(p,chisq,l2,annot,varargin{:});

if checkIndexlist(p.Results.OutputSubspace)
    output_subspace=sparse(p.Results.OutputSubspace,1:length(p.Results.OutputSubspace),ones(length(p.Results.OutputSubspace),1),pp_l2,length(p.Results.OutputSubspace));
else
    output_subspace=p.Results.OutputSubspace;
end

no_blocks=p.Results.NoJackknifeBlocks;

idx_regression=p.Results.RegressionIndices;

if isempty(p.Results.l2Weights)
    l2weights=1./max(1,l2(idx_regression,1));
else
    l2weights=p.Results.l2Weights;% 2nd-moment weights matrix
end
WW2_mat=diag(sparse(l2weights));

NoiseVar=p.Results.NoiseVar;

% Bootstrap options
run_bootstrap=p.Results.NoBootstrapIter>0;
if  run_bootstrap
    no_iters=p.Results.NoBootstrapIter;
    if p.Results.SeedBootstrap==0
        rng('shuffle');
    else
        rng(p.Results.SeedBootstrap)
    end
else
    no_iters=no_blocks;
end


% Annotations (or linear combinations thereof) for which to report estimates
if isempty(output_subspace);  output_subspace=eye(pp_l2);end
output_annot=annot*output_subspace;

% Jackknife blocks (ref SNPs)
if checkIndexlist(p.Results.JackknifeBlocks)
    jackknife_blocks=p.Results.JackknifeBlocks;
    no_blocks=length(jackknife_blocks);
else
    blocksize=floor(mm_tot/no_blocks);
    jackknife_blocks=cell(no_blocks,1);
    for jk=1:no_blocks
        jackknife_blocks{jk}=(jk-1)*blocksize+1:jk*blocksize;
    end
end

% Jackknife blocks (regression SNPs)
regression_blocks=jackknife_blocks;block_weights=jackknife_blocks;
for jk=1:no_blocks
    regression_blocks{jk}=idx_regression(jackknife_blocks{jk});
    block_weights{jk}=sum(l2weights(regression_blocks{jk}));
end

% Mixed model setup
fixedeffect_blocks=false(no_blocks,1);
if p.Results.MixedModel
    MixedModelThreshold=p.Results.MixedModelThreshold;
    
    for ii=1:no_blocks
        if any(chisq(regression_blocks{ii})>MixedModelThreshold)
            fixedeffect_blocks(ii)=true;
        end
    end
    fprintf('Treating %d out of %d blocks as fixed effects\n',...
        sum(fixedeffect_blocks),no_blocks);
    if ~run_bootstrap
        no_iters=no_iters-sum(fixedeffect_blocks);
    end
    
end

clear p varargin

% Handling LD score intercept options
if size(l2,2)==pp_l2 && isempty(NoiseVar)
    l2=[ones(mm_regression,1) l2];
end
if size(l2,2)==pp_l2+1 && any(l2(:,1)~=1)
    warning('If size(l2,2)==size(annot,2)+1, then first column of l2 should usually be all ones')
end
if size(l2,2)==pp_l2+1 && ~isempty(NoiseVar)
    warning('If size(l2,2)==size(annot,2)+1, then NoiseVar input is ignored')
end
if isscalar(NoiseVar)
    NoiseVar = ones(mm_regression,1) * NoiseVar;
end

l2_intercept=size(l2,2)-pp_l2;

pp_l2=size(l2,2);
pp_annot=size(annot,2);
pp_output=size(output_annot,2);


ell_ell=zeros(pp_l2,pp_l2,no_blocks);

annot_annot=zeros(pp_annot,pp_output,no_blocks);


ell_annot=zeros(pp_annot,pp_output);

ell_a1=zeros(pp_l2,no_blocks);
ell_mean=ell_a1;


% Compute various moments for each jackknife block
mean_weight=sum(l2weights)/no_blocks;
for jk=1:no_blocks
    ell_ell(:,:,jk)=l2(regression_blocks{jk},:)'*WW2_mat(regression_blocks{jk},regression_blocks{jk})*l2(regression_blocks{jk},:)/mean_weight;
    ell_a1(:,jk)=l2(regression_blocks{jk},:)'*WW2_mat(regression_blocks{jk},regression_blocks{jk})*chisq(regression_blocks{jk})/mean_weight;%
    ell_mean(:,jk)=mean(WW2_mat(regression_blocks{jk},regression_blocks{jk})*l2(regression_blocks{jk},:))/mean_weight;
    annot_annot(:,:,jk)=annot(jackknife_blocks{jk},:)'*output_annot(jackknife_blocks{jk},:)./max(1,sum(output_annot(jackknife_blocks{jk},:)));
    ell_annot(:,:,jk)=l2(regression_blocks{jk},l2_intercept+1:end)'*output_annot(jackknife_blocks{jk},:)./max(1,sum(output_annot(jackknife_blocks{jk},:)));
end

% Combine jackknife estimates to obtain leave-one-out regression
% coefficients

ind=cell(no_iters,1);
included_blocks=find(~fixedeffect_blocks);
for jk=1:no_iters
    if run_bootstrap
        ind{jk}=randsample(included_blocks,length(included_blocks),true);
    else
        ind{jk}=included_blocks([1:jk-1,jk+1:no_iters]);
    end
end
tau1=zeros(pp_annot,no_iters);
intercept1=zeros(no_iters,1);
beta2tot1=zeros(pp_output,no_iters);

for jk=1:no_iters
    
    %   S-LDSC
    if l2_intercept==1
        temp=mean(ell_ell(:,:,ind{jk}),3)\mean(ell_a1(:,ind{jk}),2);% l2var^-1 * l2a2
        tau1(:,jk)=temp(l2_intercept+1:end);
        intercept1(jk)=temp(1);
        
    else
        temp=mean(ell_ell(:,:,ind{jk}),3)\mean(ell_a1(:,ind{jk})-NoiseVar(ind{jk})'.*ell_mean(:,ind{jk}),2);% l2var^-1 * l2a2
        tau1(:,jk)=temp(l2_intercept+1:end);
        intercept1(jk)=mean(NoiseVar(ind{jk}));
    end
    
    % Second moments of effect-size distribution
    beta2tot1(:,jk)=tau1(:,jk)'*mean(annot_annot(:,:,ind{jk}),3);
    
end

% fixed-effect h2
fixedeffect_blocks=find(fixedeffect_blocks);
fixedeffect_h2=zeros(length(fixedeffect_blocks),1);
for jk=1:length(fixedeffect_blocks)
    fixedeffect_h2(jk)=max(chisq(regression_blocks{fixedeffect_blocks(jk)}))...
        -mean(intercept1);
end

% Heritabilities (ie M*E(beta^2))
h21_jk=beta2tot1.*sum(output_annot)';

LDSCout(1).h2=h21_jk;
LDSCout(1).tau=tau1;
LDSCout(1).intercept=intercept1;
LDSCout(1).h2est=mean(h21_jk);
LDSCout(1).h2err=std(h21_jk)*sqrt(no_blocks+1);
LDSCout(1).fixedeffecth2=fixedeffect_h2;

end




