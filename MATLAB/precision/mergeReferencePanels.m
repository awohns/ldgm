function [R, T] = mergeReferencePanels(X1, X2, T1, T2, keep_X2_only, MAF_threshold)
%mergeReferencePanels combines reference data from two LD reference panels
%and computes a meta-analyzed correlation matrix
% Input arguments:
% X1, X2: genotype matrices for each reference panel
% T1, T2: a table for each reference panel. Must contain "SNP" field for
% merging, as well as "A1" and "A2". 
% keep_X2_only: whether to keep SNPs that only appear in the second
% reference panel
% MAF_threshold: minimum MAF in the combined reference panel
% Output arguemnts:
% R: correlation matrix from the combined sample
% whichReference: for each row/column of R, whether it was present in
% reference panel 1 only ("1"), 2 only ("2"), or both ("0")

n1 = height(T1); n2 = height(T2);

% Detect if data are haploid or diploid
if max(X1(:)) > 1
    X1 = X1/2;
end
if max(X2(:)) > 1
    X2 = X2/2;
end

[~, i1, i2] = intersect(T1.SNP, T2.SNP, 'stable');

phase = mergealleles(T1.A1(i1),T1.A2(i1),T2.A1(i2),T2.A2(i2));
if mean(phase == 0) > 0.01
    warning('%f of variants had mismatched alleles, possibly indicating strandedness issue', mean(phase==0))
end
i1 = i1(phase~=0); 
i2 = i2(phase~=0); 
phase = phase(phase~=0); 

% convert to phase of T1
X2(:, i2) = (1 - phase')/2 + phase' .* X2(:, i2);

i1only = setdiff(1:n1, i1); 

if ~keep_X2_only
    i2only = setdiff(1:n2, i2); 
else
    i2only = [];
end

mm = size(X2,1);
X2only = X2(:,i2only);
X2 = [X1(:, i1); X2(:, i2)];


X2 = [X2, [X1(:,i1only); repmat(mean(X1(:,i1only)), mm, 1)]];
X2 = [X2, [repmat(mean(X2only), size(X1,1), 1); X2only]];
T = [T1(i1,:); T1(i1only,:)];
if ~keep_X2_only
    T = [T; T2(i2only,:)];
end

AF = mean(X2); MAF = min(AF, 1-AF);
T.AF = AF';
T.whichReference = [zeros(length(i1),1); ones(length(i1only),1); 2 * ones(length(i2only),1)];

if exist('MAF_threshold')
    X2 = X2(:,MAF > MAF_threshold);
    T = T(MAF > MAF_threshold,:);
end

R = corr(X2);


end