function [combinedEst] = combine_estimates(est, ref)
%combine_estimates combines h2MLE heritability estimates and calculates
% error and significance statistics
if nargin < 2
    ref = 1; % reference annotation for calculating enrichment
end

fields = fieldnames(est);
assert(length(intersect(fields,{'h2','annotSum','enrichment'}))==3)
combinedEst.h2 = sum(vertcat(est(:).h2));
combinedEst.annotSum = sum(vertcat(est(:).annotSum));
combinedEst.enrichment = combinedEst.h2*(combinedEst.annotSum(ref)/combinedEst.h2(ref))./combinedEst.annotSum;
combinedEst.params = mean(horzcat(est.params)')';
end

