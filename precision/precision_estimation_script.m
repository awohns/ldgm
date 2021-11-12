addpath(genpath('~/Documents/GitHub/ld_graph/precision'))
path_prefix='';
filename = '';

maxReps = 1e5;
convergenceTol=1e-5;

% number of SNPs in band for error estimation
bandsize = 100;

% Genotype matrix
X=load([path_prefix, filename, '.genos'])';

% Weighted adjacency matrix
import_weighted = 1;
A_weighted=importGraph([path_prefix, filename, '.adjlist'], import_weighted);

% Empty rows/columns correspond to duplicate SNPs (on same brick as
% another SNP)
SNPs = find(any(A_weighted));
X = X(:,SNPs);
if any(X(:)==-1)
    error('Missing genotypes not supported')
end
A_weighted = A_weighted(SNPs,SNPs);
[numHaplotypes, numSNPs] = size(X);
allele_freq = mean(X);
AF_threshold=0;

R = corr(X);
% [V, S] = eigs(R,1);
% R = R - V*S*V';

thresholds = [1:8];
density = zeros(length(thresholds),1); error = density; 
obj_val = density; converged = density;
banded_error = density;

precisionEstimate = speye(numSNPs);
tt=1;
for tt = tt:length(thresholds)
    A = A_weighted + speye(numSNPs) >= 1/(1 + thresholds(tt));
    if nnz(A) >= nnz(A_weighted)
        warning('No edges were discarded at threshold %f', thresholds(tt))
    end
    density(tt) = full(mean(A(:)));
    fprintf('Density at threshold %.1f: %.3f\n', thresholds(tt), density(tt))

    % fit with threshold tt, starting at precisionEstimate
    [precisionEstimate, obj_val(tt), converged(tt), data{tt}] =...
        LDPrecision(R, 'graphical_model', A, ...
        'P0', precisionEstimate, 'max_steps', maxReps,...
        'convergence_tol', convergenceTol, 'printstuff', 2);
    
    all_precision_estimates{tt} = precisionEstimate;
    
    % Plotting
    if 0
        plotting_script;drawnow
    else
        Rr = inv(precisionEstimate);
        MSE = mean((Rr(:) - R(:)).^2) / mean(R(:).^2);
        fprintf('Error at threshold %.1f: %.3f\n', thresholds(tt), MSE)

    end
    error(tt) = MSE;
    banded_error(tt) = mean((spdiags(R,1:bandsize)-spdiags(Rr,1:bandsize)).^2) / ...
        mean(spdiags(R,1:bandsize).^2);
    fprintf('Banded error: %.3f\n', banded_error(tt))
    
    
end

figure;

subplot(1,4,1)
scatter(density,error,'filled')
xlim([0 0.5]); ylim([0 max(density)])
xlabel('Density');ylabel('Error')

subplot(1,4,2)
title(filename,'interpreter','none'); hold on
scatter(thresholds,error,'filled')
ylim([0 .5])
xlabel('Threshold');ylabel('Error')

subplot(1,4,3)
scatter(thresholds,banded_error,'filled')
xlabel('Threshold');ylabel('Banded error')
ylim([0 .5])

subplot(1,4,4)
scatter(thresholds,obj_val,'filled')
xlabel('Threshold');ylabel('Objective value')
