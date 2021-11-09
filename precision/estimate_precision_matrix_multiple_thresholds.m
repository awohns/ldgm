addpath(genpath('~/Documents/GitHub/ld_graph/precision'))
path_prefix='/Users/loconnor/Dropbox/Pouria/data/';
filename = '1kg_med';

maxReps = 1e5;
convergenceTol=1e-5;

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

thresholds = [1:.5:5, 6:8];
density = zeros(length(thresholds),1); error = density; obj_val = density; converged = density;


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
    [precisionEstimate, obj_val(tt), converged(tt)] = LDPrecision(R, 'graphical_model', A, ...
        'P0',precisionEstimate,'max_steps',maxReps,...
        'convergence_tol',convergenceTol,'printstuff',true);

    % Plotting
    if 0
        plotting_script;drawnow
    else
        Rr = inv(precisionEstimate);
        MSE = mean((Rr(:) - R(:)).^2) / mean(R(:).^2);
        fprintf('Error at threshold %.1f: %.3f\n', thresholds(tt), MSE)

    end
    error(tt) = MSE;
end

figure;

subplot(1,3,1)
scatter(density,error,'filled')
xlim([0 0.5]);ylim([0 .5])
xlabel('Density');ylabel('Error')

subplot(1,3,2)
title(filename,'interpreter','none'); hold on
scatter(thresholds,error,'filled')
ylim([0 .5])
xlabel('Threshold');ylabel('Error')

subplot(1,3,3)
scatter(thresholds,obj_val,'filled')
xlabel('Threshold');ylabel('Objective value')



