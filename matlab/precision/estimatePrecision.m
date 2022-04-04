function P = estimatePrecision(genos_path, varargin)
% Estimates precision matrices for EUR superpopulation using LDGMs derived
% from union of NYGC superpopulations

p=inputParser;
% directory containing genotype matrices
addRequired(p, 'genos_dir', @isstr);

% directory containing adjacency matrices, defaults to genos_dir
addParameter(p, 'adjlist_dir', genos_path, @isstr);

% instead of using adjacency matrix, create adjacency matrix with nearest
% neighbors. Number of neighbors will be chosen based on adjacency matrix
% and path_distance_threshold
addParameter(p, 'banded_control_ldgm', 0, @isscalar);

% directory containing SNP list, defaults to genos_dir
addParameter(p, 'snplist_dir', genos_path, @isstr);

% regular expression for genotype files in genos_dir
addParameter(p, 'genos_filename', '*.genos', @isstr);

% among files matching [genos_dir genos_filename], which one to use
addParameter(p, 'genos_file_index', 1, @isscalar);

% directory for saving
addParameter(p, 'output_dir', '', @isstr);

% file containing population data for each sample in genotype matrix
addParameter(p, 'population_data_file', '', @isstr);

% name of population for which to infer precision matrix, eg 'EUR'
addParameter(p, 'population_name', 'ALL', @isstr);

% retain LDGM edges with path distance < threshold 
addParameter(p, 'path_distance_threshold', 4, @isscalar);

% L1 penalty parameter
addParameter(p, 'l1_penalty', 0.05, @isscalar);

% verbosity for gradient descent function (0,1,2)
addParameter(p, 'verbose', 1, @isscalar);

% maximum number of gradient descent steps with L1 penalty
addParameter(p, 'max_gradient_steps_penalized', 2e4, @isscalar);

% maximum number of gradient descent steps without L1 penalty
addParameter(p, 'max_gradient_steps_unpenalized', 2e5, @isscalar);

% convergence tolerance for gradient descent
addParameter(p, 'gradient_convergence_tolerance', 1e-6, @isscalar);

% minimum MAF filter
addParameter(p, 'minimum_maf', 0.05, @isscalar);

% band size for error estimation
addParameter(p, 'bandsize', 200, @isscalar);

% turns p.Results.x into just x
parse(p, genos_path, varargin{:});
fields = fieldnames(p.Results);
for k=1:numel(fields)
    line = sprintf('%s = p.Results.%s;', fields{k}, fields{k});
    eval(line);
end

% suffixes for output files
output_matrix_suffix = '.precisionMatrix';
output_stats_suffix = '.stats';
output_suffix = sprintf('.path_distance=%.1f.l1_pen=%.2f.maf=%.2f.%s',...
    path_distance_threshold,l1_penalty,minimum_maf,population_name);

% Threshold to decide if edge weights are zero
smallNumber = 1e-12;

addpath(genpath('/broad/oconnor/trees/MATLAB'))

% genotype file
genos_files = dir([genos_dir,genos_filename]);
assert(~isempty(genos_files));
genos_filename = genos_files(genos_file_index).name;
assert(strcmp(genos_filename(end-5:end),'.genos'));
filename = genos_filename(1:end-6);

% genotype matrix
X = load([genos_dir, genos_filename])';
assert(all(unique(X(:)) == [0 1]'));

% which samples to retain
if ~isempty(population_data_file) && ~strcmp(population_name,'ALL')
    T = readtable(population_data_file);
    superpops = table2cell(T(:,3));
    assert(length(superpops) == size(X,1),'Number of rows in population_data_file should match number of samples in .genos file')
    rows = cellfun(@(s)strcmp(s,population_name),superpops);
    assert(~isempty(rows))
    X = X(rows,:);
end
noSNPs = size(X,2);

% SNP list file
snplist_files = dir([snplist_dir,filename,'.snplist']);
if isempty(snplist_files)
    warning("SNP list not found in [snplist_dir,filename,'.snplist']")
    snpTable = 1:noSNPs;
else
    snpTable = readtable([snplist_dir,filename,'.snplist'],'FileType','text');
end

% adjacency file
adjlist_files = dir([adjlist_dir,filename,'.adjlist']);
assert(~isempty(adjlist_files),"Adjacency list file not found in [adjlist_dir,filename,'.adjlist']");

% adjacency matrix
A_weighted=importGraph([adjlist_dir, filename, '.adjlist'], 1, noSNPs);

% Empty rows/columns correspond to duplicate SNPs (on same brick as
% another SNP); also get rid of LF SNPs
SNPs = any(A_weighted) .* (min(mean(X), 1-mean(X))>minimum_maf) == 1;

A_weighted = A_weighted(SNPs,SNPs);
X = X(:,SNPs);
SNPs = find(SNPs);
snpTable = snpTable(SNPs,:);

[noHaplotypes, noSNPs] = size(X);

% LD correlation matrix
R = corr(X);
R_band = spdiags(R,1:bandsize);

% Initial graphical model
A = A_weighted + speye(noSNPs) >= 1/(1 + path_distance_threshold);
initialDegree = ceil(nnz(A) / (2*noSNPs));

if banded_control_ldgm
    A = spdiags(true(noSNPs,2*initialDegree+1),-initialDegree:initialDegree,noSNPs,noSNPs);
    output_suffix = [output_suffix, '.banded_control'];
end

% estimate precision matrix with L1 penalty
tic;
[precisionEstimatePenalized, ~, converged_penalized, convergence_data_penalized] =...
    LDPrecision(R, 'graphical_model', A, ...
    'P0', speye(noSNPs), 'max_steps', max_gradient_steps_penalized,...
    'convergence_tol', gradient_convergence_tolerance, 'printstuff', verbose, 'lambda', l1_penalty);
time_gd_penalized=toc;
    
% Updated graphical model
A = abs(precisionEstimatePenalized) > smallNumber;
avgDegree = nnz(A) / length(A);

% initialize at regularized estimate
precisionEstimatePenalized(~A)=0;
[~,a] = chol(precisionEstimatePenalized);
while a > 0 % precisionEstimatePenalized is not PSD
    precisionEstimatePenalized = precisionEstimatePenalized + 0.01 * speye(noSNPs);
    [~,a] = chol(precisionEstimatePenalized);
end

% Estimate precision matrix at edges selected with L1 penalty, dropping the
% penalty term
[precisionEstimate, ~, converged, convergence_data] =...
    LDPrecision(R, 'graphical_model', A, ...
    'P0', precisionEstimatePenalized, 'max_steps', max_gradient_steps_unpenalized,...
    'convergence_tol', gradient_convergence_tolerance, 'printstuff', 1, 'lambda', 0);
time_gd_unpenalized=toc;

% Error metrics
A_band = spdiags(A,1:bandsize);
Rr = inv(precisionEstimate);
Rr_band = spdiags(Rr,1:bandsize);
banded_error = mean((R_band(:) - Rr_band(:)).^2) / ...
    mean(R_band(:).^2);
banded_error_noA = mean((R_band(~A_band) - Rr_band(~A_band)).^2) / ...
    mean(R_band(~A_band).^2);
error = mean((Rr(:) - R(:)).^2) / mean(R(:).^2);
error_A = mean((R(A) - Rr(A)).^2) / mean(R(A).^2);

mkdir(output_dir);
save([output_dir,filename,output_suffix,output_stats_suffix,'.mat'],...
    '*error*','avgDegree','initialDegree','varargin','convergence_data*','converged*',...
    'noSNPs','SNPs','time*')

save([output_dir,filename,output_suffix,output_matrix_suffix,'.mat'],...
    'precisionEstimate*','snpTable','A_weighted','-v7.3')

saveGraph([output_dir,filename,output_suffix,'.adjlist'],precisionEstimate,SNPs);


