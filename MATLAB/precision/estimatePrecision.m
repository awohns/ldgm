function [precisionEstimate, R, ldSnpTable] = estimatePrecision(data_directory, varargin)
% Estimates LDGM precision matrix from an LDGM and genotype data using
% modified DP-GLASSO algorithm
%
% Input data (located at data_directory) can be one of three types,
% indicated by 'data_type' input:
%   -genotypes: an individials-by-SNPs binary matrix with filename extension
%   .genos, and columns that correspond to the rows of .snplist file (see
%   below). If this is specified, there several additional options (see below)
%   -correlation: a SNPs-by-SNPs symmetric matrix with filename extension
%   .txt, encoding the correlation between each pair of SNPs.
%   -edgelist: an edges-by-3 .edgelist file encoding the correlation between
%   each pair of SNPs with an edge in the LDGM. If any values are missing,
%   the corresponding edge of the LDGM is dropped.
%
% Differences between these options are as follows. If genotypes are
% specified, the .genos file must match the SNPs in the .snplist file.
% There is no option to specify a second .snplist file and merge the two
% sets of SNPs. On the other hand, there are several additional options
% (see below). If correlation or edgelist is specified, it is required to
% also specify LD_matrix_snplist_dir, whose rows correspond to the SNPs in
% the LD matrix (or edgelist). This will automatically be merged with the
% LDGM .snplist file, using the column named 'site_ids'. The resulting 
% precision matrix will have one row/column for each labeled brick (i.e.,
% nonredundant SNP) in the LDGM such that at least one of its SNPs is
% present in the LD matrix. Optionally, specify 'AF_column_name' together
% with 'minimum_AF'. With this option, the allele frequencies will be
% looked up in the LD matrix .snplist file (rather than the LDGM .snplist
% file), under the column with the corresponding name. Optionally, specify
% 'A1_column_name' and 'A2_column_name'. With this option, the alleles will
% be matched between A1/A2 and ancestral/derived allele, and the resulting
% precision matrix will be phased to ancestral/derived.
%
% If genotypes are specified, there are various additional options:
%   -population_data_file: a file describing the ancestry of each
%   individual (column) of the .genos file. If this is specified, it can be
%   used to restrict to a subset of these individuals by ancestry, or to
%   compute a weighted meta-analysis across ancestries
%   -population_name: which population (in population_data_file) to restrict
%   to. This can also be specified with the other two options, but it
%   only affects which column of the .snplist file is used as the allele
%   frequencies (i.e., AF_POPULATIONNAME).
%   -meta_weights: meta-analysis weights to compute LD as a weighted
%   average across populations in the population_data_file. To use this
%   option, set population_name to 'META'
%   -downsample_fraction: option to randomly downsample individuals in the
%   .genos file. If this option is used, the .stats.txt output file will
%   contain some additional columns reporting cross-validated MSE.
%   -local_ancestry_dir: not yet supported
%
% FILENAME HANDLING OPTIONS:
%   -edgelist_dir: directory containing LDGM .edgelist files. Defaults to
%   data_directory.
%   -snplist_dir: directory containing .snplist files. Defaults to
%   data_directory.
%   -output_dir: directory in which to save the output files
%   -custom_filename: if specified, this is prepended to the name of the
%   output files within output_dir
%   -data_pattern: if specified, only data files matching this
%   expression will be found
%   -data_file_index: out of the data files that are found, which one to
%   compute a precision matrix for
%   -LD_matrix_snplist_dir: only applicable if specifying an LD matrix or
%   LD correlation matrix SNP list as the data file; directory in which to
%   find a corresponding .snplist file, for merging with the LDGM .snplist
%   file
% 
% MERGING RELATED OPTIONS:
% If LD_matrix_snplist_dir is specified, there are various options for how
% to merge between the LD information and the LDGM
%  -match_snps_on_position: if true (default), there must be a column named
%  'position' in both SNP lists, and this will be used for the merge
%  (together with the allele info). If false, there must be a column named
%  'site_ids' which will be used instead.
%  -A1_column_name, A2_column_name: supply these in order to tell the
%  method which columns contain A1 and A2 (it doesn't matter which is
%  which)
%  -AF_column_name: if supplied, this allele frequency will be used instead
%  of the allele frequency specified in the LDGM SNP list
%
% METHOD PARAMETERS: (all of these can probably be left at default values)
%   -minimum_maf: minimum allele frequency, as calculated either using the
%   LDGM .snplist file (default) or the LD matrix .snplist file (if
%   AF_column_name is specified)
%   -l1_penalty: L1 penalty for graphical lasso. Method performs some
%   coordinate descent steps with the penalty, throws out edges whose
%   precision matrix entries are zero, and continues without the penalty.
%   Default value is 0.1.
%   -path_distance_threshold: LDGM entries with distance below this
%   threshold will be retained. Default value is 4.
%   -patch_paths: how to remove missing and low-frequency SNPs. Options are
%   0 (discard these SNPs, do not patch paths that are broken), 1 (patch
%   paths for missing SNPs, not for low-frequency SNPs), and 2 (patch paths
%   for both missing and low-frequency SNPs) (not recommended). Default
%   value is 0.
%   -penalized_iterations: number of coordinate descent iterations to
%   perform for each SNP in the L1-penalized graphical lasso step. Default
%   value is 5.
%   -unpenalized_iterations: number of coordinate descent iterations to
%   perform for each SNP in the unpenalized graphical lasso step. Default
%   value is 25.
%   -lasso_iterations: number of inner-loop iterations in the coordinate
%   descent procedure. Default value is 10.
%   -banded_control_ldgm: instead of using the LDGM that is specified,
%   replace it with a banded-diagonal matrix. The band size is chosen to
%   match the density of the LDGM (with the specified
%   path_distance_threshold). This should produce worse results than the
%   original LDGM.
%   -r2_control_ldgm: similar to banded_control_ldgm, except that instead
%   of a banded diagonal matrix it uses a matrix with a '1' wherever the
%   r^2 between two SNPs is high enough. This should produce much worse
%   results than using the original LDGM.
%
% OUTPUT FILES: two output files are written to output_dir directory. They
% are named [custom_filename, filename, extension] where custom_filename is an
% optional user-specified string, filename is the name of the data file
% (except without its filename extension), and extension is either
% '.precisionMatrix.edgelist' or '.stats.txt'
%   -'*.stats.txt' is a table with one row. It contains information like
%   the block size, the runtime, and the MSE. Note: if the data file is a
%   SNP list, then the MSE will be calculated from whatever entries are
%   provided in that file, and it might be unreliable if too few are
%   provided.
%   -'*.precisionMatrix.edgelist' is an edge list with each row
%   corresponding to an entry of the precision matrix. This file is
%   zero-indexed, and the rows/columns of the precision matrix correspond
%   to the 'index' column of the LDGM .snplist file.
%
% Instead of saving the output, it is also an option to return the
% precision matrix P directly.


p=inputParser;

% directory containing genotype matrices (default) or LD matrix
addRequired(p, 'data_path', @isstr);

% directory containing adjacency matrices, defaults to data_path
addParameter(p, 'edgelist_dir', data_directory, @isstr);

% directory for saving (optional, but error is thrown if this is empty and
% nargout == 0)
addParameter(p, 'output_dir', '', @isstr);

% input LD correlation matrix or genotypes? Options: 'genotypes' (default)
% or 'correlation'
addParameter(p, 'data_type', 'genotypes', @isstr);

% filename extension for data files
addParameter(p, 'data_file_extension', '', @isstr);

% among files matching [data_path data_pattern data_file_extension], which one to use
addParameter(p, 'data_file_index', 1, @isscalar);

% expression to match for data files, e.g. '*chr22_*'
addParameter(p, 'data_pattern', '*', @isstr);

% directory containing LDGM SNP list, defaults to data_path
addParameter(p, 'snplist_dir', data_directory, @isstr);

% directory containing LD matrix SNP list, defaults to data_path
addParameter(p, 'LD_matrix_snplist_dir', '', @isstr);

% Whether to match SNPs in LD_matrix_snplist_dir using position column or
% site_ids column
addParameter(p, 'match_snps_on_position', true, @isstr);

% name of LD matrix SNP list column that contains allele frequencies. If
% not specified, a column of the LDGM .snplist will be used instead.
addParameter(p, 'AF_column_name', '', @isstr);

% name of LD matrix SNP list column that contains one of the alleles
% (non-case-sensitive).
addParameter(p, 'A1_column_name', '', @isstr);

% name of LD matrix SNP list column that contains the other allele.
addParameter(p, 'A2_column_name', '', @isstr);

% extra filename field, added at beginning
addParameter(p, 'custom_filename', '', @isstr);

% retain LDGM edges with path distance < threshold
addParameter(p, 'path_distance_threshold', 4, @(x)isscalar(x) && all(x > 0,'all'));

% L1 penalty parameter
addParameter(p, 'l1_penalty', 0.1, @isscalar);

% file path containing population data for each sample in genotype matrix
addParameter(p, 'population_data_file', '', @isstr);

% name of population for which to infer precision matrix, eg 'EUR'
addParameter(p, 'population_name', 'ALL', @isstr);

% meta-analysis weights across superpopulations
addParameter(p, 'meta_weights', [], @isvector);

% file path containing local ancestry matrix. Local ancestry file should
% have same name as .genos file, ending in .localancestry instead. File
% format should be same as .genos file (SNPs x individuals), each entry
% encoding which ancestry group individual ii belongs to at position jj
addParameter(p, 'local_ancestry_dir', '', @isstr);

% number of individuals to downsample_fraction
addParameter(p, 'downsample_fraction', [], @isscalar);

% verbosity for gradient descent function (0,1,2)
addParameter(p, 'verbose', 1, @isscalar);

% minimum MAF filter
addParameter(p, 'minimum_maf', 0.01, @isscalar);

% maximum number of iterations with L1 penalty
addParameter(p, 'penalized_iterations', 5, @isscalar);

% maximum number of iterations without L1 penalty
addParameter(p, 'unpenalized_iterations', 25, @isscalar);

% lasso iterations for inner loop of coordinate descent
addParameter(p, 'lasso_iterations', 10, @isscalar);

% whether to eliminate missing SNPs by calling reduce_weighted_graph to
% patch paths that are lost, or (default) to simply remove those SNPs
addParameter(p, 'patch_paths', false, @isscalar);

% instead of using adjacency matrix, create adjacency matrix with nearest
% neighbors. Number of neighbors will be chosen based on adjacency matrix
% and path_distance_threshold
addParameter(p, 'banded_control_ldgm', 0, @isscalar);

% instead of using adjacency matrix, create adjacency matrix with most
% highly correlated neighbors. . Number of neighbors will be chosen based
% on adjacency matrix and path_distance_threshold
addParameter(p, 'r2_control_ldgm', 0, @isscalar);

% turns p.Results.x into just x
parse(p, data_directory, varargin{:});
fields = fieldnames(p.Results);
for k=1:numel(fields)
    line = sprintf('%s = p.Results.%s;', fields{k}, fields{k});
    eval(line);
end

if isempty(output_dir) && nargout == 0
    error('Please specify output_dir so that the results are saved')
end

% Add dependencies to path
if isfolder('../utility')
    addpath('../utility')
else
    warning('Expected to find a directory with utility functions, but did not. Check directory structure.')
end

if isfolder('../../MATLAB')
    addpath('../../MATLAB')
else
    warning('Expected to find a MATLAB directory, but did not. Check directory structure.')
end

% suffixes for output files
output_suffix = sprintf('.path_distance=%.1f.l1_pen=%.2f.maf=%.2f.%s',...
    path_distance_threshold,l1_penalty,minimum_maf,population_name);

% Expression for what type of data file to look for (genotypes, correlation,
% or edgelist)
if isempty(data_file_extension)
    if strcmp(data_type,'genotypes')
        data_file_extension = '.genos';
    elseif strcmp(data_type,'correlation')
        data_file_extension = '.correlation';
    elseif strcmp(data_type,'edgelist')
        data_file_extension = '.edgelist';
    else
        error('Options for data_type are genotypes, correlation, edgelist')
    end
end

% Threshold to decide if precision matrix entries are zero
smallNumber = 1e-12;

% data file
if isfolder(data_directory)
    data_files = dir([data_directory, data_pattern, data_file_extension]);
    assert(~isempty(data_files));
    filename = data_files(data_file_index).name;
    filename = filename(1:end-length(data_file_extension));
elseif isfile(data_directory)
    data_files = dir(data_directory);
    filename = data_files.name(1:end-length(data_file_extension));
    data_directory = [data_files.folder, '/'];
else
    error('No data files found');
end

% check if output file already exists
if isfile([output_dir,custom_filename,filename,output_suffix,'.stats.txt'])
    error('Output .stats file already exists')
end

% Load genotype matrix and compute correlation matrix
if strcmp(data_type,'genotypes')

    % genotype matrix
    X = load([data_directory, filename, data_file_extension])';
    assert(all(unique(X(:)) == [0 1]'), 'genos file should be binary');

    % which samples to retain
    if ~isempty(population_data_file) && ~strcmp(population_name,'ALL')
        T = readtable(population_data_file);
        superpops = table2cell(T(:,3));
        assert(length(superpops) == size(X,1),...
            'Number of rows in population_data_file should match number of samples in .genos file')
        if ~strcmp(population_name,'META')
            assert(isempty(meta_weights))
            rows = cellfun(@(s)strcmp(s,population_name),superpops);
            assert(~isempty(rows))
            X = X(rows,:);
        end
    end
    [noSamples, noSNPs] = size(X);

    % local ancestry file to handle admixture (or other structure)
    if ~isempty(local_ancestry_dir)
        assert(isempty(population_data_file),...
            'Specify either local ancestry or population labels but not both')
        local_ancestry_file = [local_ancestry_dir, data_file_extension(1:end-6), '.localancestry'];
        local_ancestry = load(local_ancestry_file)';
        noPopns = max(local_ancestry(:));
        AF = mean(X);
        for ii = 1:noSNPs
            for pop = 1:noPopns
                if any(local_ancestry(:,ii) == pop)
                    X(local_ancestry(:,ii) == pop, ii) = ...
                        X(local_ancestry(:,ii) == pop, ii) - ...
                        mean(X(local_ancestry(:,ii) == pop, ii));
                end
            end
        end
    end

    % meta-analysis across populations
    if strcmp(population_name,'META')
        superpop_names = unique(superpops);
        assert(length(superpop_names) == length(meta_weights), ...
            'meta-analysis weights should be a vector of length number of populations')
        for ii=1:length(meta_weights)
            rows = cellfun(@(s)strcmp(s,superpop_names(ii)),superpops);
            AF(ii,:) = mean(X(rows,:));
            X(rows,:) = (X(rows,:) - AF(ii,:)) * meta_weights(ii)/sum(meta_weights) * size(X,1)/sum(rows);
        end
        AF = meta_weights*AF/sum(meta_weights);
    elseif isempty(local_ancestry_dir)
        AF = mean(X);
    end
    MAF = min(AF,1-AF);

    % downsampling (without replacement)
    if ~isempty(downsample_fraction)
        assert(~strcmp(population_name,'META'),'META option incompatitble with downsampling')
        assert(downsample_fraction < 1 & downsample_fraction > 0)
        samples = randsample(sum(rows),floor(sum(rows)*downsample_fraction),false);
        X_out = X(setdiff(1:sum(rows),samples),:);
        R_out = corr(X_out);
        X = X(samples,:);
    end

    % LD correlation matrix
    R = corr(X);

elseif strcmp(data_type,'correlation')
    % LD correlation matrix
    R = load([data_directory,filename, data_file_extension]);

elseif strcmp(data_type,'edgelist')
    % LD correlation matrix
    R = readedgelist([data_directory,filename, data_file_extension]);
end

% SNP list file
if isfolder(snplist_dir)
    snplist_files = dir([snplist_dir,filename,'.snplist']);
    assert(~isempty(snplist_files),"SNP list not found in [snplist_dir,filename,'.snplist']")
    ldgmSnpTable = readtable([snplist_dir,filename,'.snplist'],'FileType','text');
elseif isfile(snplist_dir)
    ldgmSnpTable = readtable(snplist_dir,'FileType','text');
else
    error('SNP list not found')
end

noIndices = numel(unique(ldgmSnpTable.index));
[~,index_representatives] = unique(ldgmSnpTable.index);

% LDGM file
if isfolder(edgelist_dir)
    edgelist_file = dir([edgelist_dir,filename,'.ldgm.edgelist']);
    assert(~isempty(edgelist_file),"Adjacency list file not found in [edgelist_dir,filename,'.edgelist']");
elseif isfile(edgelist_dir)
    edgelist_file = dir(edgelist_dir);
    edgelist_dir = [edgelist_file.folder,'/'];
else
    error('Edge list not found')
end

replace_zeros = eps;
A_weighted = readedgelist([edgelist_dir, edgelist_file.name], noIndices, replace_zeros);

% Merge SNPs between LDGM and correlation matrix

% Allele frequencies
if isempty(AF_column_name)
    assert(any(strcmp(ldgmSnpTable.Properties.VariableNames,population_name)),...
        'No allele frequency column in LDGM SNP list found')
    AF = table2array(ldgmSnpTable(:,population_name));
    MAF = min(AF,1-AF);
end

% If no LD matrix SNP list is specified, require that LD matrix is
% already merged with SNP list
if isempty(LD_matrix_snplist_dir)
    assert(noSNPs == height(ldgmSnpTable),...
        'LD_matrix_snplist_dir not specified, and size of LD matrix does not match size of LDGM')
    ldSnpTable = [];
elseif isfolder(LD_matrix_snplist_dir)
    ldSnpTable = readtable([LD_matrix_snplist_dir,filename,'.snplist'],'FileType','text');
elseif isfile(LD_matrix_snplist_dir)
    ldSnpTable = readtable(LD_matrix_snplist_dir,'FileType','text');
else
    error('LD matrix SNP list not found')
end


% Merge SNPs between the LD SNP list and LDGM SNP list
if ~isempty(ldSnpTable)
    if match_snps_on_position
        [~, idxLDGM, idxR] = intersect(ldgmSnpTable.position, ldSnpTable.position, 'stable');
    else
        [~, idxLDGM, idxR] = intersect(ldgmSnpTable.site_ids, ldSnpTable.site_ids, 'stable');
    end
    if ~isempty(A1_column_name) && ~isempty(A2_column_name)
        phase = mergealleles(...
            table2array(ldSnpTable(idxR,A1_column_name)), ...
            table2array(ldSnpTable(idxR,A2_column_name)),...
            ldgmSnpTable.anc_alleles(idxLDGM),ldgmSnpTable.deriv_alleles(idxLDGM));
        if mean(phase == 0) > 0.5
            error('>50% of SNPs had mismatching alleles between LDGM and correlation matrix edgelists, indicating something is wrong')
        end
    end
    % phase==0 for mismatching alleles
    idxR = idxR(phase~=0); idxLDGM = idxLDGM(phase~=0);
    phase = phase(phase~=0);

    ldSnpTable = ldSnpTable(idxR,:);
    R = phase .* R(idxR,idxR) .* phase';

    if ~isempty(AF_column_name)
        AF = table2array(ldSnpTable(:,AF_column_name));
        MAF = min(AF,1-AF);
    else
        MAF = MAF(idxLDGM);
    end
else
    idxLDGM = 1:length(MAF);
end

% Restrict to common SNPs
commonSNPs = MAF > minimum_maf;
R = R(commonSNPs, commonSNPs);


% Convert merged SNPs to merged indices
[mergedIndices, representatives] = unique(ldgmSnpTable(idxLDGM(commonSNPs),:).index,'stable');
R = R(representatives, representatives);


% Reduce A_weighted to bricks with a matching SNP in the LD matrix
if patch_paths
    % Restrict LDGM to indices corresponding to a common SNP
    commonIndices = unique(ldgmSnpTable.index, 'stable');
    A_weighted = A_weighted(commonIndices + 1, commonIndices + 1);

    idx = lift(mergedIndices + 1, commonIndices + 1);
    missingIdx = setdiff(1:length(A_weighted),idx);
    if length(missingIdx) > length(idx)/3
        warning('More than 1/4 of SNPs in LDGM are missing, possibly leading to a low quality solution');
    end
    A_weighted = reduce_weighted_graph(A_weighted,...
        missingIdx, 0);
else
    A_weighted = A_weighted(mergedIndices + 1, mergedIndices + 1);
end

% Number of SNPs remaining
noSNPs = length(R);

% Initial graphical model
A = spfun(@(x)x < path_distance_threshold, A_weighted + eps * speye(noSNPs));

% banded control option
if banded_control_ldgm
    bandSize = ceil(nnz(A) / (2*noSNPs));
    A = spdiags(true(noSNPs,2*bandSize+1),-bandSize:bandSize,noSNPs,noSNPs);
    output_suffix = [output_suffix, '.banded_control'];
end

% r2 control option
if r2_control_ldgm
    assert(~banded_control_ldgm, ...
        'banded_control_ldgm and r2_control_ldgm are incompatible options')
    assert(~strcmp(data_type,'edgelist'), ...
        'data_type edgelist and r2_control_ldgm are incompatible options')
    r2 = R(:).^2;
    r2_threshold = quantile(r2, 1 - mean(A(:)));
    A = R.^2 >= r2_threshold;
    output_suffix = [output_suffix, '.r2_control'];
end

% Initial average degree of LDGM, before L1 penalized step
initialDegree = nnz(A) / noSNPs;

% estimate precision matrix with L1 penalty
tic;
precisionEstimatePenalized = precisionCoordinateDescent(A,...
    R, penalized_iterations, l1_penalty, lasso_iterations);
time_penalized=toc;

% Updated graphical model
A_l1 = abs(precisionEstimatePenalized) > smallNumber;
avgDegree = nnz(A_l1) / length(A_l1);

% initialize at regularized estimate
precisionEstimatePenalized(~A_l1)=0;
[~,a] = chol(precisionEstimatePenalized);
while a > 0 % precisionEstimatePenalized is not positive definite
    precisionEstimatePenalized = precisionEstimatePenalized + 0.01 * speye(noSNPs);
    [~,a] = chol(precisionEstimatePenalized);
end

% Estimate precision matrix at edges selected with L1 penalty, dropping the
% penalty term
tic
precisionEstimate = precisionCoordinateDescent(A_l1,...
    R, unpenalized_iterations, 0, 0, precisionEstimatePenalized);
time_unpenalized=toc;

% Error metrics
Rr = inv(precisionEstimate);
mse = mean((Rr - R).^2, 'all');

RP = R * precisionEstimate;
PR = RP';
II = eye(noSNPs);
alt_mse = mean((II - RP) .* (II - PR), 'all');

% stats table
T = table('Size',[1 0]);
T.distance_threshold = path_distance_threshold;
T.l1_pen = l1_penalty;
T.degree = avgDegree;
T.size = noSNPs;
T.initial_degree = initialDegree;
T.mse = mse;
T.alt_mse = alt_mse;
T.runtime = time_unpenalized + time_penalized;


if ~isempty(downsample_fraction)
    T.mse_P_out = mean((Rr(~A) - R_out(~A)).^2);
    T.mse_R_out = mean((R(~A) - R_out(~A)).^2);
end


% saving
if ~isempty(output_dir)
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    % assign entries of precisionEstimate to rows/cols corresponding to the
    % rows/cols of the LDGM
    P = zeros(max(ldgmSnpTable.index) + 1);
    P(mergedIndices + 1, mergedIndices + 1) = precisionEstimate;

    writetable(T, [output_dir,custom_filename,filename,output_suffix,'.stats.txt']);
    writeedgelist([output_dir,custom_filename,filename,output_suffix,'.edgelist'], P);
end


