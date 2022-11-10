function R = writeCorrelationEdgelist(data_directory, varargin)
% writeCorrelationEdgelist inputs a genotype matrix and an LDGM, and it
% outputs the correlations between pairs of SNPs with edges in the LDGM, as
% a .edgelist file.
%
% Input data (located at data_directory) is an individials-by-SNPs binary
%   matrix with filename extension
%   .genos, and columns that correspond to the rows of .snplist file (see
%   below).
%
%
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
%   -local_ancestry_dir: not yet supported
%
% FILENAME HANDLING OPTIONS:
%   -edgelist_dir: directory containing LDGM .edgelist files. Defaults to
%   data_directory.
%   snplist_dir: directory containing .snplist files. Defaults to
%   data_directory.
%   -output_dir: directory in which to save the output files
%   -custom_filename: if specified, this is prepended to the name of the
%   output files within output_dir
%   -data_pattern: if specified, only data files matching this regular
%   expression will be found
%   -data_file_index: out of the data files that are found, which one to
%   compute a precision matrix for
%
% METHOD PARAMETERS: (all of these can probably be left at default values)
%   -minimum_maf: minimum allele frequency, as calculated either using the
%   LDGM .snplist file (default) or the LD matrix .snplist file (if
%   AF_column_name is specified)
%   -path_distance_threshold: LDGM entries with distance below this
%   threshold will be retained. Default value is 4.
%
% OUTPUT FILES:
%   -'*.correlationMatrix.edgelist' is an edge list with each row
%   corresponding to an entry of the LDGM. This file is
%   zero-indexed, and the rows/columns of the precision matrix correspond
%   to the 'index' column of the LDGM .snplist file.
%
% Instead of saving the output, it is also an option to return the
% correlation matrix R directly.

p=inputParser;

% directory containing genotype matrices (default) or LD matrix
addRequired(p, 'data_path', @isstr);

% directory containing adjacency matrices, defaults to data_path
addParameter(p, 'edgelist_dir', data_directory, @isstr);

% directory for saving (optional, but error is thrown if this is empty and
% nargout == 0)
addParameter(p, 'output_dir', '', @isstr);

% filename extension for data files
addParameter(p, 'data_file_extension', '', @isstr);

% among files matching [data_path data_pattern data_file_extension], which one to use
addParameter(p, 'data_file_index', 1, @isscalar);

% expression to match for data files, e.g. '*chr22_*'
addParameter(p, 'data_pattern', '*', @isstr);

% directory containing LDGM SNP list, defaults to data_path
addParameter(p, 'snplist_dir', data_directory, @isstr);

% extra filename field, added at beginning
addParameter(p, 'custom_filename', '', @isstr);

% retain LDGM edges with path distance < threshold
addParameter(p, 'path_distance_threshold', 4, @(x)isscalar(x) && all(x > 0,'all'));

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

% minimum MAF filter
addParameter(p, 'minimum_maf', 0.01, @isscalar);


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

% suffixes for output files
output_suffix = sprintf('.path_distance=%.1f.maf=%.2f.%s',...
    path_distance_threshold,minimum_maf,population_name);

% Expression for what type of data file to look for (genotypes, correlation,
% or edgelist)
if isempty(data_file_extension)
    data_file_extension = '.genos';
end

% data file
data_files = dir([data_directory, data_pattern, data_file_extension]);
assert(~isempty(data_files));
filename = data_files(data_file_index).name;
filename = filename(1:end-length(data_file_extension));

% check if output file already exists
if isfile([output_dir,custom_filename,filename,output_suffix,'.correlationMatrix.edgelist'])
    error('Output edgelist file already exists')
end

% Load genotype matrix and compute correlation matrix

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

% SNP list file
snplist_files = dir([snplist_dir,filename,'.snplist']);
if isempty(snplist_files)
    error("SNP list not found in [snplist_dir,filename,'.snplist']")
else
    ldgmSnpTable = readtable([snplist_dir,filename,'.snplist'],'FileType','text');
    [~,index_representatives] = unique(ldgmSnpTable.index);
    noIndices = length(index_representatives);
    X = X(:,index_representatives);
    MAF = MAF(index_representatives);
end

% LDGM file
edgelist_file = dir([edgelist_dir,filename,'.ldgm.edgelist']);
assert(~isempty(edgelist_file),"Adjacency list file not found in [edgelist_dir,filename,'.edgelist']");
replace_zeros = eps;
A_weighted = readedgelist([edgelist_dir, edgelist_file.name], noIndices, replace_zeros);

% LD correlation matrix
commonSNPs = MAF > minimum_maf;
R = zeros(noIndices);
R(commonSNPs,commonSNPs) = corr(X(:, commonSNPs));

% Graphical model
A = spfun(@(x)x < path_distance_threshold, A_weighted + eps * speye(noIndices));
[ii,jj] = find(A);

% Sparse correlation matrix with entries corresponding to A
R = sparse(ii,jj,R(A),noIndices,noIndices);

% saving
if ~isempty(output_dir)
    mkdir(output_dir);

    writeedgelist([output_dir,custom_filename,filename,output_suffix,'.correlationMatrix.edgelist'],R);
end


