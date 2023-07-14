function [matrices, snplists, AF] = loadLDGMs(filepath,varargin)
%loadLDGMs reads LDGMs or LDGM precision matrices from
%specified file or directory, together with corresponding SNP lists
% 
% Input arguments:
% filepath: where to look for .edgelist files to load
% 
% Optional input arguments (name-value pairs):
%  popnNames: names of populations to be included, or arbitrary strings that must be
%  included in the filenames
%  whichBlocks: indices of which LD blocks to include, out of all of those that are found
%  snsplistpath: where to look for SNP lists, if different from filepath
%  normalizePrecision: whether to normalize precision matrices such that
%  the diagonal of their inverse is exactly 1

p=inputParser;

% directory containing genotype matrices (default) or LD matrix
addRequired(p, 'filepath', @isstr);

% names of populations to be included, or arbitrary strings that must be
% included in the filenames
addParameter(p, 'popnNames', '', @(x)isstr(x) || iscell(x));

% indices of which LD blocks to include, out of all of those that are found
addParameter(p, 'whichBlocks', 0, @isnumeric);

% directory containing SNP lists
addParameter(p, 'snplistpath', filepath, @isstr);

% whether to normalize LDGM precision matrices such that the diagonal of
% their inverse is exactly 1
addParameter(p, 'normalizePrecision', false, @isscalar);


% whether to add column containing chromosome to each SNP list
addParameter(p, 'addChromosomeColumn', false, @isscalar);

% turns p.Results.x into just x
parse(p, filepath, varargin{:});
fields = fieldnames(p.Results);
for k=1:numel(fields)
    line = sprintf('%s = p.Results.%s;', fields{k}, fields{k});
    eval(line);
end

if isstr(popnNames)
    popnNames = {popnNames};
end
if ~exist('snplistpath')
    snplistpath = filepath;
end
snplist_files = dir([snplistpath,'*','.snplist']);

if all(whichBlocks > 0)
    snplist_files = snplist_files(whichBlocks);
end
noFiles = length(snplist_files);
if noFiles == 0
    warning('No files were found')
end

snplists = cell(noFiles,1);
matrices = cell(noFiles,length(popnNames));

% Read snplists
for jj = 1:noFiles
    snplists{jj} = readtable([snplist_files(jj).folder,'/',snplist_files(jj).name],'FileType','text');
    if addChromosomeColumn
        before = extractBefore(snplist_files(jj).name, 'chr');
        chr = sscanf(snplist_files(jj).name, [before, 'chr%d']);
        snplists{jj}.chr(:) = chr;
    end
end

% Read edgelists
edgelist_files = dir([filepath, '*.edgelist']);
filedir = edgelist_files(1).folder;
filenames = {edgelist_files.name};
for ii = 1:length(popnNames)
    for jj = 1:noFiles
        pattern = extractBefore(snplist_files(jj).name,'snplist');
        file = contains(filenames, pattern) & contains(filenames,popnNames{ii});
        if sum(file) == 1
            matrices{jj,ii} = readedgelist([filedir, '/', ...
                filenames{file}], max(snplists{jj}.index)+1);
            if normalizePrecision
                incl = any(matrices{jj,ii});
                D = ones(length(matrices{jj,ii}),1);
                D(incl) = sqrt(diag(sparseinv(matrices{jj,ii}(incl,incl))));
                matrices{jj,ii} = D .* matrices{jj,ii} .* D';
            end
        elseif sum(file) == 0
            matrices{jj,ii} = [];
            warning('No matching edgelist found for at least one LD block')
        else
            error('Multiple matching edgelists found for at least one LD block')
        end
    end
end

% Get allele frequencies
if nargout >= 3
    [noBlocks,noPopns] = size(matrices);
    AF = cell(noBlocks,noPopns);
    for ii = 1:noBlocks
        % call unique() on the SNP list indices with two output arguments in
        % order to pick a representative SNP for each index
        [~,representatives] = unique(snplists{ii}.index,'stable');

        % snplists table can be sliced using the names of each column as a cell
        % array of strings
        AF(ii,:) = mat2cell(table2array(snplists{ii}(representatives,popnNames)),...
            length(representatives),ones(1,noPopns));
    end
end
end



