function [matrices, snplists] = loadLDGMs(filepath,popn_names)
%loadLDGMs reads LDGMs or LDGM precision matrices from
%specified file or directory, together with corresponding SNP lists
% It expects to find .edgelist files named [filepath, popn_name{j}, '.edgelist']
% and .snplist files named [filepath, '.snplist']. If you
% are loading data for multiple populations, specify the population names
% as a cell array of strings, e.g. {'AFR','EUR'}; otherwise, specify it
% just as a string.

if isstr(popn_names)
    popn_names = {popn_names};
end
snplist_files = dir([filepath,'*','.snplist']);
noFiles = length(snplist_files);

% Check that edgelist files exist
for ii = 1:length(popn_names)
    for jj = 1:noFiles
        assert(isfile([snplist_files(jj).folder, '/', ...
            extractBefore(snplist_files(jj).name,'snplist'), ...
            popn_names{ii}, '.edgelist']), 'edglist file not found: please put edgelist and snplist files in the same directory, and ensure that there is an edgelist file for every snplist-population pair')
    end
end

snplists = cell(noFiles,1);
matrices = cell(noFiles,length(popn_names));

% Read snplists
for jj = 1:length(snplist_files)
    snplists{jj} = readtable([snplist_files(jj).folder,'/',snplist_files(jj).name],'FileType','text');
end

% Read edgelists
for ii = 1:length(popn_names)
    for jj = 1:length(snplist_files)
        matrices{jj,ii} = readedgelist([snplist_files(jj).folder, '/', ...
            extractBefore(snplist_files(jj).name,'snplist'), ...
            popn_names{ii}, '.edgelist'], max(snplists{jj}.index)+1);
    end
end


end

