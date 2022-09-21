function [matrices, snplists] = loadLDGMs(filepath,popn_names,whichBlocks)
%loadLDGMs reads LDGMs or LDGM precision matrices from
%specified file or directory, together with corresponding SNP lists
% It expects to find .edgelist files named [filepath, '*', popn_names{j}, '.edgelist']
% and .snplist files named [filepath, '*', '.snplist']. If you
% are loading data for multiple populations, specify the population names
% as a cell array of strings, e.g. {'AFR','EUR'}; otherwise, specify it
% just as a string.
% Specify whichBlocks to load only one or a few LD blocks at a time; this
% should be an array of indices

if isstr(popn_names)
    popn_names = {popn_names};
end
snplist_files = dir([filepath,'*','.snplist']);
if nargin >= 3
    snplist_files = snplist_files(whichBlocks);
end
noFiles = length(snplist_files);

snplists = cell(noFiles,1);
matrices = cell(noFiles,length(popn_names));

% Read snplists
for jj = 1:noFiles
    snplists{jj} = readtable([snplist_files(jj).folder,'/',snplist_files(jj).name],'FileType','text');
end

% Read edgelists
for ii = 1:length(popn_names)
    for jj = 1:noFiles
        if isfile([snplist_files(jj).folder, '/', ...
                extractBefore(snplist_files(jj).name,'snplist'), ...
                popn_names{ii}, '.edgelist'])
            matrices{jj,ii} = readedgelist([snplist_files(jj).folder, '/', ...
                extractBefore(snplist_files(jj).name,'snplist'), ...
                popn_names{ii}, '.edgelist'], max(snplists{jj}.index)+1);
        else
            matrices{jj,ii} = [];
        end
    end
end


end

