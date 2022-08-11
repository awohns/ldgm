function [whichIndices, mergedSumstats, whichSNPs, sumstats_SNPs_in_snplists] = mergesnplists(snplists,sumstats,P)
%mergesnplists merges (1) a cell array of .snplist tables, one for each LD block,
%with (2) a table of summary statistics containing RSIDs.
% 
% Input arguments:
% snplists: cell array where each cell contains a snplist table.
% 
% sumstats: summary statistics table with mandatory column name
% (non-case-sensitive) SNP or RSID, to be used for merging.
%   Optionally, sumstats can also contain columns named A1/A2
%   or anc_alleles/deriv_alleles. If these columns are
%   specified, they will be merged with the anc_alleles/deriv_alleles columns
%   of snplists.
%
% P (optional): cell array of LDGM precision matrices. Any indices whose
%  corresponding rows/columns in P are empty will be ingored.
%  Alternatively, specify a logical row vector.
% 
% Ouput arguments:
% whichIndices: cell array of indices, a subset of the 'index' column of each
% snplist. Indexes which row/column of LDGMs have a corresponding SNP in
% the sumstats table.
% 
% mergedSumstats: cell array of sumstats tables, each containing a subset
% of the rows of the original sumstats table, that match the indices of
% each LD block. Note: when multiple SNPs correspond to the same row/column
% of the LDGM, the first-listed one will arbitrarily be chosen as a
% representative. 
% 
% If alleles are specified in the sumstats table, the mergedSumstats tables
% will have an extra column appended called 'phase', which indicates
% whether the alleles match (+1) or anti-match (-1). If they mismatch
% (i.e., flipping their labels doesn't make them match), then the SNP will
% be discarded.
% 
% Note: sometimes it makes this function run much faster if you go
% chromosome-by-chromosome; this is not done automatically.

assert(iscell(snplists), 'Please specify snplists as a cell array of tables')
assert(all(cellfun(@istable,snplists)), 'Please specify snplists as a cell array of tables')
assert(istable(sumstats), 'Please specify sumstats as a table')

sumstats_colnames = sumstats.Properties.VariableNames;
snpcolumn = strcmpi(sumstats_colnames,'SNP');
if sum(snpcolumn) ==0
    snpcolumn = strcmpi(sumstats_colnames,'RSID');
end
assert(sum(snpcolumn) == 1, ...
    'Please specify exactly one column in sumstats table with name SNP or RSID')

concatenated_snplists = vertcat(snplists{:});

[~, ldgm_idx, sumstats_idx] = intersect(concatenated_snplists.site_ids,...
    table2cell(sumstats(:,snpcolumn)), 'stable');

% Subset sumstats to matching SNPs
sumstats = sumstats(sumstats_idx,:);

% Indices to recover LD blocks from concatenated list
blocksizes = cellfun(@height,snplists);
cumulative_blocksizes = [0;cumsum(blocksizes)];
blocks = arrayfun(@(i,j){i+1:j}, cumulative_blocksizes(1:end-1), cumulative_blocksizes(2:end));

% Which SNPs in the snplists have a corresponding entry in the
% sumstats file
whichSNPs = false(cumulative_blocksizes(end),1);
whichSNPs(ldgm_idx) = true;
whichSNPs = cellfun(@(i){whichSNPs(i)}, blocks);

% Indices to recover merged LD blocks from concatenated list
blocksizes = cellfun(@sum,whichSNPs);
cumulative_blocksizes = [0;cumsum(blocksizes)];
blocks = arrayfun(@(i,j){i+1:j}, cumulative_blocksizes(1:end-1), cumulative_blocksizes(2:end));

% Sumstats tables merged with each snplist
mergedSumstats = cellfun(@(ii){sumstats(ii,:)}, blocks);

% Convert from SNPs to LDGM row/col indices (which can have multiple SNPs)
noBlocks = numel(snplists);
whichIndices = cell(size(snplists));
representatives = cell(size(snplists));
for ii = 1:noBlocks
    [whichIndices{ii}, representatives{ii}] = ...
        unique(snplists{ii}.index(whichSNPs{ii}) + 1);
    
    % Get rid of indices whos corresponding columns of P are empty
    if nargin > 2
        include_indices = find(any(P{ii}));
        [whichIndices{ii}, include_indices] = intersect(whichIndices{ii},include_indices);
        representatives{ii} = representatives{ii}(include_indices);
    end
    mergedSumstats{ii} = mergedSumstats{ii}(representatives{ii},:);
end

% which sumstats SNPs had a matching SNP in each LD block
sumstats_SNPs_in_snplists = cellfun(@(i){sumstats_idx(i)}, blocks);

% phasing alleles
a1column = strcmpi(sumstats_colnames,'A1');
if sum(a1column) == 0
    a1column = strcmpi(sumstats_colnames,'anc_alleles');
end
a2column = strcmpi(sumstats_colnames,'A2');
if sum(a2column) == 0
    a2column = strcmpi(sumstats_colnames,'deriv_alleles');
end
if sum(a1column)==1 && sum(a2column) == 1
    for ii = 1:noBlocks
        
        idx = find(whichSNPs{ii});
        idx = idx(representatives{ii});
        
        % +1 for matching alleles, -1 for anti-matching, 0 for mismatching
        phase = mergealleles(table2cell(mergedSumstats{ii}(:,a1column)), ...
            table2cell(mergedSumstats{ii}(:,a2column)), ...
            snplists{ii}.anc_alleles(idx),...
            snplists{ii}.deriv_alleles(idx));
                
        % assign phase field to merged sumstats
        mergedSumstats{ii}.phase = phase;
        
        % assign Z score field, phased to derived allele
        if any(strcmpi(sumstats_colnames,'Z'))
            z_col = strcmpi(sumstats_colnames,'Z');
            mergedSumstats{ii}.Z_deriv_allele = phase .* ...
                table2array(mergedSumstats{ii}(:,z_col));
        elseif any(strcmpi(sumstats_colnames,'beta'))
            beta_col = strcmpi(sumstats_colnames,'beta');
            se_col = strcmpi(sumstats_colnames,'SE');
            assert(any(se_col),...
                'Sumstats should specify either a z score (Z) or an effect size (beta) and standard error (se)')
            mergedSumstats{ii}.Z_deriv_allele = phase .* ...
                table2array(mergedSumstats{ii}(:,beta_col)) ./...
                table2array(mergedSumstats{ii}(:,se_col));
        else
            error('Sumstats should specify either a z score (Z) or an effect size (beta) and standard error (se)')
        end
        
        % get rid of mismatched alleles
        if mean(phase==0) > 0.5
            warning('In block %d, more than half of putatively matching SNPs had mismatched alleles (perhaps due to strandedness?)',ii)
        end
        mergedSumstats{ii} = mergedSumstats{ii}(phase~=0, :);
        whichIndices{ii} = whichIndices{ii}(phase~=0, :);
    end
else
    warning('Did not find allele information in sumstats file. Please ensure that alleles are matched between summary statistics and the LDGM')
end


end

