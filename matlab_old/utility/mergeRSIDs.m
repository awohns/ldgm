function [merged_rsid_cells, is_inlist_cells, idx_cells, phase] = mergeRSIDs(rsid_cells,rsid_list,a1_cells,a1_list)
%mergeRSIDs merges a cell array containing cell arrays of RSIDs with a
%single cell array of RSIDs
% 
% Input arguments:
% rsid_cells: cell array where each cell contains a (cell) array
% rsid_list: (cell) array to be merged with each element of rsid_cells
% a1_cells: allele 1, a cell array where each cell contains a (cell) array
% with strings saying which allele is A1, corresonding to rsid_cells
% a1_list: a cell array of strings saying which allele is A1, corresponding
% to rsid_list
% 
% Ouput arguments:
% merged_rsid_cells: the intersection between each element of rsid_cells
% and rsid_list, in same order as rsid_cells started
% is_inlist_cells: cell array of boolean vectors, such that
% merged_rsid_cells{j} == rsid_cells{j}(is_inlist_cells{j})
% idx_cells: cell array of indices, such that 
% merged_rsid_cells{j} == rsid_list(idx_cells{j})
% phase: cell array of vectors whose entries are +-1, where +1 means that
% a1_cells and a1_list agree, and -1 means that they disagree

if iscolumn(rsid_cells)
    rsid_cells = rsid_cells';
end
all_ldgm_rsids = vertcat(rsid_cells{:});

[~, ldgm_idx, sumstats_idx] = intersect(all_ldgm_rsids, rsid_list, 'stable');

blocksizes = cellfun(@length,rsid_cells);
cumulative_blocksizes = [0,cumsum(blocksizes)];

blocks = arrayfun(@(i,j){i+1:j}, cumulative_blocksizes(1:end-1), cumulative_blocksizes(2:end));

is_inlist_cells = false(cumulative_blocksizes(end),1);
is_inlist_cells(ldgm_idx) = true;
is_inlist_cells = cellfun(@(i){is_inlist_cells(i)}, blocks);

merged_rsid_cells = cellfun(@(rs,incl){rs(incl)}, rsid_cells, is_inlist_cells);

blocksizes = cellfun(@sum,is_inlist_cells);
cumulative_blocksizes = [0,cumsum(blocksizes)];
blocks = arrayfun(@(i,j){i+1:j}, cumulative_blocksizes(1:end-1), cumulative_blocksizes(2:end));

idx_cells = cellfun(@(i){sumstats_idx(i)}, blocks);

% phasing alleles
if nargout > 3
    assert(nargin == 4)
    if iscolumn(a1_cells)
        a1_cells = a1_cells';
    end
    a1_cells = cellfun(@(C,incl){C(incl)},a1_cells,is_inlist_cells);
    a1_list_cells = cellfun(@(incl){a1_list(incl)},idx_cells);
    phase = cellfun(@(C1,C2){-(-1).^cellfun(@strcmp,C1,C2)}, a1_cells, a1_list_cells);
end


end

