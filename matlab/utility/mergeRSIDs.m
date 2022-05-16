function [merged_rsid_cells, is_inlist_cells, idx_cells] = mergeRSIDs(rsid_cells,rsid_list)
%mergeRSIDs merges a cell array containing cell arrays of RSIDs with a
%single cell array of RSIDs
% 
% Input arguments:
% rsid_cells: cell array where each cell contains a (cell) array
% rsid_list: (cell) array to be merged with each element of rsid_cells
% 
% Ouput arguments:
% merged_rsid_cells: the intersection between each element of rsid_cells
% and rsid_list, in same order as rsid_cells started
% is_inlist_cells: cell array of boolean vectors, such that
% merged_rsid_cells{j} == rsid_cells{j}(is_inlist_cells{j})
% idx_cells: cell array of indices, such that 
% merged_rsid_cells{j} == rsid_list(idx_cells{j})

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



end

