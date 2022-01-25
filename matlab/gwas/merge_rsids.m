function [merged_rsid_cells, is_inlist_cells, idx_cells] = merge_rsids(rsid_cells,rsid_list)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

all_ldgm_rsids = vertcat(rsid_cells{:});

[~, ldgm_idx, sumstats_idx] = intersect(all_ldgm_rsids, rsid_list, 'stable');

blocksizes = cellfun(@length,rsid_cells);
cumulative_blocksizes = [0;cumsum(blocksizes)];

blocks = arrayfun(@(i,j){i+1:j}, cumulative_blocksizes(1:end-1), cumulative_blocksizes(2:end));

is_inlist_cells = false(cumulative_blocksizes(end),1);
is_inlist_cells(ldgm_idx) = true;
is_inlist_cells = cellfun(@(i){is_inlist_cells(i)}, blocks);

merged_rsid_cells = cellfun(@(rs,incl){rs(incl)}, rsid_cells, is_inlist_cells);

blocksizes = cellfun(@sum,is_inlist_cells);
cumulative_blocksizes = [0;cumsum(blocksizes)];
blocks = arrayfun(@(i,j){i+1:j}, cumulative_blocksizes(1:end-1), cumulative_blocksizes(2:end));

idx_cells = cellfun(@(i){sumstats_idx(i)}, blocks);



end

