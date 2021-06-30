"""
Utility functions
"""
import collections


def get_mut_edges(ts):
    muts_to_brick = {}
    node_edge_dict = {}
    for tree, (interval, edges_out, edges_in) in zip(ts.trees(), ts.edge_diffs()):
        for edge in edges_out:
            node_edge_dict.pop(edge.child)
        for edge in edges_in:
            node_edge_dict[edge.child] = edge.id
        for site in tree.sites():
            for mut in site.mutations:
                node = mut.node
                muts_to_brick[mut.id] = node_edge_dict[node]

    bricks_to_muts = collections.defaultdict(list)

    for mut, brick in muts_to_brick.items():
        bricks_to_muts[brick].append(mut)
    return bricks_to_muts
