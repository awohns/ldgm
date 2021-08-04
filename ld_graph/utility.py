"""
Utility functions
"""
import collections

import numpy as np


def get_mut_edges(ts):
    muts_to_brick = {}
    node_edge_dict = {}
    for tree, (_, edges_out, edges_in) in zip(ts.trees(), ts.edge_diffs()):
        for edge in edges_out:
            node_edge_dict.pop(edge.child)
        for edge in edges_in:
            node_edge_dict[edge.child] = edge.id
        for site in tree.sites():
            for mut in site.mutations:
                node = mut.node
                if node in node_edge_dict:
                    muts_to_brick[mut.id] = node_edge_dict[node]

    bricks_to_muts = collections.defaultdict(list)

    for mut, brick in muts_to_brick.items():
        bricks_to_muts[brick].append(mut)
    return bricks_to_muts


def get_brick_frequencies(ts):
    freqs = np.zeros(ts.num_edges)
    for tree, (_, _, edges_in) in zip(ts.trees(), ts.edge_diffs()):
        for edge in edges_in:
            # assert that we've never seen the brick before
            assert freqs[edge.id] == 0
            freqs[edge.id] = tree.num_samples(edge.child) / tree.num_samples()
    return freqs
