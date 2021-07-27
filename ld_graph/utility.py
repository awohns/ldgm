"""
Utility functions
"""
import collections

import numpy as np
from tqdm import tqdm


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


def get_brick_descendants(ts):
    X = np.zeros((ts.num_samples, ts.num_edges))
    node_edge_dict = {}
    for tree2, (_, edges_out, edges_in) in tqdm(zip(ts.trees(), ts.edge_diffs())):
        for edge in edges_out:
            node_edge_dict.pop(edge.child)
        for edge in edges_in:
            node_edge_dict[edge.child] = edge.id
        for key, val in node_edge_dict.items():
            for sample in tree2.samples(key):
                X[sample, val] = 1
    return X


def add_dummy_bricks(bts, epsilon=1e-6):
    # Check that the first (num_samples) nodes are all samples
    for i in range(bts.num_samples):
        assert i in bts.samples()

    node_mapping = {}
    tables = bts.dump_tables()

    tables.nodes.clear()
    for sample in bts.samples():
        tables.nodes.add_row(flags=1, time=bts.node(sample).time)

    # Add all the nodes in
    for node in bts.nodes():
        if node.id not in bts.samples():
            node_mapping[node.id] = tables.nodes.add_row(flags=0, time=node.time)
        else:
            # Add a dummy for samples
            node_mapping[node.id] = tables.nodes.add_row(
                flags=0, time=node.time + epsilon
            )
    tables.edges.clear()
    sequence_length = bts.get_sequence_length()
    # Then we add bricks in
    for dummy in bts.samples():
        tables.edges.add_row(
            left=0, right=sequence_length, parent=node_mapping[dummy], child=dummy
        )
    for edge in bts.edges():
        tables.edges.add_row(
            left=edge.left,
            right=edge.right,
            parent=node_mapping[edge.parent],
            child=node_mapping[edge.child],
        )
    # Fix the mutation nodes
    tables.mutations.clear()
    for mut in bts.mutations():
        tables.mutations.add_row(
            site=mut.site,
            node=node_mapping[mut.node],
            derived_state=mut.derived_state,
            time=mut.time,
        )
    # Then make a new brick tree sequence
    tables.sort()
    return tables.tree_sequence()
