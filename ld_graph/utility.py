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


def interval_while_leaf(ts):
    """
    Returns a dictionary, with keys equal to leaves and values a list of
    lists of intervals where the node is a leaf
    """
    leaf_spans = collections.defaultdict(list)
    for tree in ts.trees():
        for leaf in tree.leaves():
            if tree.parent(leaf) != -1:
                leaf_spans[leaf].append(tree.interval)

    uninterrupted_leaf_spans = collections.defaultdict(list)
    for val, interval in leaf_spans.items():
        cur_interval = interval[0]
        prev_right = cur_interval.right
        uninterrupted_leaf_spans[val].append(cur_interval.left)
        for cur_interval in interval[1:]:
            if cur_interval.left != prev_right:
                uninterrupted_leaf_spans[val].append(prev_right)
                uninterrupted_leaf_spans[val].append(cur_interval.left)
            prev_right = cur_interval.right
        uninterrupted_leaf_spans[val].append(cur_interval.right)
    for dummy, intervals in uninterrupted_leaf_spans.items():
        leaf_list = []
        for idx, _ in enumerate(intervals[::2]):
            leaf_list.append((intervals[idx * 2], intervals[idx * 2 + 1]))
        uninterrupted_leaf_spans[dummy] = leaf_list

    return uninterrupted_leaf_spans


def add_dummy_bricks(bts, mode="samples", epsilon="adaptive"):
    """
    Add a dummy bricks to the tree sequence, allowing each sample or leaf node
    (depending on mode) to be the *parent* node of a brick.
    Dummy nodes are *not* marked as samples.
    Convert from dummy node ids to previous ids by subtracting the number of
    leaves or samples
    (depending on mode used).
    """
    # Check that the first (num_samples) nodes are all samples
    # for i in range(bts.num_samples):
    #    assert i in bts.samples()
    if epsilon == "adaptive":
        time = bts.tables.nodes.time
        gaps = time[bts.tables.edges.parent] - time[bts.tables.edges.child]
        epsilon = np.min(np.unique(gaps)) * 0.1

    node_mapping = {}
    tables = bts.dump_tables()

    tables.nodes.clear()
    if mode == "samples":
        targets = bts.samples()
    elif mode == "leaves":
        targets = set()
        for tree in bts.trees():
            for leaf in tree.leaves():
                if tree.parent(leaf) != -1:
                    targets.add(leaf)

    # Add all the nodes in
    for node in bts.nodes():
        if node.id not in targets:
            tables.nodes.add_row(flags=0, time=node.time)
        else:
            # Add target
            tables.nodes.add_row(flags=node.flags, time=node.time + epsilon)
    # adding dummy nodes
    for target in targets:
        node_mapping[target] = tables.nodes.add_row(flags=0, time=bts.node(target).time)
    tables.edges.clear()
    if mode == "leaves":
        # Then we add bricks in
        leaf_spans = interval_while_leaf(bts)
        for dummy, intervals in leaf_spans.items():
            for interval in intervals:
                tables.edges.add_row(
                    left=interval[0],
                    right=interval[1],
                    parent=dummy,
                    child=node_mapping[dummy],
                )
    elif mode == "samples":
        sequence_length = bts.get_sequence_length()
        for target in targets:
            tables.edges.add_row(
                left=0, right=sequence_length, parent=target, child=node_mapping[target]
            )
    for edge in bts.edges():
        tables.edges.add_row(
            left=edge.left,
            right=edge.right,
            parent=edge.parent,
            child=edge.child,
        )
    # Fix the mutation nodes
    tables.mutations.clear()
    for mut in bts.mutations():
        tables.mutations.add_row(
            site=mut.site,
            node=mut.node,
            derived_state=mut.derived_state,
            time=mut.time,
        )
    # Then make a new brick tree sequence
    tables.sort()
    return tables.tree_sequence()


def check_bricked(ts):
    """
    Checks that the input tree sequence is "bricked" by tree, node, or leaf.
    Returns True if bricked by the given mode, False if not.
    """
    return "TODO"
