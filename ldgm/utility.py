"""
Utility functions
"""
import collections
import itertools
import json
import math

import numpy as np
import tskit

from . import provenance


def softmin(val1, val2):
    return -2 * np.log(math.e ** (-0.5 * val1) + math.e ** (-0.5 * val2))


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
        if bts.num_mutations > 0:
            if not np.isnan(np.min(bts.tables.mutations.time)):
                epsilon = np.min(
                    [np.min(np.unique(gaps)) * 0.1, np.min(bts.tables.mutations.time)]
                )
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
            tables.nodes.add_row(flags=0, time=node.time + epsilon)
    # adding dummy nodes
    for target in targets:
        node_mapping[target] = tables.nodes.add_row(
            flags=bts.node(target).flags, time=bts.node(target).time
        )
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


def remove_node(g, node, path_threshold, use_softmin=False):
    if g.is_directed():
        sources = [source for source, _ in g.in_edges(node)]
        targets = [target for _, target in g.out_edges(node)]
    else:
        raise ValueError

    new_edges = itertools.product(sources, targets)
    new_edges_no_self = []
    for source, target in new_edges:
        if source != target:
            combined_weight = (
                g.get_edge_data(source, node)["weight"]
                + g.get_edge_data(node, target)["weight"]
            )
            if g.has_edge(source, target):
                if use_softmin:
                    combined_weight = softmin(
                        g.get_edge_data(source, target)["weight"], combined_weight
                    )
                else:
                    combined_weight = np.minimum(
                        g.get_edge_data(source, target)["weight"], combined_weight
                    )
            if combined_weight <= path_threshold:
                new_edges_no_self.append((source, target, combined_weight))
    g.add_weighted_edges_from(new_edges_no_self)
    g.remove_node(node)
    return g


def get_provenance_dict(parameters=None):
    """
    Returns a dictionary encoding an execution of tskit conforming to the
    provenance schema.
    """
    document = {
        "schema_version": "1.0.0",
        "software": {"name": "ldgm", "version": provenance.__version__},
        "parameters": parameters,
    }
    return document


def check_bricked(ts):
    """
    Checks that ldgm.brick_ts() has been run on the input tree sequence.
    Returns True if so, False if not.
    """
    bricked = False
    for provenance in ts.provenances():
        if "brick_ts" in provenance.record:
            bricked = True

    return bricked


def prune_sites(ts, threshold):
    """
    Prune SNPs beneath a given threshold
    """
    a = ts.num_samples
    sites_to_delete = []
    for tree in ts.trees():
        for site in tree.sites():
            assert len(site.mutations) == 1
            freq = tree.num_samples(site.mutations[0].node) / a
            if freq < threshold or freq > 1 - threshold:
                sites_to_delete.append(site.id)
    return ts.delete_sites(sites_to_delete)


def identify_sites(bricked_ts):
    # Identifying sites only supports infinite sites
    assert bricked_ts.num_sites == bricked_ts.num_mutations
    bricks_to_muts = get_mut_edges(bricked_ts)
    identified_sites = collections.defaultdict(list)

    for brick, muts in bricks_to_muts.items():
        print(brick, muts)
        identified_sites[muts[0]] = muts
    return identified_sites


def return_site_list(bricked_ts, site_metadata_id=None):
    """
    Return list of SNPs in the tree sequence.
    :param TreeSequence bricked_ts: The input :class`tskit.TreeSequence`,
        from which the site list will be outputted. This tree sequence
        must have been "bricked", i.e. ldgm.brick_ts() was run on
        the tree sequence.
    :param string site_metadata_id: The site ID in the site metadata field
        in the tree sequence.
    :return: A tuple containing the index of sites in the tree sequence,
        the site ids (if specified, else None), the ancestral alleles for
        each site, the derived alleles for each site, the identified SNP
        for each site.
    :rtype: tuple of lists
    """
    assert check_bricked(bricked_ts)

    if site_metadata_id is not None:
        try:
            rsids = [
                json.loads(site.metadata)[site_metadata_id]
                for site in bricked_ts.sites()
            ]
        except json.decoder.JSONDecodeError:
            raise KeyError("Metadata is not JSON encoded")
        except KeyError:
            raise KeyError('"ID" is not a key in the metadata of this site')
    else:
        rsids = None

    identified_sites_dict = identify_sites(bricked_ts)
    identified_sites = np.full(bricked_ts.num_sites, -1)
    for site_id, site_targets in identified_sites_dict.items():
        for target in site_targets:
            identified_sites[target] = site_id
    anc_alleles = tskit.unpack_strings(
        bricked_ts.tables.sites.ancestral_state,
        bricked_ts.tables.sites.ancestral_state_offset,
    )
    derived_alleles = tskit.unpack_strings(
        bricked_ts.tables.mutations.derived_state,
        bricked_ts.tables.mutations.derived_state_offset,
    )
    assert len(anc_alleles) == len(derived_alleles) == bricked_ts.num_sites

    if rsids is not None:
        assert len(rsids) == len(anc_alleles)
    index = np.arange(0, bricked_ts.num_sites)
    return (index, rsids, anc_alleles, derived_alleles, identified_sites)
