"""
Utility functions
"""
import collections
import itertools
import json
import pandas as pd

import numpy as np
import tskit
import networkx as nx

from . import provenance


def get_mut_edges(ts):
    """
    Returns a dictionary where keys are brick IDs and values are a list of
    mutations on that brick.

    NOTE: The first mutation on a brick is used as the node ID in the LDGM
    """
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
    # Ensure that these mutations are sorted by increasing ID
    for brick, muts in bricks_to_muts.items():
        assert np.array_equal(np.sort(muts), muts)
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


def remove_node(g, node, path_threshold):
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
    for recorded_provenance in ts.provenances():
        if "brick_ts" in recorded_provenance.record:
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


def identify_bricks(bricked_ts):
    # This function only supports infinite sites
    assert bricked_ts.num_sites == bricked_ts.num_mutations
    bricks_to_muts = get_mut_edges(bricked_ts)
    identified_bricks = collections.defaultdict(list)

    for brick_id, (brick, muts) in enumerate(bricks_to_muts.items()):
        identified_bricks[brick_id] = muts
    return identified_bricks


def make_snplist(bricked_ts, site_metadata_id=None, population_dict=None):
    """
    Return list of SNPs in the tree sequence.
    """
    assert check_bricked(bricked_ts)

    return_lists = {}
    if site_metadata_id is not None:
        try:
            site_ids = [
                json.loads(site.metadata)[site_metadata_id]
                for site in bricked_ts.sites()
            ]
        except json.decoder.JSONDecodeError:
            raise KeyError("Metadata is not JSON encoded")
        except KeyError:
            raise KeyError('"ID" is not a key in the metadata of this site')
        assert len(site_ids) == bricked_ts.num_sites
        return_lists["site_ids"] = np.array(site_ids)

    identified_bricks_dict = identify_bricks(bricked_ts)
    index = np.full(bricked_ts.num_sites, -1)
    for site_id, site_targets in identified_bricks_dict.items():
        for target in site_targets:
            index[target] = site_id
    return_lists["index"] = index
    anc_alleles = tskit.unpack_strings(
        bricked_ts.tables.sites.ancestral_state,
        bricked_ts.tables.sites.ancestral_state_offset,
    )
    deriv_alleles = tskit.unpack_strings(
        bricked_ts.tables.mutations.derived_state,
        bricked_ts.tables.mutations.derived_state_offset,
    )
    return_lists["anc_alleles"] = anc_alleles
    return_lists["deriv_alleles"] = deriv_alleles
    assert len(anc_alleles) == len(deriv_alleles) == bricked_ts.num_sites

    if population_dict is not None:
        for sample_set_id, (pop_id, sample_set) in enumerate(population_dict.items()):
            return_lists[pop_id] = np.full(bricked_ts.num_sites, -1, dtype=float)
            num_samples = len(sample_set)
            for tree in bricked_ts.trees(tracked_samples=sample_set):
                for site in tree.sites():
                    assert len(site.mutations) == 1
                    for mutation in site.mutations:
                        return_lists[pop_id][site.id] = (
                            tree.num_tracked_samples(mutation.node) / num_samples
                        )
            assert np.sum(return_lists[pop_id] == -1) == 0
    return pd.DataFrame(return_lists, columns=list(return_lists.keys()))


def convert_node_ids(reduced_graph, bricked_ts):
    bricks_to_muts = get_mut_edges(bricked_ts)
    old_to_new_ids = {}
    for new_brick, (old_brick, muts) in enumerate(bricks_to_muts.items()):
        for mut in muts:
            old_to_new_ids[mut] = new_brick
    return nx.relabel_nodes(reduced_graph, mapping=old_to_new_ids)


def return_edgelist(reduced_graph):
    """
    Function which takes an LDGM (the output of `ldgm.reduce()` or `ldgm.reduce_graph`
    and returns a list of edges.
    """
    edge_weights = []
    in_nodes = []
    out_nodes = []
    for u, v in reduced_graph.edges():
        edge_weights.append(reduced_graph.get_edge_data(u, v)["weight"])
        in_nodes.append(u)
        out_nodes.append(v)
    return pd.DataFrame(
        {"from": in_nodes, "to": out_nodes, "weight": np.round(edge_weights, 4)}
    )
