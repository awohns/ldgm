"""
Functions to produce reduced graph and intermediate steps
"""
import networkx as nx

from . import brickgraph
from . import bricks
from . import reduction
from . import utility


def brick_ts(ts, threshold, add_dummy_bricks=False, progress=True):
    """
    Make a bricked tree sequence
    """
    brick = bricks.Bricks(ts, threshold, add_dummy_bricks, progress=progress)
    bricked = brick.naive_split_edges()
    return bricked


def brick_graph(brick_ts, threshold=None, use_rule_two=False, progress=True):
    """
    Make a brick graph
    """
    brick_grapher = brickgraph.BrickGraph(
        brick_ts, threshold, use_rule_two=use_rule_two, progress=progress
    )
    bricked_graph = brick_grapher.make_brick_graph()
    return bricked_graph


def reduce_graph(
    brick_graph, brick_ts, threshold, progress=True
):
    """
    Make a reduced graph from a brick graph and bricked tree sequence
    """
    snp_grapher = reduction.SNP_Graph(
        brick_graph, brick_ts, threshold, progress=progress
    )
    reduced_graph, mut_node, reach_star_sets = snp_grapher.create_reduced_graph()
    return reduced_graph, mut_node, reach_star_sets


def reduce(
    ts,
    path_threshold,
    recombination_threshold=None,
    progress=False,
):
    # Step 1: brick ts
    bts = brick_ts(ts, threshold=recombination_threshold, progress=progress)
    # Step 2 is brickgraph with no rule two
    bricked_graph = brick_graph(
        bts,
        path_threshold,
        progress=progress,
        use_rule_two=False,
    )
    # Step 3 compute reach* and create SNP-haplo graph
    reduced = reduce_graph(
        bricked_graph,
        bts,
        threshold=path_threshold,
        progress=progress,
    )
    H1 = reduced[0]
    # Step 4: brickgraph with rule two
    brickgraph_rule_two = brick_graph(
        bts, path_threshold, use_rule_two=True, progress=progress
    )
    # Step 5: reduce brickgraph created with rule two
    reduced_rule_two = reduce_graph(
        brickgraph_rule_two,
        bts,
        threshold=path_threshold,
        progress=progress,
    )
    H2 = reduced_rule_two[0]
    # Step 6: combine H1 and H2
    H_12 = nx.compose_all([H1, nx.reverse(H1), H2])
    # Remove nodes to reduce H1 and H2 to H_12
    H_12_reduced = H_12.copy()
    nodes = list(H_12_reduced.nodes())
    for node in nodes:
        if node < 0:
            H_12_reduced = utility.remove_node(H_12_reduced, node, path_threshold)
    H_12_reduced = H_12_reduced.to_undirected()

    return H_12_reduced, bts


def prune_snps(ts, threshold):
    return utility.prune_snps(ts, threshold)


def return_snp_list(ts):
    return utility.return_snp_list(ts)
