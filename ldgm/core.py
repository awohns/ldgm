"""
Functions to produce reduced graph and intermediate steps
"""
import networkx as nx
from tqdm import tqdm

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


def brick_graph(brick_ts, threshold=None, make_sibs=True, progress=True):
    """
    Make a brick graph
    """
    brick_grapher = brickgraph.BrickGraph(
        brick_ts, threshold, make_sibs=make_sibs, progress=progress
    )
    bricked_graph = brick_grapher.make_brick_graph()
    return bricked_graph


def reduce_graph(
    brick_graph, brick_ts, threshold, num_processes=1, progress=True, chunksize=100
):
    """
    Make a reduced graph from a brick graph and bricked tree sequence
    """
    snp_grapher = reduction.SNP_Graph(
        brick_graph,
        brick_ts,
        threshold,
        num_processes=num_processes,
        progress=progress,
        chunksize=chunksize,
    )
    reduced_graph, mut_node, reach_star_sets = snp_grapher.create_reduced_graph()
    return reduced_graph, mut_node, reach_star_sets


def reduce(
    ts,
    path_threshold,
    recombination_threshold=None,
    use_softmin=False,
    num_processes=1,
    chunksize=100,
    progress=False,
):
    """
    Run the entire pipeline from tree sequence to reduced graph.

    :param TreeSequence tree_sequence: The input :class:`tskit.TreeSequence`,
        which will be used to create the linkage disequilibrium graphical
        model.
    :param float path_threshold: The maximum path threshold to retain in
        the linkage disequilibrium graphical model.
    :param float recombination_threshold: Defines the minimum frequency of
        recombination events used to create bricks in the tree sequence.
        Default: None (any recombination creates bricks).
    :param bool use_softmin: If True, use the softmin function to combine
        path weights when removing nodes. If False, use the minimum of path
        weights. Default: False
    :param bool progress: Whether to display a progress bar. Default: False
    """
    # Step 1: brick ts
    bts = brick_ts(ts, threshold=recombination_threshold, progress=progress)
    # Step 2: is brickgraph with no rule two or uturns
    bricked_graph = brick_graph(
        bts,
        path_threshold,
        progress=progress,
        make_sibs=False,
    )
    # Step 3: compute reach* and create SNP-haplo graph
    H1, _, _ = reduce_graph(
        bricked_graph,
        bts,
        threshold=path_threshold,
        num_processes=num_processes,
        chunksize=chunksize,
        progress=progress,
    )
    # Step 4: brickgraph with rule two and uturns
    brickgraph_rule_two = brick_graph(
        bts, path_threshold, make_sibs=True, progress=progress
    )
    # Step 5: reduce brickgraph created with rule two
    H2, _, _ = reduce_graph(
        brickgraph_rule_two,
        bts,
        threshold=path_threshold,
        num_processes=num_processes,
        progress=progress,
        chunksize=chunksize,
    )
    # Step 6: combine H1 and H2
    H_12 = nx.compose_all([H1, nx.reverse(H1), H2])
    # Step 7: Remove haplotype vertices to reduce H_12
    H_12_reduced = H_12.copy()
    nodes = list(H_12_reduced.nodes())
    for node in tqdm(nodes, total=len(nodes), desc="Removing nodes"):
        if node < 0:
            H_12_reduced = utility.remove_node(
                H_12_reduced, node, path_threshold, use_softmin=use_softmin
            )
    H_12_reduced = H_12_reduced.to_undirected()

    return H_12_reduced, bts


def prune_snps(ts, threshold):
    return utility.prune_snps(ts, threshold)


def return_snp_list(ts):
    return utility.return_snp_list(ts)
