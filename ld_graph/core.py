"""
Functions to produce reduced graph and intermediate steps
"""
from . import brickgraph
from . import bricks
from . import reduction


def brick_ts(ts, threshold, add_dummy_bricks=False, progress=True):
    """
    Make a bricked tree sequence
    """
    brick = bricks.Bricks(ts, threshold, add_dummy_bricks, progress=progress)
    bricked = brick.naive_split_edges()
    return bricked


def brick_graph(brick_ts, threshold=None, progress=True):
    """
    Make a brick graph
    """
    brick_grapher = brickgraph.BrickGraph(brick_ts, threshold, progress=progress)
    bricked_graph = brick_grapher.make_brick_graph()
    return bricked_graph


def reduce_graph(brick_graph, brick_ts, threshold, progress=True):
    """
    Make a reduced graph from a brick graph and bricked tree sequence
    """
    snp_grapher = reduction.SNP_Graph(
        brick_graph, brick_ts, threshold, progress=progress
    )
    reduced_graph = snp_grapher.create_reduced_graph()
    return reduced_graph


def reduce(ts, threshold, brick_threshold=None, progress=True):
    bricked = brick_ts(ts, threshold=brick_threshold, progress=progress)
    bricked_graph = brick_graph(bricked, threshold, progress=progress)
    reduced_graph = reduce_graph(bricked_graph, bricked, threshold, progress=progress)
    return reduced_graph
