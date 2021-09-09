"""
Functions to produce reduced graph and intermediate steps
"""
from . import bricks
from . import bricks_graph
from . import regularization
from . import snp_graph


def brick_ts(ts, add_dummy_bricks=True):
    """
    Make a bricked tree sequence
    """
    brick = bricks.Bricks(ts, add_dummy_bricks)
    bricked = brick.naive_split_edges()
    return bricked


def time_regularize(ts, time):
    """
    Regularize a tree sequence by time
    """
    return regularization.Regularize(ts).time_regularize(time)


def frequency_regularize(ts, freq, **kwargs):
    """
    Regularize a tree sequence by frequency
    """
    return regularization.Regularize(ts).frequency_regularize(freq, **kwargs)


def brick_graph(brick_ts, threshold=None):
    """
    Make a brick graph
    """
    brick_grapher = bricks_graph.BrickGraph(brick_ts, threshold)
    bricked_graph = brick_grapher.make_brick_graph()
    return bricked_graph


def reduce_graph(brick_graph, brick_ts, threshold):
    """
    Make a reduced graph from a brick graph and bricked tree sequence
    """
    snp_grapher = snp_graph.SNP_Graph(brick_graph, brick_ts, threshold)
    reduced_graph = snp_grapher.create_reduced_graph()
    return reduced_graph


def reduce(ts, threshold):
    bricked = brick_ts(ts)
    bricked_graph = brick_graph(bricked, threshold)
    reduced_graph = reduce_graph(bricked_graph, bricked, threshold)
    return reduced_graph
