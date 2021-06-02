"""
Functions to produce reduced graph and intermediate steps
"""
from . import bricks
from . import bricks_graph
from . import snp_graph


def brick_ts(ts):
    """
    Make a bricked tree sequence
    """
    brick = bricks.Bricks(ts)
    bricked = brick.naive_split_edges()
    return bricked


def brick_graph(brick_ts):
    """
    Make a brick graph
    """
    brick_grapher = bricks_graph.BrickGraph(brick_ts)
    bricked_graph = brick_grapher.make_brick_graph()
    return bricked_graph


def reduce_graph(brick_graph, brick_ts):
    """
    Make a reduced graph from a brick graph and bricked tree sequence
    """
    snp_grapher = snp_graph.SNP_Graph(brick_graph, brick_ts)
    reduced_graph, id_to_muts = snp_grapher.create_reduced_graph()
    return reduced_graph, id_to_muts


def reduce(ts):
    bricked = brick_ts(ts)
    bricked_graph = brick_graph(bricked)
    reduced_graph, id_to_muts = reduce_graph(bricked_graph, bricked)
    return reduced_graph, id_to_muts
