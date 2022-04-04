"""
Test cases for building the brick graph
"""
import unittest

import ld_graph
import numpy as np
import msprime
import networkx as nx

from . import utility_functions


class TestNumNodes(unittest.TestCase):
    """
    Test for some of the basic functions used in tsdate
    """

    def test_num_nodes(self):
        """
        The reduced graph should have the same number of nodes as labeled bricks
        """
        ts = msprime.simulate(
            100,
            mutation_rate=1e-8,
            random_seed=13,
            recombination_rate=1e-8,
            record_full_arg=False,
            length=1e4,
            Ne=10000,
        )
        bricked = ld_graph.brick_ts(ts, threshold=None)
        number_of_labeled_bricks = len(ld_graph.utility.get_mut_edges(bricked).keys())
        assert (
            ld_graph.reduce(ts, path_threshold=100)[0].number_of_nodes()
            == number_of_labeled_bricks
        )
        bricked_graph = ld_graph.brick_graph(bricked)
        reduced_graph, _, _ = ld_graph.reduce_graph(bricked_graph, bricked, threshold=None)
        num_brick_nodes = np.sum(np.array(list(reduced_graph.nodes())) >= 0)
        assert num_brick_nodes == number_of_labeled_bricks
        max_haplotype = np.max(np.abs(np.array(list(reduced_graph.nodes()))))
        assert max_haplotype <= bricked.num_nodes


class TestReduce(unittest.TestCase):
    """
    Test reduce works on a simple example
    """

    def two_labeled_nodes(self):
        """
        Simple example with two labeled nodes

        """
        ts = utility_functions.single_tree_ts_n2_2_mutations()
        bts = ld_graph.brick_ts(ts, threshold=None, add_dummy_bricks=False)
        brick_graph = ld_graph.brick_graph(bts)
        snp_grapher = ld_graph.snp_graph.SNP_Graph(brick_graph, bts, threshold=None)
        reduced_graph = snp_grapher.create_reduced_graph()

        manual_graph = nx.Graph()
        manual_graph.add_edge(0, 1)
        assert nx.is_isomorphic(reduced_graph, manual_graph)
