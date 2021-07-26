"""
Test cases for building the brick graph
"""
import unittest

import ld_graph
import msprime


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
        bricked = ld_graph.brick_ts(ts)
        number_of_labeled_bricks = len(ld_graph.utility.get_mut_edges(bricked).keys())
        assert ld_graph.reduce(bricked).number_of_nodes() == number_of_labeled_bricks
        bricked_graph = ld_graph.brick_graph(bricked)
        reduced_graph = ld_graph.reduce_graph(bricked_graph, bricked)
        assert reduced_graph.number_of_nodes() == number_of_labeled_bricks
