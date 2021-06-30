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
            1000,
            mutation_rate=1e-8,
            random_seed=13,
            recombination_rate=1e-8,
            record_full_arg=False,
            length=1e6,
            Ne=10000,
        )
        ld_graph.reduce(ts).number_of_nodes() == len(
            ld_graph.utility.get_mut_edges(ts).keys()
        )
