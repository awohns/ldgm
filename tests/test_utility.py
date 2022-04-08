"""
Test cases for utility code
"""
import unittest

import ld_graph
import msprime
import pytest


class TestRemoveNodes(unittest.TestCase):
    """
    Test the remove_nodes() function.
    """

    def test_is_directed(self):
        """
        remove_nodes() should fail if an undirected graph is passed
        """
        ts = msprime.simulate(10, mutation_rate=1, random_seed=1)
        reduced = ld_graph.reduce(ts, path_threshold=100)
        with pytest.raises(ValueError):
            ld_graph.utility.remove_node(reduced[0], 0, path_threshold=100)
