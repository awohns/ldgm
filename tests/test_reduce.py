"""
Test cases for building the brick graph
"""
import unittest

import ldgm
import numpy as np
import msprime
import networkx as nx
import pytest

from . import utility_functions


class TestNumNodes(unittest.TestCase):
    """
    Test the reduced graph node number is correct
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
        bricked = ldgm.brick_ts(ts, threshold=None)
        number_of_labeled_bricks = len(ldgm.utility.get_mut_edges(bricked).keys())
        assert (
            ldgm.reduce(ts, path_threshold=100)[0].number_of_nodes()
            == number_of_labeled_bricks
        )
        bricked_graph = ldgm.brick_graph(bricked)
        reduced_graph, _, _ = ldgm.reduce_graph(
            bricked_graph, bricked, threshold=None
        )
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
        bts = ldgm.brick_ts(ts, threshold=None, add_dummy_bricks=False)
        brick_graph = ldgm.brick_graph(bts)
        snp_grapher = ldgm.snp_graph.SNP_Graph(brick_graph, bts, threshold=None)
        reduced_graph = snp_grapher.create_reduced_graph()

        manual_graph = nx.Graph()
        manual_graph.add_edge(0, 1)
        assert nx.is_isomorphic(reduced_graph, manual_graph)

    def test_no_mutations(self):
        """
        Test fails with no mutations
        """

        ts = msprime.simulate(10)
        with pytest.raises(ValueError):
            ldgm.reduce(ts, path_threshold=100)


class TestExamples(unittest.TestCase):
    """
    Test specific trees by hand
    """

    def test_fig1(self, num_threads=1):
        ts = utility_functions.figure_one_example()
        reduced = ldgm.reduce(ts, path_threshold=100)
        assert nx.is_connected(reduced[0])
        assert reduced[0].number_of_edges() == 6
        return reduced[0]

    def test_supplementary(self, num_threads=1):
        ts = utility_functions.supplementary_example()
        reduced = ldgm.reduce(ts, path_threshold=100)
        edges = reduced[0].edges()
        assert (0, 1) in edges
        assert (1, 3) in edges
        assert (0, 2) in edges
        assert (0, 3) in edges
        assert (2, 3) in edges
        assert (1, 2) not in edges
        assert len(edges) == 5
        return reduced[0]

    def test_triangle(self, num_threads=1):
        ts = utility_functions.triangle_example()
        reduced = ldgm.reduce(ts, path_threshold=100)
        assert nx.is_connected(reduced[0])
        assert reduced[0].number_of_edges() == 3
        return reduced[0]

    def test_multithreaded(self):
        for test in [self.test_fig1, self.test_supplementary, self.test_triangle]:
            singlethreaded = test(num_threads=1)
            multithreaded = test(num_threads=2)
            assert nx.is_isomorphic(singlethreaded, multithreaded)
