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
        bricked = ldgm.brick_ts(ts, recombination_freq_threshold=None)
        number_of_labeled_bricks = len(ldgm.utility.get_mut_edges(bricked).keys())
        assert (
            ldgm.make_ldgm(ts, path_weight_threshold=100)[0].number_of_nodes()
            == number_of_labeled_bricks
        )
        bricked_graph = ldgm.brick_haplo_graph(bricked)
        reduced_graph = ldgm.reduce_graph(
            bricked_graph, bricked, path_weight_threshold=None
        )
        num_brick_nodes = np.sum(np.array(list(reduced_graph.nodes())) >= 0)
        assert num_brick_nodes == number_of_labeled_bricks
        max_haplotype = np.max(np.abs(np.array(list(reduced_graph.nodes()))))
        assert max_haplotype <= bricked.num_nodes


class TestMakeLdgm(unittest.TestCase):
    """
    Test ldgm.make_ldgm() works on a simple example
    """

    def two_labeled_nodes(self):
        """
        Simple example with two labeled nodes

        """
        ts = utility_functions.single_tree_ts_n2_2_mutations()
        bts = ldgm.brick_ts(
            ts,
            recombination_freq_threshold=None,
        )
        brick_haplo_graph = ldgm.brick_haplo_graph(bts)
        snp_grapher = ldgm.snp_graph.SNP_Graph(
            brick_haplo_graph, bts, path_weight_threshold=None
        )
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
            ldgm.make_ldgm(ts, path_weight_threshold=100)


class TestExamples(unittest.TestCase):
    """
    Test specific trees by hand
    """

    def test_fig1(self, num_processes=1):
        ts = utility_functions.figure_one_example()
        reduced, _ = ldgm.make_ldgm(
            ts, path_weight_threshold=100, num_processes=num_processes
        )
        assert nx.is_connected(reduced)
        assert reduced.number_of_edges() == 6
        return reduced

    def test_supplementary(self, num_processes=1):
        ts = utility_functions.supplementary_example()
        reduced, _ = ldgm.make_ldgm(
            ts, path_weight_threshold=100, num_processes=num_processes
        )
        edges = reduced.edges()
        assert (0, 1) in edges
        assert (1, 3) in edges
        assert (0, 2) in edges
        assert (0, 3) in edges
        assert (2, 3) in edges
        assert (1, 2) not in edges
        assert len(edges) == 5
        return reduced

    def test_triangle(self, num_processes=1):
        ts = utility_functions.triangle_example()
        reduced, _ = ldgm.make_ldgm(
            ts, path_weight_threshold=100, num_processes=num_processes
        )
        assert nx.is_connected(reduced)
        assert reduced.number_of_edges() == 3
        return reduced

    def test_multithreaded(self):
        for test in [self.test_fig1, self.test_supplementary, self.test_triangle]:
            singlethreaded = test(num_processes=1)
            multithreaded = test(num_processes=2)
            assert nx.is_isomorphic(singlethreaded, multithreaded)
        example_ts = msprime.simulate(
            25,
            mutation_rate=1e-8,
            recombination_rate=1e-8,
            Ne=10000,
            length=1e4,
            random_seed=1,
        )
        singlethreaded, _ = ldgm.make_ldgm(
            example_ts, path_weight_threshold=100, num_processes=1
        )
        multithreaded, _ = ldgm.make_ldgm(
            example_ts, path_weight_threshold=100, num_processes=5
        )
        assert nx.is_isomorphic(singlethreaded, multithreaded)


class TestPathWeightThreshold(unittest.TestCase):
    """
    Test no edges are greater than the given path_weight_threshold
    """

    def test_path_weight_threshold(self, num_processes=1):
        example_ts = msprime.simulate(
            50,
            mutation_rate=1e-8,
            recombination_rate=1e-8,
            Ne=10000,
            length=1e4,
            random_seed=1,
        )
        reduced_2, _ = ldgm.make_ldgm(example_ts, path_weight_threshold=2)
        edge_weights_2 = [
            reduced_2.get_edge_data(u, v)["weight"] for u, v in reduced_2.edges()
        ]
        assert np.max(edge_weights_2) < 2
        reduced_4, _ = ldgm.make_ldgm(example_ts, path_weight_threshold=4)
        edge_weights_4 = [
            reduced_4.get_edge_data(u, v)["weight"] for u, v in reduced_4.edges()
        ]
        assert np.max(edge_weights_4) < 4
        reduced_8, _ = ldgm.make_ldgm(example_ts, path_weight_threshold=8)
        edge_weights_8 = [
            reduced_8.get_edge_data(u, v)["weight"] for u, v in reduced_8.edges()
        ]
        assert np.max(edge_weights_8) < 8
