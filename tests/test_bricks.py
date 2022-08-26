"""
Test cases for building a brick tree sequence
"""
import unittest

import ldgm
import pytest

from . import utility_functions


class TestBrickTreeSequence(unittest.TestCase):
    def test_no_bricks_added(self):
        for ts in [
            utility_functions.single_tree_ts_n2(),
            utility_functions.two_tree_ts(),
            utility_functions.two_tree_ts_with_unary_n3(),
        ]:
            bts = ldgm.brick_ts(ts, recombination_freq_threshold=None)
            assert ts.num_edges == bts.num_edges
            assert ts.num_nodes == bts.num_nodes

    def test_figure_one_example(self):
        ts = utility_functions.figure_one_example()
        bts = ldgm.brick_ts(ts, recombination_freq_threshold=None)
        assert bts.num_edges == 14
        assert ts.num_nodes == bts.num_nodes

    def test_supplementary_example(self):
        ts = utility_functions.supplementary_example()
        bts = ldgm.brick_ts(ts, recombination_freq_threshold=None)
        assert bts.num_edges == ts.num_edges + 4
        assert ts.num_nodes == bts.num_nodes
        # Frequency of recombination should be 1 / 5, so rec passes this threshold
        bts = ldgm.brick_ts(ts, recombination_freq_threshold=(1 / 6))
        assert bts.num_edges == ts.num_edges + 4
        assert ts.num_nodes == bts.num_nodes
        # Frequency of recombination should be 1 / 5, so rec does not pass this threshold
        bts = ldgm.brick_ts(ts, recombination_freq_threshold=(1 / 5))
        assert bts.num_edges == ts.num_edges
        assert ts.num_nodes == bts.num_nodes
