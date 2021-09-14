"""
Test cases for building a brick tree sequence
"""
import unittest

import ld_graph
import pytest

from . import utility_functions


class TestBrickTreeSequence(unittest.TestCase):
    def test_modes(self):
        ts = utility_functions.single_tree_ts_n2()
        bricks = ld_graph.bricks.Bricks(ts)
        with pytest.raises(Exception):
            bricks.naive_split_edges(mode="wrong")