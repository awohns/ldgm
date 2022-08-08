"""
Test cases for building a brick tree sequence
"""
import unittest

import ldgm
import pytest

from . import utility_functions


class TestBrickTreeSequence(unittest.TestCase):
    def test_modes(self):
        ts = utility_functions.single_tree_ts_n2()
        bricks = ldgm.bricks.Bricks(ts, recombination_freq_threshold=None)
        with pytest.raises(Exception):
            bricks.naive_split_edges(mode="wrong")
