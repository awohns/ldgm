"""
Test cases for utility code
"""
import unittest

import ldgm
import msprime
import pytest
import numpy as np


class TestDummyBricks(unittest.TestCase):
    """
    Tests adding dummy bricks
    """

    def test_add_dummy_bricks(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=1)
        bts = ldgm.brick_ts(ts, threshold=None, add_dummy_bricks=True)
        assert bts.num_nodes == ts.num_nodes + ts.num_samples


class TestPruneSnps(unittest.TestCase):
    """
    Test that after pruning SNPs, no SNPs below the given frequency remain
    """

    def test_removing_edges(self):
        ts = msprime.simulate(
            100,
            mutation_rate=1e-8,
            recombination_rate=1e-8,
            Ne=10000,
            length=1e5,
            random_seed=1,
        )
        genos = ts.genotype_matrix()
        assert np.any(np.sum(genos, axis=1) == 1)
        pruned_ts = ldgm.prune_snps(ts, threshold=0.02)
        pruned_genos = pruned_ts.genotype_matrix()
        assert pruned_genos.shape[0] < genos.shape[0]
        assert ~np.any(np.sum(genos, axis=0) <= 1)


class TestRemoveNodes(unittest.TestCase):
    """
    Test the remove_nodes() function.
    """

    def test_is_directed(self):
        """
        remove_nodes() should fail if an undirected graph is passed
        """
        ts = msprime.simulate(10, mutation_rate=1, random_seed=1)
        reduced = ldgm.reduce(ts, path_threshold=100)
        with pytest.raises(ValueError):
            ldgm.utility.remove_node(reduced[0], 0, path_threshold=100)
