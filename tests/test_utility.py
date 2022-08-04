"""
Test cases for utility code
"""
import unittest

import ldgm
import msprime
import pytest
import numpy as np
import json


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
        pruned_ts = ldgm.prune_sites(ts, threshold=0.02)
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


class TestCheckBricked(unittest.TestCase):
    """
    Test the check_bricked() function.
    """

    def test_check_bricked(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=1)
        bricked = ldgm.brick_ts(ts, threshold=None)
        assert ldgm.utility.check_bricked(bricked)


class TestReturnSNPList(unittest.TestCase):
    """
    Test the return_site_list() function.
    """

    def test_return_site_list(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=1)
        bricked = ldgm.brick_ts(ts, threshold=None)
        with pytest.raises(KeyError):
            ldgm.return_site_list(bricked, site_metadata_id="ID")
        (
            index,
            rsids,
            anc_alleles,
            derived_alleles,
        ) = ldgm.return_site_list(bricked)
        assert np.array_equal(index, np.array([0, 0, 2, 3, 0, 5, 6, 7, 8]))
        assert rsids == None
        assert np.array_equal(anc_alleles, np.full(bricked.num_sites, "0"))
        assert np.array_equal(derived_alleles, np.full(bricked.num_sites, "1"))

    def test_return_site_list_metadata(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=1)
        bricked = ldgm.brick_ts(ts, threshold=None)
        tables = bricked.dump_tables()
        sites_table = tables.sites.copy()
        tables.sites.clear()
        for site in sites_table:
            tables.sites.add_row(
                position=site.position,
                ancestral_state=site.ancestral_state,
                metadata=json.dumps({"ID": "an_id"}).encode(),
            )
        bricked_w_ids = tables.tree_sequence()
        with pytest.raises(KeyError):
            ldgm.return_site_list(bricked_w_ids, site_metadata_id="wrong_ID")
        (
            index,
            rsids,
            anc_alleles,
            derived_alleles,
        ) = ldgm.return_site_list(bricked_w_ids, site_metadata_id="ID")
        assert np.array_equal(index, np.array([0, 0, 2, 3, 0, 5, 6, 7, 8]))
        assert np.array_equal(rsids, np.full(bricked_w_ids.num_sites, "an_id"))
        assert np.array_equal(anc_alleles, np.full(bricked_w_ids.num_sites, "0"))
        assert np.array_equal(derived_alleles, np.full(bricked_w_ids.num_sites, "1"))
