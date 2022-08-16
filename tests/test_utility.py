"""
Test cases for utility code
"""
import unittest

import ldgm
import msprime
import pytest
import numpy as np
import json

from . import utility_functions


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
        reduced = ldgm.make_ldgm(ts, path_weight_threshold=100)
        with pytest.raises(ValueError):
            ldgm.utility.remove_node(reduced[0], 0, path_threshold=100)


class TestCheckBricked(unittest.TestCase):
    """
    Test the check_bricked() function.
    """

    def test_check_bricked(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=1)
        bricked = ldgm.brick_ts(ts, recombination_freq_threshold=None)
        assert ldgm.utility.check_bricked(bricked)


class TestMakeSnpList(unittest.TestCase):
    """
    Test the make_snplist() function.
    """

    def test_make_snplist(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=1)
        bricked = ldgm.brick_ts(ts, recombination_freq_threshold=None)
        with pytest.raises(KeyError):
            ldgm.make_snplist(bricked, site_metadata_id="ID")
        results = ldgm.make_snplist(bricked)
        assert np.array_equal(results["index"], np.array([0, 0, 1, 2, 0, 3, 4, 5, 6]))
        assert np.array_equal(results["anc_alleles"], np.full(bricked.num_sites, "0"))
        assert np.array_equal(results["deriv_alleles"], np.full(bricked.num_sites, "1"))

    def test_make_snplist_metadata(self):
        ts = msprime.simulate(10, mutation_rate=1, random_seed=1)
        bricked = ldgm.brick_ts(ts, recombination_freq_threshold=None)
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
            ldgm.make_snplist(bricked_w_ids, site_metadata_id="wrong_ID")
        results = ldgm.make_snplist(bricked_w_ids, site_metadata_id="ID")
        assert np.array_equal(results["index"], np.array([0, 0, 1, 2, 0, 3, 4, 5, 6]))
        assert np.array_equal(
            results["site_ids"], np.full(bricked_w_ids.num_sites, "an_id")
        )
        assert np.array_equal(
            results["anc_alleles"], np.full(bricked_w_ids.num_sites, "0")
        )
        assert np.array_equal(
            results["deriv_alleles"], np.full(bricked_w_ids.num_sites, "1")
        )

    def test_population_frequencies(self):
        ts = msprime.simulate(
            population_configurations=[
                msprime.PopulationConfiguration(2),
                msprime.PopulationConfiguration(2),
            ],
            migration_matrix=[[0, 1], [1, 0]],
            recombination_rate=1,
            mutation_rate=1,
            random_seed=3,
        )
        bricked = ldgm.brick_ts(ts, recombination_freq_threshold=None)
        results = ldgm.make_snplist(
            bricked, population_dict={"pop0": [0, 1], "pop1": [2, 3]}
        )
        assert np.array_equal(results["index"], np.array([0, 1, 2, 3, 3, 3, 4]))
        assert np.array_equal(results["anc_alleles"], np.full(bricked.num_sites, "0"))
        assert np.array_equal(results["deriv_alleles"], np.full(bricked.num_sites, "1"))
        assert np.array_equal(
            results["pop0"],
            np.array([1, 0, 1, 0, 0, 0, 0]),
        )
        assert np.array_equal(
            results["pop1"],
            np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]),
        )

    def test_three_population_frequencies(self):
        ts = msprime.simulate(
            population_configurations=[
                msprime.PopulationConfiguration(2),
                msprime.PopulationConfiguration(2),
                msprime.PopulationConfiguration(3),
            ],
            migration_matrix=[[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
            recombination_rate=0.1,
            mutation_rate=0.1,
            random_seed=3,
        )
        bricked = ldgm.brick_ts(ts, recombination_freq_threshold=None)
        results = ldgm.make_snplist(
            bricked, population_dict={"pop0": [0, 1], "pop1": [2, 3], "pop2": [4, 5, 6]}
        )
        assert np.array_equal(results["index"], np.array([0, 0, 1, 0, 2, 3, 4]))
        assert np.array_equal(results["anc_alleles"], np.full(bricked.num_sites, "0"))
        assert np.array_equal(results["deriv_alleles"], np.full(bricked.num_sites, "1"))
        assert np.array_equal(
            results["pop0"], np.array([0.5, 0.5, 0, 0.5, 0, 0.5, 0.5])
        )
        assert np.array_equal(results["pop1"], np.array([0, 0, 0.5, 0, 0, 0, 0]))
        assert np.array_equal(results["pop2"], np.array([0, 0, 0, 0, 2 / 3, 1 / 3, 0]))


class TestReturnEdgelist(unittest.TestCase):
    """
    Test the return_edgelist() function.
    """

    ts = utility_functions.multiple_snps_branch()
    reduced, _ = ldgm.make_ldgm(ts, path_weight_threshold=100)
    edgelist = ldgm.return_edgelist(reduced)
    assert np.array_equal(edgelist["from"], np.array([2, 2, 1, 1, 0]))
    assert np.array_equal(edgelist["to"], np.array([3, 0, 0, 3, 3]))
    assert np.array_equal(
        edgelist["weight"], np.array([0.9808, 0.8109, 0.9808, 0.8109, 0.0000])
    )
