"""
Test cases for the command line interface for tsdate.
"""
import pathlib
import tempfile
from unittest import mock

import ldgm
import ldgm.cli as cli
import msprime
import networkx as nx


class TestLdgmArgParser:
    """
    Tests for the tsdate argument parser.
    """

    infile = "tmp.trees"
    output = "output"

    def test_default_values(self):
        with mock.patch("ldgm.cli.setup_logging"):
            parser = cli.ldgm_cli_parser()
            args = parser.parse_args(["reduce", self.infile, self.output, "4"])
        assert args.tree_sequence == self.infile
        assert args.output == self.output


class TestEndToEnd:
    """
    Class to test input to CLI outputs dated tree sequences.
    """

    def compare_python_api(self, input_ts):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_filename = pathlib.Path(tmpdir) / "input.trees"
            input_ts.dump(input_filename)
            output_filename = pathlib.Path(tmpdir) / "output.edgelist"
            full_cmd = "reduce " + str(input_filename) + f" {output_filename} " + " 4"
            cli.ldgm_main(full_cmd.split())
            output_network = nx.read_weighted_edgelist(output_filename, nodetype=int)
        snp_graph, bts = ldgm.make_ldgm(input_ts, 4)
        em = nx.algorithms.isomorphism.numerical_edge_match("weight", 1)
        print(output_network.edges(), snp_graph.edges())
        assert nx.is_isomorphic(output_network, snp_graph, edge_match=em)

    def test_compare_python_api(self):
        input_ts = msprime.simulate(
            50,
            Ne=10000,
            mutation_rate=1e-8,
            recombination_rate=1e-8,
            length=2e4,
            random_seed=10,
        )
        self.compare_python_api(input_ts)
