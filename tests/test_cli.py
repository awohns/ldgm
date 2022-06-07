"""
Test cases for the command line interface for tsdate.
"""
from unittest import mock

import ldgm.cli as cli


class TestLdgmArgParser:
    """
    Tests for the tsdate argument parser.
    """

    infile = "tmp.trees"
    output = "output"

    def test_default_values(self):
        with mock.patch("ldgm.cli.setup_logging"):
            parser = cli.ldgm_cli_parser()
            args = parser.parse_args(["reduce", self.infile, self.output])
        assert args.tree_sequence == self.infile
        assert args.output == self.output
