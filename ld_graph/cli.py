"""
Command line interface for ld_graph.
"""
import argparse
import sys
import networkx as nx
import logging

import tskit

import ld_graph

logger = logging.getLogger(__name__)
log_format = "%(asctime)s %(levelname)s %(message)s"


def error_exit(message):
    """
    Exit with the specified error message, setting error status.
    """
    sys.exit("{}: {}".format(sys.argv[0], message))


def setup_logging(args):
    log_level = "WARN"
    if args.verbosity > 0:
        log_level = "INFO"
    if args.verbosity > 1:
        log_level = "DEBUG"
    logging.basicConfig(level=log_level, format=log_format)


def ld_graph_cli_parser():
    top_parser = argparse.ArgumentParser(
        description="This is the command line interface for ld_graph, a tool to create \
                graphical models of SNP dependencies from tree sequences."
    )
    top_parser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {ld_graph.__version__}"
    )

    subparsers = top_parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    parser = subparsers.add_parser(
        "reduce",
        help=("Takes a tree sequence and returns a SNP graph."),
    )

    parser.add_argument(
        "tree_sequence",
        help="The path and name of the input tree sequence from which \
                        we create the SNP graph.",
    )
    parser.add_argument(
        "output",
        help="The path and name of output file where the \
                        SNP graphical model will be stored.",
    )
    # TODO
    # parser.add_argument(
    #    "-p", "--progress", action="store_true", help="Show progress bar."
    # )
    # parser.add_argument(
    #    "-v", "--verbosity", type=int, default=0, help="How much verbosity to output."
    # )
    parser.set_defaults(runner=run_reduce)

    return top_parser


def run_reduce(args):
    try:
        ts = tskit.load(args.tree_sequence)
    except tskit.FileFormatError as ffe:
        error_exit(f"Error loading '{args.tree_sequence}: {ffe}")
    snp_graph = ld_graph.reduce(
        ts,
    )
    nx.readwrite.edgelist.write_edgelist(snp_graph, args.output)


def ld_graph_main(arg_list=None):
    parser = ld_graph_cli_parser()
    args = parser.parse_args(arg_list)
    setup_logging(args)
    args.runner(args)
