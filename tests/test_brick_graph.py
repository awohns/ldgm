"""
Test cases for building the brick graph
"""
import io
import unittest

import ld_graph
import numpy as np
import tskit

from . import utility_functions


class TestExampleTrees(unittest.TestCase):
    def verify(self, ts):
        bts = ld_graph.brick_ts(ts, add_dummy_bricks=True)
        g = ld_graph.brick_graph(bts)
        self.check_rule_0(bts, g)
        self.check_rule_1(bts, g)
        self.check_in_out_nodes(g)

    def check_rule_0(self, brick_ts, brick_graph):
        """
        Check down nodes are connected to up nodes at unlabeled bricks
        With a weight equal to log(1).
        Also checks that labeled bricks do not have an up or down node
        in the graphical model
        """
        graphical_model_edges = brick_graph.edges()
        graphical_model_nodes = brick_graph.nodes()
        unlabeled_bricks = []

        bricks_to_muts = ld_graph.utility.get_mut_edges(brick_ts)
        labeled_bricks = list(bricks_to_muts.keys())

        for edge in brick_ts.edges():
            if edge.id not in labeled_bricks:
                unlabeled_bricks.append(edge.id)

        for unlabeled_brick in unlabeled_bricks:
            reindexed_brick = unlabeled_brick * 4
            assert (reindexed_brick + 1, reindexed_brick + 0) in graphical_model_edges

        # Check that up down nodes do not exist for labeled bricks
        for labeled_brick in labeled_bricks:
            reindexed_brick = labeled_brick * 4
            assert reindexed_brick not in graphical_model_nodes
            assert reindexed_brick + 1 not in graphical_model_nodes

        # Check that in out nodes do not exist for unlabeled bricks
        for unlabeled_brick in unlabeled_bricks:
            reindexed_brick = unlabeled_brick * 4
            assert reindexed_brick + 2 not in graphical_model_nodes
            assert reindexed_brick + 3 not in graphical_model_nodes

    def check_rule_1(self, brick_ts, brick_graph):
        """
        Check that parent child bricks are connected
        Up node of child should be connected to up of parent
        Down of parent should be connected to down of child
        """
        graphical_model_edges = brick_graph.edges()
        unlabeled_bricks = []

        bricks_to_muts = ld_graph.utility.get_mut_edges(brick_ts)
        labeled_bricks = list(bricks_to_muts.keys())

        for edge in brick_ts.edges():
            if edge.id not in labeled_bricks:
                unlabeled_bricks.append(edge.id)

        node_edge_dict = {}
        for tree, (_, edges_out, edges_in) in zip(
            brick_ts.trees(), brick_ts.edge_diffs()
        ):
            for edge in edges_out:
                node_edge_dict.pop(edge.child)
            for edge in edges_in:
                node_edge_dict[edge.child] = edge.id
            for node in tree.nodes():
                if tree.parent(node) != -1 and tree.parent(node) != tree.root:
                    reindex_brick = 4 * node_edge_dict[node]
                    reindex_brick_parent = 4 * node_edge_dict[tree.parent(node)]
                    print(reindex_brick, reindex_brick_parent, graphical_model_edges)
                    if (
                        node_edge_dict[node] in unlabeled_bricks
                        and node_edge_dict[tree.parent(node)] in unlabeled_bricks
                    ):
                        # Up of child to up of parent
                        assert (
                            reindex_brick + 0,
                            reindex_brick_parent + 0,
                        ) in graphical_model_edges
                        # Down of parent to down of child
                        assert (
                            reindex_brick_parent + 1,
                            reindex_brick + 1,
                        ) in graphical_model_edges
                    elif (
                        node_edge_dict[node] in unlabeled_bricks
                        and node_edge_dict[tree.parent(node)] in labeled_bricks
                    ):
                        # Up of child to out of parent
                        assert (
                            reindex_brick + 0,
                            reindex_brick_parent + 2,
                        ) in graphical_model_edges
                        # In of parent to down of child
                        assert (
                            reindex_brick_parent + 3,
                            reindex_brick + 1,
                        ) in graphical_model_edges
                    elif (
                        node_edge_dict[node] in labeled_bricks
                        and node_edge_dict[tree.parent(node)] in unlabeled_bricks
                    ):
                        # In of child to up of parent
                        assert (
                            reindex_brick + 3,
                            reindex_brick_parent + 0,
                        ) in graphical_model_edges
                        # Down of parent to out of child
                        assert (
                            reindex_brick_parent + 1,
                            reindex_brick + 2,
                        ) in graphical_model_edges
                    elif (
                        node_edge_dict[node] in labeled_bricks
                        and node_edge_dict[tree.parent(node)] in labeled_bricks
                    ):
                        # In of child to out of parent
                        assert (
                            reindex_brick + 3,
                            reindex_brick_parent + 2,
                        ) in graphical_model_edges
                        # In of parent to out of child
                        assert (
                            reindex_brick_parent + 3,
                            reindex_brick + 2,
                        ) in graphical_model_edges
                    else:
                        raise ValueError

    def check_in_out_nodes(self, brick_graph):
        nodes = np.array(list(brick_graph.nodes()))
        l_in = nodes[nodes % 4 == 2]
        l_out = nodes[nodes % 4 == 3]
        for edge in brick_graph.edges():
            assert edge[0] not in l_in
            assert edge[1] not in l_out

    def test_examples(self):
        for (
            _,
            val,
        ) in (
            utility_functions.__dict__.items()
        ):  # iterate through every module's attributes
            if callable(val):  # check if callable (normally functions)
                self.verify(val())


class TestDummyBricks(unittest.TestCase):
    """
    Test that dummy bricks make expected connections
    """

    """
    Two example tree sequences to test dummy bricks
    """

    def dummy_brick_test_one(self):
        r"""
        Minimal example where a different reduced graph occurs with and without
        dummy bricks.
           1   |  2
           x 0 |  x 1
           0   |  0
        """
        nodes = io.StringIO(
            """\
        id      is_sample   time
        0       1           0
        1       0           1
        2       0           1
        """
        )
        edges = io.StringIO(
            """\
        left    right   parent  child
        0       0.5     1       0
        0.5     1       2       0
        """
        )
        sites = io.StringIO(
            """\
        position    ancestral_state
        0.2         0
        0.8         0
        """
        )
        mutations = io.StringIO(
            """\
        site    node    derived_state
        0       0       1
        1       0       1
        """
        )
        return tskit.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
        )

    def dummy_brick_test_two(self):
        r"""
        Small example where a different reduced graph occurs with and without dummy
        bricks.
           3   |  3
           x 1 |  |
           1   |  2
           x 0 |  x 2
           0   |  0
        """
        nodes = io.StringIO(
            """\
        id      is_sample   time
        0       1           0
        1       0           1
        2       0           1
        3       0           2
        """
        )
        edges = io.StringIO(
            """\
        left    right   parent  child
        0       0.5     1       0
        0.5     1       2       0
        0       0.5     3       1
        0.5     1       3       2
        """
        )
        sites = io.StringIO(
            """\
        position    ancestral_state
        0.2         0
        0.3         0
        0.8         0
        """
        )
        mutations = io.StringIO(
            """\
        site    node    derived_state
        0       0       1
        1       1       1
        2       0       1
        """
        )
        return tskit.load_text(
            nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
        )

    def test_dummy_example_one(self):
        ts = self.dummy_brick_test_one()
        bts_no_dummy = ld_graph.brick_ts(ts, add_dummy_bricks=False)
        brick_graph_no_dummy = ld_graph.brick_graph(bts_no_dummy)
        reduced_graph_no_dummy = ld_graph.reduce_graph(
            brick_graph_no_dummy, bts_no_dummy, threshold=100
        )
        bts_dummy = ld_graph.brick_ts(ts, add_dummy_bricks=True)
        brick_graph_dummy = ld_graph.brick_graph(bts_dummy)
        reduced_graph_dummy = ld_graph.reduce_graph(
            brick_graph_dummy, bts_dummy, threshold=100
        )
        assert reduced_graph_no_dummy.number_of_edges() == 0
        assert reduced_graph_dummy.number_of_edges() == 1
        assert (0, 1) in reduced_graph_dummy.edges()

    def test_dummy_example_two(self):
        ts = self.dummy_brick_test_two()
        bts_no_dummy = ld_graph.brick_ts(ts, add_dummy_bricks=False)
        brick_graph_no_dummy = ld_graph.brick_graph(bts_no_dummy)
        reduced_graph_no_dummy = ld_graph.reduce_graph(
            brick_graph_no_dummy, bts_no_dummy, threshold=100
        )
        bts_dummy = ld_graph.brick_ts(ts, add_dummy_bricks=True)
        brick_graph_dummy = ld_graph.brick_graph(bts_dummy)
        reduced_graph_dummy = ld_graph.reduce_graph(
            brick_graph_dummy, bts_dummy, threshold=100
        )
        assert reduced_graph_no_dummy.number_of_edges() == 1
        assert reduced_graph_dummy.number_of_edges() == 2
        assert (0, 1) in reduced_graph_no_dummy.edges()
        assert (0, 2) in reduced_graph_dummy.edges()
