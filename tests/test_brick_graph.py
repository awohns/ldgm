"""
Test cases for building the brick graph
"""
import io
import unittest

import ld_graph
import numpy as np
import pytest
import tskit

from . import utility_functions


class TestExampleTrees(unittest.TestCase):
    def verify(self, ts):
        bts = ld_graph.brick_ts(ts, threshold=None, add_dummy_bricks=False)
        g = ld_graph.brick_graph(bts)
        self.check_rule_0(bts, g)
        self.check_rule_1(bts, g)
        self.check_out_nodes(g)

    def check_rule_0(self, brick_ts, brick_graph):
        """
        Check each brick is connected to its child haplotype
        With a weight equal to log(1).
        """
        graphical_model_nodes = brick_graph.nodes()
        unlabeled_bricks = []

        bricks_to_muts = ld_graph.utility.get_mut_edges(brick_ts)
        labeled_bricks = list(bricks_to_muts.keys())

        for edge in brick_ts.edges():
            if edge.id not in labeled_bricks:
                unlabeled_bricks.append(edge.id)
        # Check that out nodes exist for labeled bricks
        for labeled_brick in labeled_bricks:
            reindexed_brick = labeled_brick * 8
            assert reindexed_brick + 4 in graphical_model_nodes

        # Check that out nodes do not exist for unlabeled bricks
        for unlabeled_brick in unlabeled_bricks:
            reindexed_brick = unlabeled_brick * 8
            assert reindexed_brick + 4 not in graphical_model_nodes

    def return_node(self, brick_id, node_type):
        if node_type == "up_before":
            return brick_id * 8 + 0
        if node_type == "up_after":
            return brick_id * 8 + 1
        if node_type == "down_before":
            return brick_id * 8 + 2
        if node_type == "down_after":
            return brick_id * 8 + 3
        if node_type == "out":
            return brick_id * 8 + 4
        if node_type == "haplo before":
            return brick_id * 8 + 5
        if node_type == "haplo after":
            return brick_id * 8 + 6

    def check_rule_1(self, brick_ts, brick_graph):
        """
        Check that parent child bricks are connected
        Up after node of child should be connected to up after of parent
        Down after of parent should be connected to down after of child
        Up before of child to up before of parent (if child is unlabeled)
        Up before of child to up after of parent (if child is labeled)
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
                    brick = node_edge_dict[node]
                    brick_parent = node_edge_dict[tree.parent(node)]

                    # Up after of child to up after of parent
                    assert (
                        self.return_node(brick, "up_after"),
                        self.return_node(brick_parent, "up_after"),
                    ) in graphical_model_edges
                    # Down after of parent to down after of child
                    assert (
                        self.return_node(brick_parent, "down_after"),
                        self.return_node(brick, "down_after"),
                    ) in graphical_model_edges

                    # Up before of child to up before of parent (if child is unlabeled)
                    if brick in unlabeled_bricks:
                        assert (
                            self.return_node(brick, "up_before"),
                            self.return_node(brick_parent, "up_before"),
                        ) in graphical_model_edges
                    # Else if child is labeled, connect up before of child to up
                    # after of parent
                    elif brick in labeled_bricks:
                        assert (
                            self.return_node(brick, "up_before"),
                            self.return_node(brick_parent, "up_after"),
                        ) in graphical_model_edges

                    # If parent is unlabeled, down before of parent to down
                    # before of child
                    if brick_parent in unlabeled_bricks:
                        assert (
                            self.return_node(brick_parent, "down_before"),
                            self.return_node(brick, "down_before"),
                        ) in graphical_model_edges
                    # Else if parent is labeled, connect down before of parent
                    # to down after of chidl
                    elif brick_parent in labeled_bricks:
                        assert (
                            self.return_node(brick_parent, "down_before"),
                            self.return_node(brick, "down_after"),
                        ) in graphical_model_edges

                    # If parent brick is labeled, make out connections
                    if brick_parent in labeled_bricks:
                        # out of parent to down before of child
                        assert (
                            self.return_node(brick_parent, "out"),
                            self.return_node(brick, "down_before"),
                        ) in graphical_model_edges

                    # If child brick is labeled
                    if brick in labeled_bricks:
                        # Out of child to up before of parent
                        assert (
                            self.return_node(brick, "out"),
                            self.return_node(brick_parent, "up_before"),
                        ) in graphical_model_edges

    def check_out_nodes(self, brick_graph):
        # Check that no edge goes into an out node
        nodes = np.array(list(brick_graph.nodes()))
        l_out = nodes[nodes % 8 == 4]
        for edge in brick_graph.edges():
            assert edge[1] not in l_out

    def test_examples(self):
        for (
            name,
            val,
        ) in (
            utility_functions.__dict__.items()
        ):  # iterate through every module's attributes
            if name == "single_tree_ts_n2_dangling":
                with pytest.raises(ValueError):
                    self.verify(val())
            elif name == "two_tree_ts_with_unary_n3":
                with pytest.raises(ZeroDivisionError):
                    self.verify(val())
            elif name == "two_tree_ts_n2_part_dangling":
                with pytest.raises((ZeroDivisionError, ValueError)):
                    self.verify(val())
            elif callable(val):  # check if callable (normally functions)
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

    @pytest.mark.skip
    def test_dummy_example_one(self):
        ts = self.dummy_brick_test_one()
        bts_no_dummy = ld_graph.brick_ts(ts, threshold=None, add_dummy_bricks=False)
        brick_graph_no_dummy = ld_graph.brick_graph(bts_no_dummy)
        reduced_graph_no_dummy, _ = ld_graph.reduce_graph(
            brick_graph_no_dummy, bts_no_dummy, threshold=100
        )
        bts_dummy = ld_graph.brick_ts(ts, threshold=None, add_dummy_bricks=True)
        brick_graph_dummy = ld_graph.brick_graph(bts_dummy)
        reduced_graph_dummy, _ = ld_graph.reduce_graph(
            brick_graph_dummy, bts_dummy, threshold=100
        )
        assert reduced_graph_no_dummy.number_of_edges() == 0
        assert reduced_graph_dummy.number_of_edges() == 1
        assert (0, 1) in reduced_graph_dummy.edges()

    @pytest.mark.skip
    def test_dummy_example_two(self):
        ts = self.dummy_brick_test_two()
        bts_no_dummy = ld_graph.brick_ts(ts, threshold=None, add_dummy_bricks=False)
        brick_graph_no_dummy = ld_graph.brick_graph(bts_no_dummy)
        reduced_graph_no_dummy, _ = ld_graph.reduce_graph(
            brick_graph_no_dummy, bts_no_dummy, threshold=100
        )
        bts_dummy = ld_graph.brick_ts(ts, threshold=None, add_dummy_bricks=True)
        brick_graph_dummy = ld_graph.brick_graph(bts_dummy)
        reduced_graph_dummy, _ = ld_graph.reduce_graph(
            brick_graph_dummy, bts_dummy, threshold=100
        )
        assert reduced_graph_no_dummy.number_of_edges() == 1
        assert reduced_graph_dummy.number_of_edges() == 2
        assert (0, 1) in reduced_graph_no_dummy.edges()
        assert (0, 2) in reduced_graph_dummy.edges()
