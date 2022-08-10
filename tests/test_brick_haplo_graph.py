"""
Test cases for building the brick graph
"""
import io
import unittest

import ldgm
import numpy as np
import pytest
import msprime
import tskit

from . import utility_functions


class TestExampleTrees(unittest.TestCase):
    def verify(self, ts):
        bts = ldgm.brick_ts(ts, recombination_freq_threshold=None)
        g = ldgm.brick_haplo_graph(bts)
        self.check_rule_0(bts, g)
        self.check_rule_1(bts, g)
        self.check_out_nodes(g)

    def check_rule_0(self, brick_ts, brick_haplo_graph):
        """
        Check each brick is connected to its child haplotype
        With a weight equal to log(1).
        """
        graphical_model_nodes = brick_haplo_graph.nodes()
        unlabeled_bricks = []

        bricks_to_muts = ldgm.utility.get_mut_edges(brick_ts)
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

    def check_rule_1(self, brick_ts, brick_haplo_graph):
        """
        Check that parent child bricks are connected
        Up after node of child should be connected to up after of parent
        Down after of parent should be connected to down after of child
        Up before of child to up before of parent (if child is unlabeled)
        Up before of child to up after of parent (if child is labeled)
        """
        graphical_model_edges = brick_haplo_graph.edges()
        unlabeled_bricks = []

        bricks_to_muts = ldgm.utility.get_mut_edges(brick_ts)
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

    def check_out_nodes(self, brick_haplo_graph):
        # Check that no edge goes into an out node
        nodes = np.array(list(brick_haplo_graph.nodes()))
        l_out = nodes[nodes % 8 == 4]
        for edge in brick_haplo_graph.edges():
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

    def test_triangle_brickgraph(self):
        ts = utility_functions.triangle_example()
        bts = ldgm.brick_ts(ts)
        brick_haplo_graph_wo_sibs = ldgm.brick_haplo_graph(bts, make_sibs=False)
        brick_haplo_graph_w_sibs = ldgm.brick_haplo_graph(bts, make_sibs=True)
        edges_wo_sibs = list(brick_haplo_graph_wo_sibs.edges())
        edges_w_sibs = list(brick_haplo_graph_w_sibs.edges())
        for name, edges in {"nosibs": edges_wo_sibs, "sibs": edges_w_sibs}.items():
            # Rule 0 labeled bricks
            assert ((8 * 0) + 2, (8 * 0) + 7) in edges
            assert ((8 * 0) + 3, (8 * 0) + 7) in edges
            assert ((8 * 0) + 4, (8 * 0) + 6) in edges
            assert ((8 * 1) + 2, (8 * 1) + 7) in edges
            assert ((8 * 1) + 3, (8 * 1) + 7) in edges
            assert ((8 * 1) + 4, (8 * 1) + 6) in edges
            assert ((8 * 2) + 2, (8 * 2) + 7) in edges
            assert ((8 * 2) + 3, (8 * 2) + 7) in edges
            assert ((8 * 2) + 4, (8 * 2) + 6) in edges
            # Rule 0 unlabeled brick
            assert ((8 * 3) + 2, (8 * 3) + 6) in edges
            assert ((8 * 3) + 3, (8 * 3) + 7) in edges
            # Rule One Brick 0 to 2
            assert ((8 * 0) + 0, (8 * 2) + 1) in edges
            assert ((8 * 0) + 1, (8 * 2) + 1) in edges
            assert ((8 * 0) + 4, (8 * 2) + 0) in edges
            assert ((8 * 2) + 2, (8 * 0) + 3) in edges
            assert ((8 * 2) + 3, (8 * 0) + 3) in edges
            assert ((8 * 2) + 4, (8 * 0) + 2) in edges
            # Rule One Brick 1 to 2
            assert ((8 * 1) + 0, (8 * 2) + 1) in edges
            assert ((8 * 1) + 1, (8 * 2) + 1) in edges
            assert ((8 * 1) + 4, (8 * 2) + 0) in edges
            assert ((8 * 2) + 2, (8 * 1) + 3) in edges
            assert ((8 * 2) + 3, (8 * 1) + 3) in edges
            assert ((8 * 2) + 4, (8 * 1) + 2) in edges

            if name == "nosibs":
                assert len(edges) == 23
            else:
                # Rule One Uturn connections
                assert ((8 * 0) + 4, (8 * 2) + 5) in edges
                assert ((8 * 2) + 5, (8 * 0) + 2) in edges
                assert ((8 * 1) + 4, (8 * 2) + 5) in edges
                assert ((8 * 2) + 5, (8 * 1) + 2) in edges
                # Rule Two Brick 0 to 1
                assert ((8 * 0) + 0, (8 * 1) + 3) in edges
                assert ((8 * 0) + 1, (8 * 1) + 3) in edges
                assert ((8 * 0) + 4, (8 * 1) + 2) in edges
                assert ((8 * 1) + 0, (8 * 0) + 3) in edges
                assert ((8 * 1) + 1, (8 * 0) + 3) in edges
                assert ((8 * 1) + 4, (8 * 0) + 2) in edges
                # Rule Two Brick 2 to 3
                assert ((8 * 2) + 0, (8 * 3) + 3) in edges
                assert ((8 * 2) + 1, (8 * 3) + 3) in edges
                assert ((8 * 2) + 4, (8 * 3) + 2) in edges
                assert ((8 * 3) + 0, (8 * 2) + 2) in edges
                assert ((8 * 3) + 1, (8 * 2) + 3) in edges

                assert len(edges) == 38


class TestEdgeWeightThreshold(unittest.TestCase):
    """
    Test that no edges larger than the given edge weight threshold
    end up in the brick haplo graph.
    """

    def test_edge_weight_threshold(self):
        ts = msprime.simulate(
            100,
            mutation_rate=1e-8,
            recombination_rate=1e-8,
            Ne=10000,
            length=1e5,
            random_seed=1,
        )
        bts = ldgm.brick_ts(ts)
        brick_haplo_graph_2 = ldgm.brick_haplo_graph(bts, edge_weight_threshold=2)
        edge_weights_2 = [
            brick_haplo_graph_2.get_edge_data(u, v)["weight"]
            for u, v in brick_haplo_graph_2.edges()
        ]
        assert np.max(edge_weights_2) < 2
        brick_haplo_graph_4 = ldgm.brick_haplo_graph(bts, edge_weight_threshold=4)
        edge_weights_4 = [
            brick_haplo_graph_4.get_edge_data(u, v)["weight"]
            for u, v in brick_haplo_graph_4.edges()
        ]
        assert np.max(edge_weights_4) < 4
        brick_haplo_graph_8 = ldgm.brick_haplo_graph(bts, edge_weight_threshold=8)
        edge_weights_8 = [
            brick_haplo_graph_8.get_edge_data(u, v)["weight"]
            for u, v in brick_haplo_graph_8.edges()
        ]
        assert np.max(edge_weights_8) < 8
