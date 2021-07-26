"""
Test cases for building the brick graph
"""
import unittest

import ld_graph

from . import utility_functions


class TestExampleTrees(unittest.TestCase):
    def verify(self, ts):
        bts = ld_graph.brick_ts(ts)
        g = ld_graph.brick_graph(bts)
        self.check_rule_0(bts, g)
        self.check_rule_1(bts, g)

    def check_rule_0(self, brick_ts, brick_graph):
        """
        Check left and right are connected to up node at unlabeled bricks
        Also checks that labeled bricks do not have an up, down, left or right node
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
            reindexed_brick = unlabeled_brick * 6
            assert (reindexed_brick + 2, reindexed_brick + 0) in graphical_model_edges
            assert (reindexed_brick + 3, reindexed_brick + 0) in graphical_model_edges

        # Check that up down left right do not exist for labeled bricks
        for labeled_brick in labeled_bricks:
            reindexed_brick = labeled_brick * 6
            assert reindexed_brick not in graphical_model_nodes
            assert reindexed_brick + 1 not in graphical_model_nodes
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
                    reindex_brick = 6 * node_edge_dict[node]
                    reindex_brick_parent = 6 * node_edge_dict[tree.parent(node)]
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
                            reindex_brick_parent + 5,
                        ) in graphical_model_edges
                        # In of parent to down of child
                        assert (
                            reindex_brick_parent + 4,
                            reindex_brick + 1,
                        ) in graphical_model_edges
                    elif (
                        node_edge_dict[node] in labeled_bricks
                        and node_edge_dict[tree.parent(node)] in unlabeled_bricks
                    ):
                        # In of child to up of parent
                        assert (
                            reindex_brick + 4,
                            reindex_brick_parent + 0,
                        ) in graphical_model_edges
                        # Down of parent to out of child
                        assert (
                            reindex_brick_parent + 1,
                            reindex_brick + 5,
                        ) in graphical_model_edges
                    elif (
                        node_edge_dict[node] in labeled_bricks
                        and node_edge_dict[tree.parent(node)] in labeled_bricks
                    ):
                        # In of child to out of parent
                        assert (
                            reindex_brick + 4,
                            reindex_brick_parent + 5,
                        ) in graphical_model_edges
                        # In of parent to out of child
                        assert (
                            reindex_brick_parent + 4,
                            reindex_brick + 5,
                        ) in graphical_model_edges
                    else:
                        raise ValueError

    def test_examples(self):
        for (
            name,
            val,
        ) in (
            utility_functions.__dict__.items()
        ):  # iterate through every module's attributes
            if callable(val):  # check if callable (normally functions)
                print(name)
                self.verify(val())
