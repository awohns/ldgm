"""
Create a brick graph
"""
import itertools

import networkx as nx
import numpy as np
from tqdm import tqdm

from . import utility


class BrickGraph:
    def __init__(self, bricked_ts, threshold):
        self.bricked_ts = bricked_ts
        self.threshold = threshold
        self.brick_graph = nx.DiGraph()
        self.freqs = utility.get_brick_frequencies(self.bricked_ts)

    def find_odds(self, brick):
        return self.freqs[brick] / (1 - self.freqs[brick])

    def log_odds(self, odds):
        if odds != 1:
            return np.log(odds) * -1
        else:
            return np.log(odds)

    def add_edge_threshold(self, from_node, to_node, weight):
        if weight <= 0:
            weight = 0
        if self.threshold is not None:
            if weight < self.threshold:
                self.brick_graph.add_edge(from_node, to_node, weight=weight)
        else:
            self.brick_graph.add_edge(from_node, to_node, weight=weight)

    # make an argument for in and out here, so we know how to split if it's labeled
    def up_vertex(self, edge_id, in_out):
        odds = self.find_odds(edge_id)
        if edge_id in self.labeled_bricks:
            if in_out == "in":
                return 4 * edge_id + 2, odds
            else:
                return 4 * edge_id + 3, odds
        else:
            return 4 * edge_id, odds

    def down_vertex(self, edge_id, in_out):
        odds = self.find_odds(edge_id)
        if edge_id in self.labeled_bricks:
            if in_out == "in":
                return 4 * edge_id + 2, odds
            else:
                return 4 * edge_id + 3, odds
        else:
            return 4 * edge_id + 1, odds

    def rule_one(self, edge, children, node_edge_dict, roots):
        """
        Rule 1:
        Connect focal brick (indexed by edge.child) to its child bricks
        """
        focal_node = edge.child
        for child in children:
            assert node_edge_dict[child] != node_edge_dict[focal_node]
            # Up of child to up of parent
            child_label, child_odds = self.up_vertex(node_edge_dict[child], "out")
            parent_label, parent_odds = self.up_vertex(node_edge_dict[focal_node], "in")
            weight = self.log_odds(child_odds / parent_odds)
            self.add_edge_threshold(child_label, parent_label, weight)

            # Down of parent to down of child
            parent_label, parent_odds = self.down_vertex(
                node_edge_dict[focal_node], "out"
            )
            child_label, child_odds = self.down_vertex(node_edge_dict[child], "in")
            weight = self.log_odds(child_odds / parent_odds)
            self.add_edge_threshold(parent_label, child_label, weight)

        # Connect focal brick to its parent brick
        if edge.parent not in roots and focal_node not in roots:
            assert node_edge_dict[focal_node] != node_edge_dict[edge.parent]
            child_label, child_odds = self.up_vertex(node_edge_dict[focal_node], "out")
            parent_label, parent_odds = self.up_vertex(
                node_edge_dict[edge.parent], "in"
            )
            weight = self.log_odds(child_odds / parent_odds)
            self.add_edge_threshold(child_label, parent_label, weight)

            parent_label, parent_odds = self.down_vertex(
                node_edge_dict[edge.parent], "out"
            )
            child_label, child_odds = self.down_vertex(node_edge_dict[focal_node], "in")
            weight = self.log_odds(child_odds / parent_odds)
            self.add_edge_threshold(parent_label, child_label, weight)

    def rule_two(self, edge, siblings, node_edge_dict):
        # Rule 2: Connect focal brick to its siblings
        if len(siblings) > 1:
            for pair in itertools.combinations(siblings, 2):
                left_brick_up, left_odds = self.up_vertex(
                    node_edge_dict[pair[0]], "out"
                )
                right_brick_down, right_odds = self.down_vertex(
                    node_edge_dict[pair[1]], "in"
                )
                weight = self.log_odds(left_odds * right_odds)
                self.add_edge_threshold(left_brick_up, right_brick_down, weight)

                right_brick_up, left_odds = self.up_vertex(
                    node_edge_dict[pair[1]], "out"
                )
                left_brick_down, right_odds = self.down_vertex(
                    node_edge_dict[pair[0]], "in"
                )
                weight = self.log_odds(left_odds * right_odds)
                self.add_edge_threshold(right_brick_up, left_brick_down, weight)

    def make_connections(self, edge, tree2, node_edge_dict, index, prev_edge_dict):
        roots = tree2.roots
        children = tree2.children(edge.child)
        siblings = tree2.children(edge.parent)
        self.rule_one(edge, children, node_edge_dict, roots)
        self.rule_two(edge, siblings, node_edge_dict)

    def make_brick_graph(self):
        """
        Given a "bricked" ts, connect bricks according to three rules:
        1. Connect parent and child bricks
        2. Connect sibling bricks
        3. Connect bricks which have two different adjacent parents
        """
        bricks_to_muts = utility.get_mut_edges(self.bricked_ts)
        self.labeled_bricks = np.array(list(bricks_to_muts.keys()))
        # The number of labeled bricks should be less than or equal to the number of SNPs
        # as well as the number of bricks
        assert len(self.labeled_bricks) <= self.bricked_ts.num_mutations
        assert len(self.labeled_bricks) <= self.bricked_ts.num_edges
        bricks = np.arange(0, self.bricked_ts.num_edges)
        assert len(bricks) == self.bricked_ts.num_edges
        self.unlabeled_bricks = bricks[~np.isin(bricks, self.labeled_bricks)]
        # There should be no overlap between labeled and unlabeled bricks
        assert len(np.intersect1d(self.labeled_bricks, self.unlabeled_bricks)) == 0
        # The number of unlabeled bricks should be greater than the number of edges minus
        # the number of SNPs
        assert len(self.unlabeled_bricks) >= (
            self.bricked_ts.num_edges - self.bricked_ts.num_mutations
        ), (len(bricks), self.bricked_ts.num_mutations, len(self.unlabeled_bricks))

        node_edge_dict = {}

        # Rule Zero
        # For unlabeled nodes: connect down to up (within a brick)
        for brick in self.unlabeled_bricks:
            self.brick_graph.add_edge(4 * brick + 1, 4 * brick, weight=self.log_odds(1))

        for index, (tree2, (_, edges_out, edges_in)) in tqdm(
            enumerate(zip(self.bricked_ts.trees(), self.bricked_ts.edge_diffs())),
            desc="Brick graph: iterate over edges",
            total=self.bricked_ts.num_trees,
        ):
            prev_edge_dict = node_edge_dict.copy()

            for edge in edges_out:
                node_edge_dict.pop(edge.child)
            for edge in edges_in:
                node_edge_dict[edge.child] = edge.id

            for edge in edges_in:
                self.make_connections(
                    edge,
                    tree2,
                    node_edge_dict,
                    index,
                    prev_edge_dict,
                )

        return self.brick_graph
