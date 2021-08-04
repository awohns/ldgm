"""
Create a brick graph
"""
import itertools

import networkx as nx
import numpy as np
from tqdm import tqdm

from . import utility


class BrickGraph:
    def __init__(self, bricked_ts):
        self.bricked_ts = bricked_ts
        self.brick_graph = nx.DiGraph()
        self.freqs = utility.get_brick_frequencies(self.bricked_ts)

    def find_odds(self, brick):
        return self.freqs[brick] / (1 - self.freqs[brick])

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
            child_brick, child_odds = self.up_vertex(node_edge_dict[child], "out")
            parent_brick, parent_odds = self.up_vertex(node_edge_dict[focal_node], "in")
            self.brick_graph.add_edge(
                child_brick, parent_brick, weight=child_odds / parent_odds
            )

            parent_brick, parent_odds = self.down_vertex(
                node_edge_dict[focal_node], "out"
            )
            child_brick, child_odds, self.down_vertex(node_edge_dict[child], "in")
            self.brick_graph.add_edge(
                parent_brick, child_brick, weight=child_odds / parent_odds
            )

        # Connect focal brick to its parent brick
        if edge.parent not in roots and focal_node not in roots:
            assert node_edge_dict[focal_node] != node_edge_dict[edge.parent]
            child_brick, child_odds = self.up_vertex(node_edge_dict[focal_node], "out")
            parent_brick, parent_odds = self.up_vertex(
                node_edge_dict[edge.parent], "in"
            )
            self.brick_graph.add_edge(
                child_brick, parent_brick, weight=child_odds / parent_odds
            )

            parent_brick, parent_odds = self.down_vertex(
                node_edge_dict[edge.parent], "out"
            )
            child_brick, child_odds, self.down_vertex(node_edge_dict[focal_node], "in")
            self.brick_graph.add_edge(
                parent_brick, child_brick, weight=child_odds / parent_odds
            )

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
                self.brick_graph.add_edge(
                    left_brick_up, right_brick_down, weight=left_odds * right_odds
                )

                right_brick_up, left_odds = self.up_vertex(
                    node_edge_dict[pair[1]], "out"
                )
                left_brick_down, right_odds = self.down_vertex(
                    node_edge_dict[pair[0]], "in"
                )
                self.brick_graph.add_edge(
                    right_brick_up, left_brick_down, weight=left_odds * right_odds
                )

    def make_connections(self, edge, tree2, node_edge_dict, index, prev_edge_dict):
        roots = tree2.roots
        children = tree2.children(edge.child)
        siblings = tree2.children(edge.parent)
        self.rule_one(edge, children, node_edge_dict, roots)
        # Two other ways (child1*child2)/parent freq > threshold
        # Only do it pairwise for children of higher frequency
        # (slightly different bc uses subset of pairs)
        # if tree2.num_samples(edge.parent)/tree2.num_samples() >= self.use_rule_two:
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
            self.brick_graph.add_edge(4 * brick + 1, 4 * brick, weight=1)

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

        # Total number of nodes should be less than (2 * number of labeled nodes) +
        # (4 * number of unlabeled nodes)
        # TODO: check why this isn't equal
        assert self.brick_graph.number_of_nodes() <= (
            2 * len(self.labeled_bricks)
        ) + 4 * len(self.unlabeled_bricks), (
            self.brick_graph.number_of_nodes(),
            (2 * len(self.labeled_bricks)),
            4 * len(self.unlabeled_bricks),
        )
        return self.brick_graph
