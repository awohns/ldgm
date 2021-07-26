"""
Create a brick graph
"""
import itertools

import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

from . import utility


class BrickGraph:
    def __init__(self, bricked_ts, use_rule_two=True):
        self.bricked_ts = bricked_ts
        self.from_to_set = set()
        self.use_rule_two = use_rule_two

    # make an argument for in and out here, so we know how to split if it's labeled
    def up_vertex(self, edge_id, in_out):
        if edge_id in self.labeled_bricks:
            if in_out == "in":
                return 6 * edge_id + 4
            else:
                return 6 * edge_id + 5
        else:
            return 6 * edge_id

    def down_vertex(self, edge_id, in_out):
        if edge_id in self.labeled_bricks:
            if in_out == "in":
                return 6 * edge_id + 4
            else:
                return 6 * edge_id + 5
        else:
            return 6 * edge_id + 1

    def left_vertex(self, edge_id, in_out):
        if edge_id in self.labeled_bricks:
            if in_out == "in":
                return 6 * edge_id + 4
            else:
                return 6 * edge_id + 5
        else:
            return 6 * edge_id + 2

    def right_vertex(self, edge_id, in_out):
        if edge_id in self.labeled_bricks:
            if in_out == "in":
                return 6 * edge_id + 4
            else:
                return 6 * edge_id + 5
        else:
            return 6 * edge_id + 3

    def rule_one(self, edge, children, node_edge_dict, roots):
        """
        Rule 1:
        Connect focal brick (indexed by edge.child) to its child bricks
        """
        focal_node = edge.child
        for child in children:
            assert node_edge_dict[child] != node_edge_dict[focal_node]
            self.from_to_set.add(
                (
                    self.up_vertex(node_edge_dict[child], "in"),
                    self.up_vertex(node_edge_dict[focal_node], "out"),
                )
            )
            self.from_to_set.add(
                (
                    self.down_vertex(node_edge_dict[focal_node], "in"),
                    self.down_vertex(node_edge_dict[child], "out"),
                )
            )
        # Connect focal brick to its parent brick
        if edge.parent not in roots and focal_node not in roots:
            assert node_edge_dict[focal_node] != node_edge_dict[edge.parent]
            self.from_to_set.add(
                (
                    self.up_vertex(node_edge_dict[focal_node], "in"),
                    self.up_vertex(node_edge_dict[edge.parent], "out"),
                )
            )
            self.from_to_set.add(
                (
                    self.down_vertex(node_edge_dict[edge.parent], "in"),
                    self.down_vertex(node_edge_dict[focal_node], "out"),
                )
            )

    def rule_two(self, edge, siblings, node_edge_dict):
        # Rule 2: Connect focal brick to its siblings
        if len(siblings) > 1:
            if len(siblings) > 1:
                for item in itertools.combinations(siblings, 2):
                    self.from_to_set.add(
                        (
                            self.up_vertex(node_edge_dict[item[0]], "in"),
                            self.down_vertex(node_edge_dict[item[1]], "out"),
                        )
                    )
                    self.from_to_set.add(
                        (
                            self.up_vertex(node_edge_dict[item[1]], "in"),
                            self.down_vertex(node_edge_dict[item[0]], "out"),
                        )
                    )

    def rule_three(self, edge, children, node_edge_dict, prev_edge_dict):
        """
        Rule 3: Connect focal brick to other bricks which share a child haplotype
        across a recombination
        """
        if edge.child in prev_edge_dict:
            # r,d of left parent to r of right parent
            self.from_to_set.add(
                (
                    self.right_vertex(prev_edge_dict[edge.child], "in"),
                    self.right_vertex(node_edge_dict[edge.child], "out"),
                )
            )
            self.from_to_set.add(
                (
                    self.down_vertex(prev_edge_dict[edge.child], "in"),
                    self.right_vertex(node_edge_dict[edge.child], "out"),
                )
            )
            # l,d of right parent to l of left parent
            self.from_to_set.add(
                (
                    self.left_vertex(node_edge_dict[edge.child], "in"),
                    self.left_vertex(prev_edge_dict[edge.child], "out"),
                )
            )
            self.from_to_set.add(
                (
                    self.down_vertex(node_edge_dict[edge.child], "in"),
                    self.left_vertex(prev_edge_dict[edge.child], "out"),
                )
            )

    def make_connections(self, edge, tree2, node_edge_dict, index, prev_edge_dict):
        roots = tree2.roots
        children = tree2.children(edge.child)
        siblings = tree2.children(edge.parent)
        self.rule_one(edge, children, node_edge_dict, roots)
        if self.use_rule_two:
            self.rule_two(edge, siblings, node_edge_dict)
        if index != 0:
            self.rule_three(edge, siblings, node_edge_dict, prev_edge_dict)

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
        self.G = None
        self.l_in = []
        self.l_out = []

        node_edge_dict = {}

        # Rule Zero
        # For unlabeled nodes: connect left to up, right to up (within a brick)
        for brick in self.unlabeled_bricks:
            self.from_to_set.add((6 * brick + 2, 6 * brick))
            self.from_to_set.add((6 * brick + 3, 6 * brick))

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

        # Create networkx graph
        df = pd.Datamuts_to_merge_dict = {
            "from": [cur_set[0] for cur_set in self.from_to_set],
            "to": [cur_set[1] for cur_set in self.from_to_set],
        }
        self.G = nx.from_pandas_edgelist(df, "from", "to", create_using=nx.DiGraph())
        # Total number of nodes should be less than (2 * number of labeled nodes) +
        # (4 * number of unlabeled nodes)
        # TODO: check why this isn't equal
        assert self.G.number_of_nodes() <= (2 * len(self.labeled_bricks)) + 4 * len(
            self.unlabeled_bricks
        ), (
            self.G.number_of_nodes(),
            (2 * len(self.labeled_bricks)),
            4 * len(self.unlabeled_bricks),
        )
        return self.G
