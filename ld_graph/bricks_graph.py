"""
Create a brick graph
"""
import itertools

import networkx as nx
import numpy as np
from tqdm import tqdm

from . import utility


class BrickGraph:
    def __init__(self, bricked_ts, use_rule_three):
        self.bricked_ts = bricked_ts
        self.brick_graph = nx.Graph()
        self.use_rule_three = use_rule_three
        self.freqs = utility.get_brick_frequencies(self.bricked_ts)

    # make an argument for in and out here, so we know how to split if it's labeled
    def labeled_vertex(self, edge_id):
        if edge_id in self.labeled_bricks:
            return 2 * edge_id + 1
        else:
            return 2 * edge_id

    def rule_one(self, edge, children, node_edge_dict, roots):
        """
        Rule 1:
        Connect focal brick (indexed by edge.child) to its child bricks
        """
        focal_node = edge.child
        for child in children:
            assert node_edge_dict[child] != node_edge_dict[focal_node]
            # get parent node and child node
            child_brick = self.labeled_vertex(node_edge_dict[child])
            parent_brick = self.labeled_vertex(node_edge_dict[focal_node])
            self.brick_graph.add_edge(parent_brick, child_brick)

        # Connect focal brick to its parent brick
        if edge.parent not in roots and focal_node not in roots:
            assert node_edge_dict[focal_node] != node_edge_dict[edge.parent]
            child_brick = self.labeled_vertex(node_edge_dict[focal_node])
            parent_brick = self.labeled_vertex(node_edge_dict[edge.parent])
            self.brick_graph.add_edge(child_brick, parent_brick)

    def rule_two(self, edge, siblings, node_edge_dict):
        # Rule 2: Connect focal brick to its siblings
        if len(siblings) > 1:
            for pair in itertools.combinations(siblings, 2):
                left_brick = self.labeled_vertex(node_edge_dict[pair[0]])
                right_brick = self.labeled_vertex(node_edge_dict[pair[1]])
                self.brick_graph.add_edge(left_brick, right_brick)

    def rule_three(self, edge, siblings, node_edge_dict, prev_edge_dict):
        # Rule 3: Connect left and right parents which share a child
        if edge.child in prev_edge_dict:
            left_parent = self.labeled_vertex(prev_edge_dict[edge.child])
            right_parent = self.labeled_vertex(node_edge_dict[edge.child])
            self.brick_graph.add_edge(left_parent, right_parent)

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
        if index != 0 and self.use_rule_three is True:
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

        node_edge_dict = {}

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

        # Total number of nodes should be less than number of bricks * 2
        assert self.brick_graph.number_of_nodes() <= self.bricked_ts.num_edges * 2
        return self.brick_graph
