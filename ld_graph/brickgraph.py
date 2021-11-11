"""
Create a brick graph
"""
import itertools

import networkx as nx
import numpy as np
from tqdm import tqdm

from . import utility


class BrickGraph:
    """
    Each brick has six nodes:
    0: up before
    1: up after
    2: down before
    3: down after
    4: in
    5: out
    """

    def __init__(self, bricked_ts, threshold, progress=True):
        self.bricked_ts = bricked_ts
        self.threshold = threshold
        self.brick_graph = nx.DiGraph()
        self.freqs = utility.get_brick_frequencies(self.bricked_ts)
        self.progress = progress

    def find_odds(self, brick):
        return self.freqs[brick] / (1 - self.freqs[brick])

    def log_odds(self, odds):
        # TODO check why for numerical reasons we do this
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

    def vertex(self, edge_id, down, after, out):
        vertex_id = 6 * edge_id
        if edge_id in self.labeled_bricks:
            if out:
                vertex_id += 1
            return vertex_id + 4
        if down:
            vertex_id += 2
        if after:
            vertex_id += 1
        return vertex_id

    def connect_vertices(
        self,
        edge_id_a,
        edge_id_b,
        child=None,
        down_a=False,
        down_b=False,
        after_a=False,
        after_b=False,
        reverse_odds=False,
        combine_odds="division",
    ):
        """
        Up and before are default for both verticies
        """
        vertex_a = self.vertex(edge_id_a, out=True, after=after_a, down=down_a)
        vertex_b = self.vertex(edge_id_b, out=False, after=after_b, down=down_b)
        odds_a = self.find_odds(edge_id_a)
        odds_b = self.find_odds(edge_id_b)
        if reverse_odds:
            odds_c = odds_a
            odds_a = odds_b
            odds_b = odds_c
        if combine_odds == "division":
            weight = self.log_odds(odds_a / odds_b)
        elif combine_odds == "multiply":
            weight = self.log_odds(odds_a * odds_b)
        elif combine_odds == "self":
            weight = self.log_odds(1)
        else:
            raise ValueError("Incorrect combine_odds method")
        self.add_edge_threshold(vertex_a, vertex_b, weight)

    def rule_one(self, edge, children, roots):
        """
        Rule 1:
        Connect focal brick (indexed by edge.child) to its child bricks
        """
        focal_node = edge.child
        for child in children:
            assert self.node_edge_dict[child] != self.node_edge_dict[focal_node]
            # Up before of child to up before of parent
            self.connect_vertices(
                self.node_edge_dict[child], self.node_edge_dict[focal_node]
            )

            # Up after of child to up after of parent
            self.connect_vertices(
                self.node_edge_dict[child],
                self.node_edge_dict[focal_node],
                after_a=True,
                after_b=True,
            )

            # Down before of parent to down before of child
            self.connect_vertices(
                self.node_edge_dict[focal_node],
                self.node_edge_dict[child],
                down_a=True,
                down_b=True,
                reverse_odds=True,
            )

            # Down after of parent to down after of child
            self.connect_vertices(
                self.node_edge_dict[focal_node],
                self.node_edge_dict[child],
                after_a=True,
                after_b=True,
                down_a=True,
                down_b=True,
                reverse_odds=True,
            )

        # Connect focal brick to its parent brick, making same connections as above
        if edge.parent not in roots and focal_node not in roots:
            assert self.node_edge_dict[focal_node] != self.node_edge_dict[edge.parent]
            self.connect_vertices(
                self.node_edge_dict[focal_node], self.node_edge_dict[edge.parent]
            )
            self.connect_vertices(
                self.node_edge_dict[focal_node],
                self.node_edge_dict[edge.parent],
                after_a=True,
                after_b=True,
            )
            self.connect_vertices(
                self.node_edge_dict[edge.parent],
                self.node_edge_dict[focal_node],
                down_a=True,
                down_b=True,
                reverse_odds=True,
            )
            self.connect_vertices(
                self.node_edge_dict[edge.parent],
                self.node_edge_dict[focal_node],
                after_a=True,
                after_b=True,
                down_a=True,
                down_b=True,
                reverse_odds=True,
            )

    def rule_two(self, edge, siblings):
        # Rule 2: Connect focal brick to its siblings
        if len(siblings) > 1:
            for pair in itertools.combinations(siblings, 2):
                # Up before left to down before right
                self.connect_vertices(
                    self.node_edge_dict[pair[0]],
                    self.node_edge_dict[pair[1]],
                    down_b=True,
                    combine_odds="multiply",
                )
                # Up before right to down before left
                self.connect_vertices(
                    self.node_edge_dict[pair[1]],
                    self.node_edge_dict[pair[0]],
                    down_b=True,
                    combine_odds="multiply",
                )
                # Up after left to down after right
                self.connect_vertices(
                    self.node_edge_dict[pair[0]],
                    self.node_edge_dict[pair[1]],
                    after_a=True,
                    after_b=True,
                    down_b=True,
                    combine_odds="multiply",
                )
                # Up after right to down after left
                self.connect_vertices(
                    self.node_edge_dict[pair[1]],
                    self.node_edge_dict[pair[0]],
                    after_a=True,
                    after_b=True,
                    down_b=True,
                    combine_odds="multiply",
                )

    def rule_three(self, edges, child):
        """
        Connect left before down to right after up and right before down to left after up
        """
        if len(edges) > 1:
            for pair in itertools.combinations(edges, 2):
                # Down before left to up after right
                self.connect_vertices(
                    pair[0],
                    pair[1],
                    child=child,
                    down_a=True,
                    after_b=True,
                    combine_odds="self",
                )
                # Down before right to up after left
                self.connect_vertices(
                    pair[1],
                    pair[0],
                    child=child,
                    down_a=True,
                    after_a=True,
                    combine_odds="self",
                )

    def make_connections(self, edge, tree2, index):
        roots = tree2.roots
        children = tree2.children(edge.child)
        siblings = tree2.children(edge.parent)
        self.rule_one(edge, children, roots)
        if self.threshold is not None:
            if self.log_odds(self.find_odds(edge.id) ** 2) < self.threshold:
                self.rule_two(edge, siblings)
        else:
            self.rule_two(edge, siblings)

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

        self.node_edge_dict = {}

        # Rule Zero
        # For unlabeled nodes: connect down before to up after (within a brick)
        for brick in self.unlabeled_bricks:
            self.connect_vertices(
                brick,
                brick,
                down_a=True,
                after_b=True,
                combine_odds="self",
            )

        for index, (tree2, (_, edges_out, edges_in)) in tqdm(
            enumerate(zip(self.bricked_ts.trees(), self.bricked_ts.edge_diffs())),
            desc="Brick graph: iterate over trees",
            total=self.bricked_ts.num_trees,
            disable=not self.progress,
        ):
            self.prev_edge_dict = self.node_edge_dict.copy()

            for edge in edges_out:
                self.node_edge_dict.pop(edge.child)
            for edge in edges_in:
                self.node_edge_dict[edge.child] = edge.id
            for edge in edges_in:
                self.make_connections(edge, tree2, index)

        return self.brick_graph
