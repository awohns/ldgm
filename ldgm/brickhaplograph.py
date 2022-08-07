"""
Create a brick graph
"""
import itertools

import networkx as nx
import numpy as np
from tqdm import tqdm

from . import utility


class BrickHaploGraph:
    """
    Each brick has six or four vertices, depending on whether it is
    labeled or unlabeled, respectively.
    brick_id + 0: up before
    brick_id + 1: up after
    brick_id + 2: down before
    brick_id + 3: down after
    brick_id + 4: out
    brick_id + 5: uturn

    Haplotypes have two vertices. node_id refers to the ID of the
    node in the input tree sequence:
    node_id + 6: haplo before
    node_id + 7: haplo after

    To find vertices in the brick-haplotype for a brick or haplotype:
    brick/haplo_index * 8 + desired_node_type_id
    """

    def __init__(
        self, bricked_ts, edge_weight_threshold, make_sibs=False, progress=True
    ):
        self.bricked_ts = bricked_ts
        self.edge_weight_threshold = edge_weight_threshold
        self.brick_graph = nx.DiGraph()
        self.freqs = utility.get_brick_frequencies(self.bricked_ts)
        self.progress = progress
        self.make_sibs = make_sibs

    def find_odds(self, brick):
        if self.freqs[brick] == 1:
            raise ZeroDivisionError("Cannot have a brick with frequency 1")
        if self.freqs[brick] == 0:
            raise ValueError("Cannot have brick with frequency 0")
        odds = self.freqs[brick] / (1 - self.freqs[brick])
        assert odds != 0, (odds, self.freqs[brick])
        return odds

    def log_odds(self, odds):
        if odds == 0:
            raise ValueError("Cannot have odds of zero")
        if odds != 1:
            return np.log(odds) * -1
        else:
            return np.log(odds)

    def add_edge_threshold(self, from_node, to_node, weight):
        if self.edge_weight_threshold is not None:
            if weight < self.edge_weight_threshold:
                self.brick_graph.add_edge(from_node, to_node, weight=weight)
        else:
            self.brick_graph.add_edge(from_node, to_node, weight=weight)

    def vertex(self, identifier, down, after, out, uturn, haplo):
        """
        Function to convert haplo/brick ID to desired vertex id.
        identifier is the haplotype (node) or brick id from the tree sequence.
        Down is a bool, if true, refer to the down vertex.
        After is a bool, referring to after vertex
        Out indicates out vertex of labeled brick
        Haplo indicates if identifier refers to a haplotype
        """
        vertex_id = 8 * identifier
        if out:
            if identifier in self.labeled_bricks:
                return vertex_id + 4
            else:
                raise ValueError("vertex not a labeled brick")
        if uturn:
            return vertex_id + 5
        if down:
            vertex_id += 2
            if haplo:
                raise ValueError("Haplos cannot be down")
        if after:
            vertex_id += 1
        if haplo:
            vertex_id += 6
        return vertex_id

    def connect_vertices(
        self,
        id_a,
        id_b,
        child=None,
        out=False,
        down_a=False,
        down_b=False,
        after_a=False,
        after_b=False,
        uturn_a=False,
        uturn_b=False,
        haplo=False,
        reverse_odds=False,
        combine_odds="rule_one",
    ):
        """
        Up and before are default for both verticies
        """
        # Both vertices cannot be uturns
        assert ~(uturn_a is True and uturn_b is True)
        vertex_a = self.vertex(
            id_a, out=out, after=after_a, down=down_a, uturn=uturn_a, haplo=False
        )
        vertex_b = self.vertex(
            id_b, out=False, after=after_b, down=down_b, uturn=uturn_b, haplo=haplo
        )
        odds_a = self.find_odds(id_a)
        odds_b = self.find_odds(id_b)
        if combine_odds == "rule_one":
            if odds_b == 0:
                raise ValueError("odds of sink brick are 0")
            weight = abs(self.log_odds(odds_a / odds_b))
        elif combine_odds == "rule_two":
            weight = abs(self.log_odds(odds_a * odds_b))
        elif combine_odds == "haplo":
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

        def do_rule_one(parent, child):
            labeled_parent = self.node_edge_dict[parent] in self.labeled_bricks
            labeled_child = self.node_edge_dict[child] in self.labeled_bricks
            assert self.node_edge_dict[child] != self.node_edge_dict[parent]

            # Up after of child to up after of parent
            self.connect_vertices(
                self.node_edge_dict[child],
                self.node_edge_dict[parent],
                after_a=True,
                after_b=True,
            )

            # Down after of parent to down after of child
            self.connect_vertices(
                self.node_edge_dict[parent],
                self.node_edge_dict[child],
                after_a=True,
                after_b=True,
                down_a=True,
                down_b=True,
            )

            # Up before of child to EITHER up before of parent (if child is unlabeled)
            # or up after of parent (if child is labeled)
            self.connect_vertices(
                self.node_edge_dict[child],
                self.node_edge_dict[parent],
                after_b=labeled_child,
            )

            # Down before of parent to EITHER down before of child
            # (if parent is unlabeled) or down after of child (if parent is labeled)
            self.connect_vertices(
                self.node_edge_dict[parent],
                self.node_edge_dict[child],
                down_a=True,
                down_b=True,
                after_b=labeled_parent,
            )

            # If parent brick is labeled, make out and uturn connections
            if labeled_parent:
                # Out of parent to down before of child
                self.connect_vertices(
                    self.node_edge_dict[parent],
                    self.node_edge_dict[child],
                    out=True,
                    down_b=True,
                )
                if self.make_sibs:
                    # uturn of labeled parent to down before of ANY child
                    self.connect_vertices(
                        self.node_edge_dict[parent],
                        self.node_edge_dict[child],
                        uturn_a=True,
                        down_b=True,
                    )
            if labeled_child:
                # out of child to up before of parent
                self.connect_vertices(
                    self.node_edge_dict[child],
                    self.node_edge_dict[parent],
                    out=True,
                )
                if labeled_parent:
                    # out of labeled child to uturn of labeled parent
                    if self.make_sibs:
                        self.connect_vertices(
                            self.node_edge_dict[child],
                            self.node_edge_dict[parent],
                            out=True,
                            uturn_b=True,
                        )
            elif labeled_parent:
                # up before of UNLABELED CHILD to uturn of labeled parent
                if self.make_sibs:
                    self.connect_vertices(
                        self.node_edge_dict[child],
                        self.node_edge_dict[parent],
                        uturn_b=True,
                    )

        for child in children:
            do_rule_one(focal_node, child)

        # Connect focal brick to its parent brick, making same connections as above
        if edge.parent not in roots and focal_node not in roots:
            assert self.node_edge_dict[focal_node] != self.node_edge_dict[edge.parent]
            do_rule_one(edge.parent, focal_node)

    def rule_two(self, edge, siblings):
        # Rule 2: Connect focal brick to its siblings
        if len(siblings) > 1:
            for pair_combo in itertools.combinations(siblings, 2):
                for pair in (
                    (pair_combo[0], pair_combo[1]),
                    (pair_combo[1], pair_combo[0]),
                ):
                    assert pair[0] != pair[1]

                    labeled_0 = self.node_edge_dict[pair[0]] in self.labeled_bricks
                    # Up after left to down after right
                    self.connect_vertices(
                        self.node_edge_dict[pair[0]],
                        self.node_edge_dict[pair[1]],
                        after_a=True,
                        after_b=True,
                        down_b=True,
                        combine_odds="rule_two",
                    )

                    # If left is labeled, left up before to right down after.
                    # If left is not labeled, we do left up before to right down before.
                    self.connect_vertices(
                        self.node_edge_dict[pair[0]],
                        self.node_edge_dict[pair[1]],
                        down_b=True,
                        after_b=labeled_0,
                        combine_odds="rule_two",
                    )

                    # If brick 0 in sibling pair is labeled, make out connections
                    if labeled_0:
                        # Out to down before
                        self.connect_vertices(
                            self.node_edge_dict[pair[0]],
                            self.node_edge_dict[pair[1]],
                            out=True,
                            down_b=True,
                            combine_odds="rule_two",
                        )

    def make_connections(self, edge, tree2, index):
        roots = tree2.roots
        children = tree2.children(edge.child)
        siblings = tree2.children(edge.parent)
        self.rule_one(edge, children, roots)
        if self.make_sibs:
            if self.edge_weight_threshold is not None:
                if (
                    self.log_odds(self.find_odds(edge.id) ** 2)
                    < self.edge_weight_threshold
                ):
                    self.rule_two(edge, siblings)
            else:
                self.rule_two(edge, siblings)

    def make_brick_graph(self):
        """
        Given a "bricked" ts, connect bricks according to three rules:
        0. Connect each brick to its child haplotype
        1. Connect parent and child bricks
        2. Connect sibling bricks
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
        # Connect each brick to its child haplotype
        for brick in self.bricked_ts.edges():
            labeled_brick = brick.id in self.labeled_bricks
            # If brick is labeled: down before to after of child haplotype
            # Else if brick is unlabeled, down before of brick to before
            # of child haplotype
            self.connect_vertices(
                brick.id,
                brick.child,
                down_a=True,
                after_b=labeled_brick,
                haplo=True,
                combine_odds="haplo",
            )
            # Down after of brick to after of child haplotype
            self.connect_vertices(
                brick.id,
                brick.child,
                down_a=True,
                after_a=True,
                after_b=True,
                haplo=True,
                combine_odds="haplo",
            )
            if labeled_brick:
                # Out of brick to before of child haplotype
                self.connect_vertices(
                    brick.id,
                    brick.child,
                    out=True,
                    haplo=True,
                    combine_odds="haplo",
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

        # Sanity checks
        # No connections from a haplotype vertex pointing outwards
        # TODO check vertices 5 and 6 never start an edge
        # Vertex 7 is never used
        # Out is never second position in an edge
        return self.brick_graph
