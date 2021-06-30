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
    def __init__(
        self,
        bricked_ts,
    ):
        self.bricked_ts = bricked_ts

    # make an argument for in and out here, so we know how to split if it's labeled
    def up_vertex(self, edge_id, in_out):
        if edge_id in self.labeled_bricks:
            if in_out == "in":
                self.labeled_nodes.append(6 * edge_id + 4)
                return 6 * edge_id + 4
            else:
                self.labeled_nodes.append(6 * edge_id + 5)
                return 6 * edge_id + 5
        else:
            return 6 * edge_id

    def down_vertex(self, edge_id, in_out):
        if edge_id in self.labeled_bricks:
            if in_out == "in":
                self.labeled_nodes.append(6 * edge_id + 4)
                return 6 * edge_id + 4
            else:
                self.labeled_nodes.append(6 * edge_id + 5)
                return 6 * edge_id + 5
        else:
            return 6 * edge_id + 1

    def left_vertex(self, edge_id, in_out):
        if edge_id in self.labeled_bricks:
            if in_out == "in":
                self.labeled_nodes.append(6 * edge_id + 4)
                return 6 * edge_id + 4
            else:
                self.labeled_nodes.append(6 * edge_id + 5)
                return 6 * edge_id + 5
        else:
            return 6 * edge_id + 2

    def right_vertex(self, edge_id, in_out):
        if edge_id in self.labeled_bricks:
            if in_out == "in":
                self.labeled_nodes.append(6 * edge_id + 4)
                return 6 * edge_id + 4
            else:
                self.labeled_nodes.append(6 * edge_id + 5)
                return 6 * edge_id + 5
        else:
            return 6 * edge_id + 3

    def make_connections(
        self, edge, tree2, max_root, node_edge_dict, from_to_set, index, prev_edge_dict
    ):
        # Rule 1: connect parent brick to its child bricks
        for child in tree2.children(edge.child):
            if node_edge_dict[child] != node_edge_dict[edge.child]:
                from_to_set.add(
                    (
                        self.up_vertex(node_edge_dict[child], "in"),
                        self.up_vertex(node_edge_dict[edge.child], "out"),
                    )
                )
                from_to_set.add(
                    (
                        self.down_vertex(node_edge_dict[edge.child], "in"),
                        self.down_vertex(node_edge_dict[child], "out"),
                    )
                )
        if edge.parent != max_root and edge.child != max_root:
            if node_edge_dict[edge.child] != node_edge_dict[edge.parent]:
                from_to_set.add(
                    (
                        self.up_vertex(node_edge_dict[edge.child], "in"),
                        self.up_vertex(node_edge_dict[edge.parent], "out"),
                    )
                )
                from_to_set.add(
                    (
                        self.down_vertex(node_edge_dict[edge.parent], "in"),
                        self.down_vertex(node_edge_dict[edge.child], "out"),
                    )
                )
        # Rule 2: Connect sibling bricks to one another
        children = tree2.children(edge.parent)
        if len(children) > 1:
            if len(children) > 1:
                for item in itertools.combinations(children, 2):
                    from_to_set.add(
                        (
                            self.up_vertex(node_edge_dict[item[0]], "in"),
                            self.down_vertex(node_edge_dict[item[1]], "out"),
                        )
                    )
                    from_to_set.add(
                        (
                            self.up_vertex(node_edge_dict[item[1]], "in"),
                            self.down_vertex(node_edge_dict[item[0]], "out"),
                        )
                    )

        # Rule 3: Connect parent bricks sharing child haplotype across a recombination
        if index != 0:
            if edge.child in prev_edge_dict:
                # r,d of left parent to r of right parent
                from_to_set.add(
                    (
                        self.right_vertex(prev_edge_dict[edge.child], "in"),
                        self.right_vertex(node_edge_dict[edge.child], "out"),
                    )
                )
                from_to_set.add(
                    (
                        self.down_vertex(prev_edge_dict[edge.child], "in"),
                        self.right_vertex(node_edge_dict[edge.child], "out"),
                    )
                )
                # l,d of right parent to l of left parent
                from_to_set.add(
                    (
                        self.left_vertex(node_edge_dict[edge.child], "in"),
                        self.left_vertex(prev_edge_dict[edge.child], "out"),
                    )
                )
                from_to_set.add(
                    (
                        self.down_vertex(node_edge_dict[edge.child], "in"),
                        self.left_vertex(prev_edge_dict[edge.child], "out"),
                    )
                )
        return from_to_set

    def make_brick_graph(self):
        """
        Given a "bricked" ts, connect bricks according to three rules:
        1. Connect parent and child bricks
        2. Connect sibling bricks
        3. Connect bricks which have two different adjacent parents
        """
        bricks_to_muts = utility.get_mut_edges(self.bricked_ts)
        self.labeled_bricks = list(bricks_to_muts.keys())
        bricks = np.arange(0, self.bricked_ts.num_nodes)
        self.unlabeled_bricks = bricks[~np.isin(bricks, self.labeled_bricks)]
        self.labeled_nodes = []
        self.G = None
        self.l_in = []
        self.l_out = []

        from_to_set = set()
        node_edge_dict = {}

        times = self.bricked_ts.tables.nodes.time
        # Rule Zero
        # For unlabeled nodes: connect left to up, right to up (within a brick)
        for brick in self.unlabeled_bricks:
            from_to_set.add((6 * brick + 2, 6 * brick))
            from_to_set.add((6 * brick + 3, 6 * brick))

        for index, (tree2, (_, edges_out, edges_in)) in tqdm(
            enumerate(zip(self.bricked_ts.trees(), self.bricked_ts.edge_diffs())),
            desc="Brick graph: iterate over edges",
            total=self.bricked_ts.num_trees,
        ):
            prev_edge_dict = node_edge_dict.copy()
            roots = tree2.roots
            max_root = roots[np.argmax(times[roots])]

            for edge in edges_out:
                node_edge_dict.pop(edge.child)
            for edge in edges_in:
                node_edge_dict[edge.child] = edge.id

            for edge in edges_in:
                from_to_set = self.make_connections(
                    edge,
                    tree2,
                    max_root,
                    node_edge_dict,
                    from_to_set,
                    index,
                    prev_edge_dict,
                )

        # Create networkx graph
        df = pd.Datamuts_to_merge_dict = {
            "from": [cur_set[0] for cur_set in from_to_set],
            "to": [cur_set[1] for cur_set in from_to_set],
        }
        self.G = nx.from_pandas_edgelist(df, "from", "to", create_using=nx.DiGraph())
        self.labeled_nodes = np.unique(self.labeled_nodes)
        return self.G
