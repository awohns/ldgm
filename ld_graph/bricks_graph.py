"""
Create a brick graph
"""
import numpy as np

import pandas as pd
import networkx as nx


class BrickGraph:
    def __init__(
        self,
        bricked_ts,
    ):
        self.bricked_ts = bricked_ts

    def make_brick_graph(self):
        """
        Given a "bricked" ts, connect bricks according to three rules:
        1. Connect parent and child bricks
        2. Connect sibling bricks
        3. Connect bricks which have two different adjacent parents
        """
        bricked_ts = self.bricked_ts
        from_to_set = set()
        node_edge_dict = {}

        for index, (tree2, (interval, edges_out, edges_in)) in enumerate(
            zip(bricked_ts.trees(), bricked_ts.edge_diffs())
        ):
            prev_edge_dict = node_edge_dict.copy()

            for edge in edges_out:
                node_edge_dict.pop(edge.child)
            for edge in edges_in:
                node_edge_dict[edge.child] = edge.id

            for edge in edges_in:
                # Rule 1: connect parent brick to its child bricks
                for child in tree2.children(edge.child):
                    if node_edge_dict[child] != node_edge_dict[edge.child]:
                        from_to_set.add(
                            (node_edge_dict[child], node_edge_dict[edge.child])
                        )
                if edge.parent != tree2.root and edge.child != tree2.root:
                    if node_edge_dict[edge.child] != node_edge_dict[edge.parent]:
                        from_to_set.add(
                            (
                                node_edge_dict[edge.child],
                                node_edge_dict[edge.parent],
                            )
                        )
                # Rule 2: Connect sibling bricks to one another
                children = tree2.children(edge.parent)
                if len(children) > 1:
                    for child in children:
                        for child_1 in children[1:]:
                            if node_edge_dict[child] != node_edge_dict[child_1]:
                                from_to_set.add(
                                    (node_edge_dict[child], node_edge_dict[child_1])
                                )
                # Rule 3: Connect bricks that share a child haplotype across a recombination
                if index != 0:
                    if edge.child in prev_edge_dict:
                        from_to_set.add(
                            (node_edge_dict[edge.child], prev_edge_dict[edge.child])
                        )

            tree1 = tree2.copy()

        # Create networkx graph
        df = pd.Datamuts_to_merge_dict = {
            "from": [cur_set[0] for cur_set in from_to_set],
            "to": [cur_set[1] for cur_set in from_to_set],
        }
        G = nx.from_pandas_edgelist(df, "from", "to")
        return G
