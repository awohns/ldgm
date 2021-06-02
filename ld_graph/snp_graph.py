"""
Create a reduced SNP graph
"""
import numpy as np

import collections
import networkx as nx


class SNP_Graph:
    def __init__(self, brick_graph, brick_ts):
        self.brick_graph = brick_graph
        self.brick_ts = brick_ts

        # Tree sequence must contain mutations
        if brick_ts.num_mutations == 0:
            raise ValueError("Tree sequence must contain mutations")

        muts_to_brick = {}
        node_edge_dict = {}
        for tree, (interval, edges_out, edges_in) in zip(
            brick_ts.trees(), brick_ts.edge_diffs()
        ):
            for edge in edges_out:
                node_edge_dict.pop(edge.child)
            for edge in edges_in:
                node_edge_dict[edge.child] = edge.id
            for site in tree.sites():
                for mut in site.mutations:
                    node = mut.node
                    muts_to_brick[mut.id] = node_edge_dict[node]

        bricks_to_muts = collections.defaultdict(list)
        id_to_muts = {}
        bricks_to_id = {}

        for mut, brick in muts_to_brick.items():
            bricks_to_muts[brick].append(mut)

        for index, (brick, muts) in enumerate(bricks_to_muts.items()):
            id_to_muts[index] = muts
            bricks_to_id[brick] = index

        self.labelled_nodes = list(bricks_to_muts.keys())
        self.unlabelled_nodes = list(
            np.arange(0, brick_ts.num_edges)[
                ~np.isin(np.arange(0, brick_ts.num_edges), self.labelled_nodes)
            ],
        )
        self.bricks_to_muts = bricks_to_muts
        self.bricks_to_id = bricks_to_id
        self.id_to_muts = id_to_muts

    def create_reduced_graph(self):
        H = self.brick_graph.subgraph(self.unlabelled_nodes)
        reduced_graph = nx.Graph(self.brick_graph.subgraph(self.labelled_nodes))
        for c in nx.connected_components(H):
            b = nx.node_boundary(self.brick_graph, c, self.labelled_nodes)
            for i in b:
                for j in b:
                    if i != j:
                        reduced_graph.add_edge(i, j)
        reduced_graph = nx.relabel_nodes(reduced_graph, self.bricks_to_id, copy=True)
        return reduced_graph, self.id_to_muts
