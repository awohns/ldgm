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
        brick_mut_list = collections.defaultdict(list)
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
        brick_to_muts = {brick: mut for mut, brick in muts_to_brick.items()}
        for mut, brick in muts_to_brick.items():
            brick_mut_list[brick].append(mut)
        self.labelled_nodes = list(brick_to_muts.keys())
        self.unlabelled_nodes = list(
            np.arange(0, brick_ts.num_edges)[
                ~np.isin(np.arange(0, brick_ts.num_edges), list(brick_to_muts.keys()))
            ],
        )
        self.brick_to_muts = brick_to_muts

    def create_reduced_graph(self):
        H = self.brick_graph.subgraph(self.labelled_nodes)
        reduced_graph = nx.Graph(self.brick_graph.subgraph(self.unlabelled_nodes))
        for c in nx.connected_components(H):
            b = nx.node_boundary(self.brick_graph, c, self.labelled_nodes)
            for i in b:
                for j in b:
                    if i != j:
                        reduced_graph.add_edge(i, j)
        reduced_graph = nx.relabel_nodes(reduced_graph, self.brick_to_muts, copy=True)
        return reduced_graph
