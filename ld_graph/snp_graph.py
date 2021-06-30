"""
Create a reduced SNP graph
"""
import collections

import networkx as nx
import numpy as np
from tqdm import tqdm

from . import utility


class SNP_Graph:
    def __init__(self, brick_graph, brick_ts):
        self.brick_graph = brick_graph
        self.brick_ts = brick_ts

        # Tree sequence must contain mutations
        if brick_ts.num_mutations == 0:
            raise ValueError("Tree sequence must contain mutations")

        bricks_to_muts = utility.get_mut_edges(brick_ts)

        id_to_muts = {}
        bricks_to_id = {}

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

        nodes = np.array(list(brick_graph.nodes()))
        self.l_out = nodes[nodes % 6 == 4]
        self.l_in = nodes[nodes % 6 == 5]

    def create_reduced_graph(self):
        desc = collections.defaultdict()
        R_desc = nx.Graph()
        R_desc.add_nodes_from(self.l_out)
        R_desc.add_nodes_from(self.l_in)
        for u in tqdm(self.l_out, total=len(self.l_out), desc="Graph reduction"):
            desc[u] = nx.descendants(self.brick_graph, u)
            for v in desc[u]:
                if v in self.l_in:
                    R_desc.add_edge(u, v)
        nodes = list(R_desc.nodes())
        equivalent_nodes = []
        vals = np.stack([self.l_out, self.l_out + 1]).transpose()
        for index, i in tqdm(enumerate(vals)):
            if i[0] in nodes and i[1] in nodes:
                equivalent_nodes.append(vals[index])
        R_identified = nx.quotient_graph(R_desc, equivalent_nodes)

        def node_data(list_of_lists):
            node_data = {}
            for _, i in enumerate(list_of_lists):
                node_data[frozenset({i[0], i[1]})] = int(i[0] / 6)
            return node_data

        R_identified_relabelled = nx.relabel_nodes(
            R_identified, node_data(equivalent_nodes), copy=True
        )
        return R_identified_relabelled, self.id_to_muts
