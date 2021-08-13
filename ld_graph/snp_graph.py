"""
Create a reduced SNP graph
"""
import networkx as nx
import numpy as np
from tqdm import tqdm

from . import utility


class SNP_Graph:
    def __init__(self, brick_graph, brick_ts, threshold):
        self.brick_graph = brick_graph
        self.threshold = threshold

        # Tree sequence must contain mutations
        if brick_ts.num_mutations == 0:
            raise ValueError("Tree sequence must contain mutations")

        # Dictionary with keys = brick ids, values = mutation ids
        self.bricks_to_muts = utility.get_mut_edges(brick_ts)

        # Dictionary with keys = Brick ordered id, value = mutation on brick
        id_to_muts = {}
        # Dictionary with keys = Brick id, values = mutation on brick
        bricks_to_id = {}

        for index, (brick, muts) in enumerate(self.bricks_to_muts.items()):
            id_to_muts[index] = muts
            bricks_to_id[brick] = index

    def create_reduced_graph(self):
        nodes = np.array(list(self.brick_graph.nodes()))
        l_in = nodes[nodes % 4 == 2]
        l_out = nodes[nodes % 4 == 3]

        path_weights = {}
        for u in tqdm(l_out):
            for v in l_in:
                try:
                    path = nx.dijkstra_path(self.brick_graph, u, v, weight="weight")
                    weight_sum = 0
                    for i in np.arange(len(path) - 1):
                        weight_sum += self.brick_graph[path[i]][path[i + 1]]["weight"]
                    path_weights[(u, v)] = weight_sum
                except nx.NetworkXNoPath:
                    continue

        # Now define the graph
        R = nx.Graph()
        R.add_nodes_from([self.bricks_to_muts[node // 4][0] for node in l_out])
        assert len(l_in) == len(l_out)

        for path, path_weight in path_weights.items():
            u = path[0]
            v = path[1]
            if path_weight < self.threshold:
                R.add_edge(
                    self.bricks_to_muts[u // 4][0], self.bricks_to_muts[v // 4][0]
                )
        return R
