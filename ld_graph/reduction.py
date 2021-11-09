"""
Create a reduced SNP graph
"""
import networkx as nx
import numpy as np
from tqdm import tqdm

from . import utility


class SNP_Graph:
    def __init__(self, brick_graph, brick_ts, threshold, progress=True):
        self.brick_graph = brick_graph
        self.threshold = threshold
        self.progress = progress

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

        # Dictionary with keys = mutation id, values = graph node id
        self.mut_node = np.zeros(brick_ts.num_mutations)
        for _, val in self.bricks_to_muts.items():
            for v in val:
                self.mut_node[v] = val[0]

    def create_reduced_graph(self):
        nodes = np.array(list(self.brick_graph.nodes()))
        # l_in = nodes[nodes % 6 == 4]
        l_out = nodes[nodes % 6 == 5]

        R = nx.Graph()
        R.add_nodes_from([self.bricks_to_muts[node // 6][0] for node in l_out])
        # assert len(l_in) == len(l_out)

        for u in tqdm(
            l_out,
            desc="Reduce graph: iterate over out nodes",
            disable=not self.progress,
        ):
            length = nx.single_source_dijkstra_path_length(
                self.brick_graph, u, cutoff=self.threshold, weight="weight"
            )
            for key, value in length.items():
                if key % 6 == 4 and key // 6 != u // 6:
                    R.add_edge(
                        self.bricks_to_muts[u // 6][0],
                        self.bricks_to_muts[key // 6][0],
                        weight=value,
                    )

        return R, self.mut_node
