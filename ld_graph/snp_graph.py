"""
Create a reduced SNP graph
"""
import networkx as nx
import numpy as np

from . import utility


class SNP_Graph:
    def __init__(self, brick_graph, brick_ts):
        self.brick_graph = brick_graph

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

    def create_reduced_graph(self):
        C = nx.condensation(self.brick_graph)
        # H = nx.transitive_closure_dag(C, topo_order=None)

        nodes = np.array(list(self.brick_graph.nodes()))
        l_in = nodes[nodes % 6 == 4]
        l_out = nodes[nodes % 6 == 5]

        # Condensation creates graph with connected components (cc)
        # we want mapping of cc to original nodes
        cc_l_in = []
        cc_l_out = []
        mapping = C.graph["mapping"]
        for out_node in l_in:
            cc_l_in.append(mapping[out_node])
        for out_node in l_out:
            cc_l_out.append(mapping[out_node])

        # Get reverse mapping: keys are condensed graph nodes
        # values are brick graph nodes
        reverse_mapping = {value: key for key, value in mapping.items()}

        # Reverse map condensed graph node id to *first* mutation on associated brick
        cc_to_mut = {}
        for in_node, out_node in zip(cc_l_in, cc_l_out):
            cc_to_mut[in_node] = self.bricks_to_muts[reverse_mapping[in_node] // 6][0]
            cc_to_mut[out_node] = self.bricks_to_muts[reverse_mapping[out_node] // 6][0]

        # Find descendants of out nodes
        R = nx.Graph()
        R.add_nodes_from([cc_to_mut[node] for node in cc_l_out])
        cc_l_in = set(cc_l_in)
        assert len(cc_l_in) == len(cc_l_out)

        print("Finding descendents of each node..")
        for u in cc_l_out:
            for v in nx.descendants(C, u):
                if v in cc_l_in:
                    R.add_edge(cc_to_mut[u], cc_to_mut[v])

        return R
