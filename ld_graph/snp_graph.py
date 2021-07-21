"""
Create a reduced SNP graph
"""
import networkx as nx
import numpy as np

from . import utility


class SNP_Graph:
    def __init__(self, brick_graph, brick_ts, identify_in_out=False):
        self.brick_graph = brick_graph
        self.identify_in_out = identify_in_out

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
        print("Condensing graph...")
        C = nx.condensation(self.brick_graph)
        print("Finding transitive closure...")
        H = nx.transitive_closure_dag(C, topo_order=None)

        nodes = np.array(list(self.brick_graph.nodes()))
        l_out = nodes[nodes % 6 == 4]
        l_in = nodes[nodes % 6 == 5]

        # Condensation creates graph with connected components (cc)
        # we want mapping of cc to original nodes
        cc_l_out = []
        cc_l_in = []
        mapping = C.graph["mapping"]
        for out_node in l_out:
            cc_l_out.append(mapping[out_node])
        for out_node in l_in:
            cc_l_in.append(mapping[out_node])
        cc_l_out = np.array(cc_l_out)
        cc_l_in = np.array(cc_l_in)

        # Also correct the id_to_muts with mapping
        mapped_id_to_muts = {}
        for node_id in self.id_to_muts:
            mapped_id_to_muts[node_id] = mapping[node_id * 6]

        R = H.subgraph(np.concatenate([cc_l_in, cc_l_out]))
        if self.identify_in_out:
            nodes = list(R.nodes())
            equivalent_nodes = []
            for i in l_out:
                equivalent_nodes.append((mapping[i], mapping[i + 1]))
            print("Finding quotient graph...")
            R_identified = nx.quotient_graph(R, equivalent_nodes)

            def node_data(list_of_lists):
                node_data = {}
                for _, i in enumerate(list_of_lists):
                    node_data[frozenset({i[0], i[1]})] = int(i[0] / 6)
                return node_data

            R_identified_relabelled = nx.relabel_nodes(
                R_identified, node_data(equivalent_nodes), copy=True
            )
            return R_identified_relabelled, mapped_id_to_muts
        else:
            return R, mapped_id_to_muts
