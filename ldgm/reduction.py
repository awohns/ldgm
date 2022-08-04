"""
Create a reduced SNP graph
"""
import collections
import multiprocessing

import networkx as nx
import numpy as np
from tqdm import tqdm

from . import utility


def find_reach_set(params):
    out_node = params[0]
    brick_graph = params[1]
    threshold = params[2]
    bricks_to_muts = params[3]
    reach_star_set = []
    new_edges = []
    removed_brick_graph = brick_graph.subgraph(
        [
            node
            for node in brick_graph.nodes()
            if node
            not in [
                out_node - 4,
                out_node - 3,
                out_node - 2,
                out_node - 1,
                out_node + 1,
            ]
        ]
    )
    reach_set = nx.single_source_dijkstra_path_length(
        removed_brick_graph, out_node, cutoff=threshold, weight="weight"
    )
    for vertex, weight in reach_set.items():
        brick_haplo_id = vertex // 8
        # that it is either labeled brick or haplotype
        is_haplo = vertex % 8 == 6
        # If brick id is in labeled list and not a haplotype,
        # it's a labeled brick
        is_labeled_brick = (
            brick_haplo_id in list(bricks_to_muts.keys())
        ) and not is_haplo
        if is_labeled_brick:
            # Make sure the vertex in the reach set is before node
            if vertex % 8 == 0 or vertex % 8 == 2:
                # Check that the reach set does NOT contain after nodes
                # for that brick
                if (brick_haplo_id * 8 + 1) not in reach_set and (
                    brick_haplo_id * 8 + 3
                ) not in reach_set:
                    # Make connection from SNP to SNP or SNP to haplotype
                    reach_star_set.append(vertex)
                    new_edges.append(
                        (
                            bricks_to_muts[out_node // 8][0],
                            bricks_to_muts[brick_haplo_id][0],
                            weight,
                        )
                    )
                    new_edges.append(
                        (
                            bricks_to_muts[brick_haplo_id][0],
                            bricks_to_muts[out_node // 8][0],
                            weight,
                        )
                    )
        elif is_haplo:
            if vertex % 8 == 6 and (brick_haplo_id * 8 + 7) not in reach_set:
                # Use -brick_haplo_id - 1 to avoid haplotype and brick id
                # collision
                reach_star_set.append(vertex)
                new_edges.append(
                    (
                        bricks_to_muts[out_node // 8][0],
                        -brick_haplo_id - 1,
                        weight,
                    )
                )
    return out_node, reach_star_set, new_edges


class SNP_Graph:
    def __init__(
        self,
        brick_graph,
        brick_ts,
        threshold,
        num_processes=1,
        chunksize=100,
        progress=True,
    ):
        self.brick_graph = brick_graph
        self.brick_ts = brick_ts
        self.threshold = threshold
        self.num_processes = num_processes
        self.chunksize = chunksize
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
        l_out = nodes[nodes % 8 == 4]
        assert len(l_out) <= self.brick_ts.num_sites

        # Make the reduced graph and add in all the SNPs corresponding to labeled bricks
        R = nx.DiGraph()
        # NOTE: The first mutation on a brick (the lowest ID) is used as the node ID
        # in the LDGM
        R.add_nodes_from([self.bricks_to_muts[node // 8][0] for node in l_out])

        reach_star_sets = collections.defaultdict(list)
        l_out_zipped = zip(
            l_out,
            (self.brick_graph for _ in range(len(l_out))),
            (self.threshold for _ in range(len(l_out))),
            (self.bricks_to_muts for _ in range(len(l_out))),
        )

        if self.num_processes == 1:
            for params in tqdm(
                l_out_zipped,
                disable=not self.progress,
                desc="Reduce graph: iterate over out nodes",
            ):
                out_node, reach_star_set, new_edges = find_reach_set(params)
                for new_edge in new_edges:
                    R.add_edge(new_edge[0], new_edge[1], weight=new_edge[2])
                reach_star_sets[out_node] = reach_star_set

        else:

            with multiprocessing.Pool(processes=self.num_processes) as pool:

                for row in tqdm(
                    pool.imap_unordered(
                        find_reach_set, l_out_zipped, chunksize=self.chunksize
                    ),
                    total=len(l_out),
                    disable=not self.progress,
                    desc="Reduce graph: iterate over out nodes",
                ):
                    for new_edge in row[2]:
                        R.add_edge(new_edge[0], new_edge[1], weight=new_edge[2])
                    reach_star_sets[row[0]] = row[1]

        return R, self.mut_node, reach_star_sets
