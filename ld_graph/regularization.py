"""
Regularize a tree sequence
"""
import networkx as nx
import numpy as np
import tskit


class Regularize:
    def __init__(
        self,
        ts,
    ):
        self.ts = ts

    def time_regularize(self, time):
        """
        Returns a tree sequence in which all the topology and mutation
        information after the specified time is remove. That is, the
        samples in the returned tree sequence will be the lineages
        extant at the specified time.
        """
        tables = self.ts.dump_tables()
        t = tables.nodes.time

        tables.edges.clear()
        samples = []
        edge_map = {}
        for edge in self.ts.edges():
            if t[edge.child] > time:
                tables.edges.add_row(edge.left, edge.right, edge.parent, edge.child)
            elif t[edge.child] <= time and t[edge.parent] > time:
                key = (edge.child, edge.parent)
                if key in edge_map:
                    # If the same parent/child combination exists, then we should map
                    # to the same ancestor node.
                    u = edge_map[key]
                else:
                    u = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=time)
                    samples.append(u)
                    edge_map[key] = u
                tables.edges.add_row(edge.left, edge.right, edge.parent, u)
        tables.sort()
        tables.simplify(samples)
        return tables.tree_sequence()


def frequency_regularize(self, frequency, **kwargs):
    """
    Takes a tree sequence and removes nodes which are less than the given frequency.

    :param TreeSequence tree_sequence: The input :class:`tskit.TreeSequence` to be
        regularized.
    :param float frequency: The maximum frequency of nodes to be retained in the tree
        sequence. All nodes with frequency less than or equal to `frequency` will
        be removed from the tree sequence.
    :param \\**kwargs: All further keyword arguments are passed to the ``tskit.simplify``
        The keyword argument `samples` cannot be passed.
    :return: A tree sequence regularized by removing nodes based on their frequency.
    :rtype: tskit.TreeSequence
    """
    if "samples" in kwargs:
        raise ValueError("Cannot pass 'samples' to simplify()")
    edges = np.array(
        [
            (parent, child)
            for child, parent in zip(
                self.ts.tables.edges.child, self.ts.tables.edges.parent
            )
        ]
    )
    G = nx.DiGraph()
    G.add_edges_from(edges)
    G_rev = G.reverse()
    paths = np.zeros(G_rev.number_of_nodes())
    for node in nx.topological_sort(G_rev):
        if G.out_degree(node) == 0:
            paths[node] = 1
        else:
            for edge in G.out_edges(node):
                paths[node] += paths[edge[1]]
    keep_nodes = np.where(
        np.isin(
            np.arange(0, self.ts.num_nodes),
            np.where(paths > (self.ts.num_samples * frequency))[0],
        )
    )[0]
    # Now we make a new tree sequence with only the keep_nodes
    # to use ts.simplify(), we need to make all nodes samples
    tables = self.ts.dump_tables()
    tables.nodes.set_columns(
        flags=np.ones(self.ts.num_nodes).astype("uint32"),
        time=self.ts.tables.nodes.time,
    )
    all_samples_ts = tables.tree_sequence()
    pruned_ts = all_samples_ts.simplify(samples=keep_nodes, **kwargs)
    assert pruned_ts.num_nodes == len(
        keep_nodes
    )  # make sure we have right number of nodes
    return pruned_ts
