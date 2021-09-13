"""
Regularize a tree sequence
"""
import networkx as nx
import numpy as np
import tskit

from . import utility


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

    def frequency_regularize(
        self, site_frequency, path_threshold, add_dummy_bricks=True, **kwargs
    ):
        """
        Takes a tree sequence and removes sites which are below a given frequency and
        bricks which are below a frequency calculated as site_frequency *
        e^(-0.5 * path_threshold).

        :param TreeSequence tree_sequence: The input :class:`tskit.TreeSequence` to be
            regularized.
        :param float snp_frequency: The maximum frequency of sites to be retained in the
            tree sequence. All sites with frequency less than or equal to `frequency`
            will be removed from the tree sequence.
        :param float path_threshold: The path_threshold which will be used to create the
            brick graph.
        :param bool add_dummy_bricks If True, adds dummy bricks above leaf nodes.
        :param \\**kwargs: All further keyword arguments are passed to
            ``tskit.simplify``. The keyword argument `samples` cannot be passed.
        :return: A tree sequence regularized by removing nodes based on their frequency.
        :rtype: tskit.TreeSequence
        :return: A numpy array whose ``u``-th element is the ID of the node in the
            regularized tree sequence that corresponds to node ``u`` in the original
            tree sequence., or ``tskit.NULL`` (-1) if ``u`` is no longer present
            in the regularized tree sequence.
        :rtype: numpy array
        :return: A numpy array whose ``u``-th element is the ID of the site in the
            regularized tree sequence that corresponds to site ``u`` in the original
            tree sequence., or ``tskit.NULL`` (-1) if ``u`` is no longer present in
            the regularized tree sequence.
        :rtype: numpy array
        """
        if "samples" in kwargs:
            raise ValueError("Cannot pass 'samples' to simplify()")
        brick_frequency = site_frequency * np.exp(-0.5 * path_threshold)
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
                np.where(paths > (self.ts.num_samples * brick_frequency))[0],
            )
        )[0]
        # Now iterate over sites and remove those with less than snp_frequency
        delete_sites = []
        site_map = np.full(self.ts.num_sites, -1)
        site_increment = 0
        for tree in self.ts.trees():
            for site in tree.sites():
                oldest_mut = site.mutations[0]
                for mut in site.mutations:
                    if self.ts.node(mut.node).time > self.ts.node(oldest_mut.node).time:
                        oldest_mut = mut
                if (tree.num_samples(mut.node) / self.ts.num_samples) < site_frequency:
                    delete_sites.append(site.id)
                else:
                    site_map[site.id] = site_increment
                    site_increment += 1
        # Now we make a new tree sequence with only the keep_nodes
        # to use ts.simplify(), we need to make all nodes samples
        tables = self.ts.dump_tables()
        tables.nodes.set_columns(
            flags=np.ones(self.ts.num_nodes).astype("uint32"),
            time=self.ts.tables.nodes.time,
        )
        tables.delete_sites(delete_sites)
        all_samples_ts = tables.tree_sequence()
        pruned_ts, node_map = all_samples_ts.simplify(
            samples=keep_nodes, map_nodes=True, **kwargs
        )
        assert pruned_ts.num_nodes == len(
            keep_nodes
        )  # make sure we have right number of nodes
        # tables = pruned_ts.dump_tables()
        # tables.nodes.clear()
        # # Make leaves samples
        # targets = set()
        # for tree in pruned_ts.trees():
        #    for leaf in tree.leaves():
        #        if tree.parent(leaf) != -1:
        #            targets.add(leaf)
        # for target in targets:
        #    tables.nodes.add_row(flags=1, time=pruned_ts.node(target).time)

        # # Add all the nodes in
        # for node in pruned_ts.nodes():
        #    if node.id not in targets:
        #        tables.nodes.add_row(flags=0, time=node.time)
        # tables.sort()
        # pruned_ts = tables.tree_sequence()
        if add_dummy_bricks:
            pruned_ts = utility.add_dummy_bricks(pruned_ts, mode="leaves")
            # Adjust node mapping when dummy bricks are added
            node_map[node_map != -1] = node_map[node_map != -1] + pruned_ts.num_samples
        return pruned_ts, node_map, site_map
