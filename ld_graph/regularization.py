"""
Regularize a tree sequence
"""
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
