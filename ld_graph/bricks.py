"""
Create a "bricked" tree sequence
"""
import numpy as np
import tskit
from tqdm import tqdm


class Bricks:
    def __init__(
        self,
        ts,
    ):
        self.ts = ts

    def check_straddle_edge_adjust(
        self, tree_left, parent, child, left, right, i, i_child
    ):
        straddle_edge = np.where(
            np.logical_and(
                np.logical_and(
                    np.logical_and(parent == i, child == i_child), left < tree_left
                ),
                right > tree_left,
            )
        )[0]
        assert len(straddle_edge) <= 1
        # Adjust straddling edge
        return straddle_edge

    def get_edge_to_add(self, tables, straddle_edge, left_interval):
        right = tables.edges.right
        left = tables.edges.left
        parent = tables.edges.parent
        child = tables.edges.child
        prev_right = right[straddle_edge[0]]
        right[straddle_edge[0]] = left_interval
        return left, right, parent, child, prev_right

    def bifurcate_edge(self, tables, tree2, parent_node, child_node):
        parent = tables.edges.parent
        left = tables.edges.left
        right = tables.edges.right
        child = tables.edges.child
        straddle_edge = self.check_straddle_edge_adjust(
            tree2.interval.left,
            parent,
            child,
            left,
            right,
            parent_node,
            child_node,
        )
        if len(straddle_edge) > 0:
            (
                left,
                right,
                parent,
                child,
                prev_right,
            ) = self.get_edge_to_add(tables, straddle_edge, tree2.interval.left)
            tables.edges.clear()
            tables.edges.set_columns(left=left, right=right, parent=parent, child=child)
            # Add a new edge for nodes to be bifurcated
            tables.edges.add_row(
                left=tree2.interval.left,
                right=prev_right,
                parent=parent_node,
                child=child_node,
            )
            return tables

    def brick_ts(self):
        """
        Create bricks given a tree sequence. Bricks are an "augmented edge table" such that each
        edge is bifurcated when any of its descendants recombine.
        """
        tables = self.ts.dump_tables()
        tables.edges.clear()
        tree1 = None
        loop_sizes = []
        for index, (tree2, (interval, edges_out, edges_in)) in tqdm(
            enumerate(zip(self.ts.trees(), self.ts.edge_diffs())),
            total=self.ts.num_trees,
        ):

            # At the first tree, add all edges in
            if index == 0:
                for edge in edges_in:
                    tables.edges.add_row(edge.left, edge.right, edge.parent, edge.child)
            # If not the first tree, check for nodes which are both coming in and going out as children
            else:
                children_in = [edge.child for edge in edges_in]
                children_out = [edge.child for edge in edges_out]

                # Nodes both coming in and going out need to be split into two nodes
                overlap = np.array(edges_in)[
                    np.isin(np.array(children_in), np.array(children_out))
                ]

                for edge in edges_in:
                    tables.edges.add_row(edge.left, edge.right, edge.parent, edge.child)

                for u in overlap:

                    i_child = u.child
                    j_child = u.child
                    i = tree1.parent(i_child)
                    j = tree2.parent(j_child)
                    loop_size = 0
                    while i != j:

                        ti = tree1.get_time(i)
                        tj = tree2.get_time(j)
                        if ti < tj and tree1.root != i:  # walk up left tree
                            i_child = i
                            i = tree1.get_parent(i)
                            # if edge from i_child to i straddles both trees, bifurcate
                            self.bifurcate_edge(tables, tree2, i, i_child)
                        elif ti > tj and tree2.root != j:  # walk up right tree
                            j_child = j
                            j = tree2.get_parent(j)
                            # if edge from j_child to j straddles both trees, bifurcate
                            self.bifurcate_edge(tables, tree2, j, j_child)
                        else:
                            # When we've reached the MRCA, no longer need to resolve edges
                            break
            tree1 = tree2.copy()
        tables.sort()

        bricked_ts = tables.tree_sequence()
        return bricked_ts

    def ld_graph(self, ts):
        """
        Given a "bricked" ts, connect bricks according to three rules:
        1. Connect parent and child bricks
        2. Connect sibling bricks
        3. Connect bricks which have two different adjacent parents
        """
        from_to_set = set()
        node_edge_dict = {}

        for index, (tree2, (interval, edges_out, edges_in)) in tqdm(
            enumerate(zip(ts.trees(), ts.edge_diffs())), total=ts.num_trees
        ):
            prev_edge_dict = node_edge_dict.copy()

            for edge in edges_out:
                node_edge_dict.pop(edge.child)
            for edge in edges_in:
                node_edge_dict[edge.child] = edge.id

            for edge in edges_in:
                node = edge.child
                # Rule 1: connect parent brick to its child bricks
                for child in tree2.children(node):
                    from_to_set.add((node_edge_dict[child], node_edge_dict[node]))
                    # Rule 3: Connect parent bricks of a recombined child brick
                    if (prev_edge_dict[child] == node_edge_dict[child]) and (
                        tree1.parent(child) in prev_edge_dict
                        and tree2.parent(child) in node_edge_dict
                    ):
                        if (
                            prev_edge_dict[tree1.parent(child)]
                            != node_edge_dict[tree2.parent(child)]
                        ):
                            from_to_set.add(
                                (
                                    node_edge_dict[tree2.parent(child)],
                                    prev_edge_dict[tree1.parent(child)],
                                )
                            )

                while node != -1:
                    # Rule 1: Connect parent bricks to child bricks
                    if tree2.parent(node) != tree2.root and node != tree2.root:
                        from_to_set.add(
                            (node_edge_dict[node], node_edge_dict[tree2.parent(node)])
                        )
                    # Rule 2: Connect sibling bricks to one another
                    children = tree2.children(node)
                    if len(children) > 1:
                        for child in children:
                            for child_1 in children[1:]:
                                from_to_set.add(
                                    (node_edge_dict[child], node_edge_dict[child_1])
                                )

                    node = tree2.parent(node)
            tree1 = tree2.copy()

        # Create networkx graph
        df = pd.Datamuts_to_merge_dict = {
            "from": [cur_set[0] for cur_set in from_to_set],
            "to": [cur_set[1] for cur_set in from_to_set],
        }
        G = nx.from_pandas_edgelist(df, "from", "to")
        return G
