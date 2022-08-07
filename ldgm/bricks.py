"""
Contains code to create a "bricked" tree sequence, where edges with
have different sets of descendants are bifurcated.
"""
import tskit
from tqdm import tqdm
import json

from . import utility


class Bricks:
    def __init__(
        self,
        ts,
        recombination_freq_threshold=None,
        add_dummy_bricks=True,
        progress=True,
    ):
        self.ts = ts
        self.add_dummy_bricks = add_dummy_bricks
        if recombination_freq_threshold is None:
            rec_threshold = 0
        self.rec_threshold = rec_threshold
        self.progress = progress

    def bifurcate_edge(self, edge_child, interval, tables, current_edges):
        """
        Bifurcate a given edge if it straddles a breakpoint
        Returns an amended edge table and dictionary of edges in current tree
        """
        if edge_child in current_edges:
            overlap_edge = current_edges[edge_child]
            if interval.left != overlap_edge.left:
                tables.edges.add_row(
                    left=overlap_edge.left,
                    right=interval.left,
                    parent=overlap_edge.parent,
                    child=overlap_edge.child,
                )
                current_edges.pop(edge_child)
                current_edges[edge_child] = tskit.Edge(
                    left=interval.left,
                    right=overlap_edge.right,
                    parent=overlap_edge.parent,
                    child=overlap_edge.child,
                )
        return tables, current_edges

    def naive_split_edges(self, mode="leaf"):
        ts = self.ts

        if mode not in ["tree", "node", "leaf"]:
            raise ValueError("Unrecognised mode")

        tables = ts.dump_tables()
        tables.edges.clear()  # build up a new edge table

        trees = ts.trees()
        tree = next(trees)
        edge_diffs = ts.edge_diffs()
        _, edges_out, edges_in = next(edge_diffs)

        current_edges = {}  # keep track of edges in current tree

        # At first tree, add all edges in
        for edge in edges_in:
            current_edges[edge.child] = edge

        prev_tree = tree.copy()

        # Look at edges in and out for subsequent trees
        for tree, (interval, edges_out, edges_in) in tqdm(
            zip(trees, edge_diffs),
            desc="Brick tree sequence: iterate over trees",
            total=ts.num_trees - 1,
            disable=not self.progress,
        ):
            # Add edges coming out to new edge table
            for edge in edges_out:
                modified_edge = current_edges[edge.child]
                tables.edges.add_row(
                    left=modified_edge.left,
                    right=modified_edge.right,
                    parent=modified_edge.parent,
                    child=modified_edge.child,
                )
                current_edges.pop(edge.child)
            # Bifurcate edges based on edges_in
            for edge in edges_in:
                current_edges[edge.child] = edge
                if mode == "leaf":
                    right = edge.parent
                    left = prev_tree.parent(edge.child)
                    if (
                        tree.num_samples(edge.child) / ts.num_samples
                        > self.rec_threshold
                    ):
                        while right != left and right != -1 and left != -1:
                            tr = tree.get_time(right)
                            tl = prev_tree.get_time(left)
                            if tr < tl:
                                tables, current_edges = self.bifurcate_edge(
                                    right, interval, tables, current_edges
                                )
                                right = tree.parent(right)
                            elif tr > tl:
                                tables, current_edges = self.bifurcate_edge(
                                    left, interval, tables, current_edges
                                )
                                left = prev_tree.parent(left)
                            else:
                                tables, current_edges = self.bifurcate_edge(
                                    right, interval, tables, current_edges
                                )
                                right = tree.parent(right)
                elif mode == "node":
                    parent = edge.parent
                    while parent != tree.root:
                        tables, current_edges = self.bifurcate_edge(
                            parent, interval, tables, current_edges
                        )
                        parent = tree.parent(parent)

            if mode == "tree":
                cur_current_edges = current_edges.copy()
                for edge_child, edge in cur_current_edges.items():
                    if edge.left != interval[0]:
                        tables, current_edges = self.bifurcate_edge(
                            edge_child, interval, tables, current_edges
                        )
            prev_tree = tree.copy()

        # Any edges still in current_edges at the last tree are added
        for _, edge in current_edges.items():
            tables.edges.add_row(
                left=edge.left, right=edge.right, parent=edge.parent, child=edge.child
            )
        tables.sort()
        tables.provenances.add_row(
            record=json.dumps(utility.get_provenance_dict({"command": "brick_ts"}))
        )
        new_ts = tables.tree_sequence()
        if self.add_dummy_bricks:
            new_ts = utility.add_dummy_bricks(new_ts)
        return new_ts
