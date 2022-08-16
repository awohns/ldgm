"""
Functions to produce reduced graph and intermediate steps
"""
import networkx as nx
from tqdm import tqdm

from . import brickhaplograph
from . import bricks
from . import reduction
from . import utility


def brick_ts(ts, recombination_freq_threshold=None, progress=True):
    """
    Take an input tree sequence and bifurcate edges to create "bricks" that
    have the same set of descendants at every position.

    :param tskit.TreeSequence ts: The input :class:`tskit.TreeSequence`, which has
        not had `ldgm.brick_ts()` run on it.
    :param float recombination_threshold: The minimum frequency above which bricks
        should be created by "bifurcating" edges. For example, if
        `recombination_freq_threshold` is 0.01, an edge with a different number of
        descendant samples to the left and right of a recombination event (marginal
        tree change) is bifurcated to create two bricks *if* the frequency of
        the edge's child node is >1% in the right hand marginal tree.
        If None, all edges with differing numbers of descendants to the left
        and right of a recombination event are bifurcated. Default: None
    :param bool progress: Whether to display a progress bar. Default: False
    :return: A bricked version of the input tree sequence
    :rtype: tskit.TreeSequence
    """
    assert utility.check_bricked(ts) is False
    brick = bricks.Bricks(
        ts,
        recombination_freq_threshold=recombination_freq_threshold,
        progress=progress,
    )
    bricked = brick.naive_split_edges()
    return bricked


def brick_haplo_graph(
    bricked_ts, edge_weight_threshold=None, make_sibs=True, progress=True
):
    """
    Takes a "bricked" tree sequence and creates a "brick-haplotype graph". This
    graph consists of vertices for labeled bricks (which have a mutation),
    unlabeled bricks (without a mutation) and haplotypes. For each labeled brick,
    we define a cluster of six vertices: Up-Before, Down-Before, Up-After, Down-After,
    U-Turn, and Out. For each unlabeled brick, we define a cluster of
    four vertices: Up-Before, Down-Before, Up-After, and Down-After. Haplotypes have
    a cluster of two vertices: Before and After. The union of these clusters form
    the vertex set of the brick-haplotype graph. Edges in the brick-haplotype graph
    are used to find the conditional dependencies between SNPs in  in the
    `ldgm.reduce_graph()` function

    Vertex labels determined as the ID of a brick (edge) in the input brick tree sequence
    plus an integer, which determines which the vertex's category:

    labeled bricks:
    brick_id + 0
    brick_id + 1
    brick_id + 2
    brick_id + 3
    brick_id + 4

    labeled bricks:
    brick_id + 0
    brick_id + 1
    brick_id + 2
    brick_id + 3
    brick_id + 4


    :param tskit.TreeSequence bricked_ts: The input :class:`tskit.TreeSequence`, which
        has been bricked (i.e. ``ldgm.brick_ts()`` has been run on it).
    :param float edge_weight_threshold: The
    :param bool make_sibs: Whether to run rule two, connecting siblings.
        Default: True
    :param bool progress: Whether to display a progress bar. Default: False
    :return: A ``Networkx`` graph encoding conditional dependence between
        labeled bricks.
    :rtype: networkx.DiGraph
    """
    assert utility.check_bricked(bricked_ts) is True
    brick_grapher = brickhaplograph.BrickHaploGraph(
        bricked_ts, edge_weight_threshold, make_sibs=make_sibs, progress=progress
    )
    bricked_graph = brick_grapher.make_brick_graph()
    return bricked_graph


def reduce_graph(
    brick_haplo_graph,
    brick_ts,
    path_weight_threshold,
    num_processes=1,
    progress=True,
    chunksize=100,
):
    """
    Make a reduced graph from a brick-haplo graph and bricked tree sequence

    :param networkx.DiGraph brick_haplo_graph: An input brick-haplo graph,
        outputted by `ldgm.brick_haplo_graph()`.
    :param tskit.TreeSequence brick_ts: The input :class:`tskit.TreeSequence`, which
        has been bricked (i.e. ``ldgm.brick_ts()`` has been run on it). This
        must be the tree sequence which was passed to `ldgm.brick_haplo_graph()`.
    :param float path_weight_threshold: The maximum path weight between nodes in
        the brick-haplo graph which should be retained in the LDGM. Lower
        `path_weight_threshold` values correspond to more regularization. The
        supplementary note of Nowbandegani et al. (2022) defines these path
        weights. Values of 4-6 have led to good results with data from the
        1000 Genomes Project.
    :param bool progress: Whether to display a progress bar. Default: False
    :return: An LDGM
    :rtype: networkx.DiGraph
    """
    assert utility.check_bricked(brick_ts)
    snp_grapher = reduction.SNP_Graph(
        brick_haplo_graph,
        brick_ts,
        path_weight_threshold=path_weight_threshold,
        num_processes=num_processes,
        progress=progress,
        chunksize=chunksize,
    )
    reduced_graph = snp_grapher.create_reduced_graph()
    return reduced_graph


def make_ldgm(
    ts,
    path_weight_threshold,
    recombination_freq_threshold=None,
    num_processes=1,
    chunksize=100,
    progress=False,
):
    """
    Take a tree sequence and produce an LDGM. The ``path_weight_threshold`` and
    ``recombination_freq_threshold`` parameters allow differing levels of
    regularization to be applied to this process.

    :param tskit.TreeSequence tree_sequence: The input :class:`tskit.TreeSequence`,
        which will be used to create the linkage disequilibrium graphical
        model.
    :param float path_weight_threshold: The maximum path weight between nodes in
        the brick-haplo graph which should be retained in the LDGM. Lower
        `path_weight_threshold` values correspond to more regularization. The
        supplementary note of Nowbandegani et al. (2022) defines these path
        weights. Values of 4-6 have led to good results with data from the 1000
        Genomes Project.
    :param float recombination_freq_threshold: The minimum frequency above which bricks
        should be created by "bifurcating" edges. For example, if
        `recombination_freq_threshold` is 0.01, an edge with a different number of
        descendant samples to the left and right of a recombination event (marginal
        tree change) is bifurcated to create two bricks *if* the frequency of
        the edge's child node is greater than 1% in the right hand marginal tree.
        If None, all edges with differing numbers of descendants to the left
        and right of a recombination event are bifurcated. Default: None
    :param int num_processes: The number of threads to use. Default: 1
    :param int chunksize: If using multiple threads, the algorithm will chop
        the dijkstra search step (the rate-limiting step of the algorithm)
        into chunks of nodes. The ``chunksize`` parameter determines the
        approximate size of chunks. Good results were observed with values
        around 100. Default: 100
    :param bool progress: Whether to display a progress bar. Default: False
    :return: A tuple of an LDGM and a bricked tree sequence.
    :rtype: (networkx.DiGraph, tskit.TreeSequence)
    """
    # Step 1: brick ts
    bts = brick_ts(
        ts, recombination_freq_threshold=recombination_freq_threshold, progress=progress
    )
    # Step 2: is brickhaplograph with no rule two or uturns
    bricked_graph = brick_haplo_graph(
        bts,
        path_weight_threshold,
        progress=progress,
        make_sibs=False,
    )
    # Step 3: compute reach* and create SNP-haplo graph
    H1 = reduce_graph(
        bricked_graph,
        bts,
        path_weight_threshold=path_weight_threshold,
        num_processes=num_processes,
        chunksize=chunksize,
        progress=progress,
    )
    # Step 4: brickhaplograph with rule two and uturns
    brickhaplograph_rule_two = brick_haplo_graph(
        bts, path_weight_threshold, make_sibs=True, progress=progress
    )
    # Step 5: reduce brickhaplograph created with rule two
    H2 = reduce_graph(
        brickhaplograph_rule_two,
        bts,
        path_weight_threshold=path_weight_threshold,
        num_processes=num_processes,
        progress=progress,
        chunksize=chunksize,
    )
    # Step 6: combine H1 and H2
    H_12 = nx.compose_all([H1, nx.reverse(H1), H2])
    # Step 7: Remove haplotype vertices to reduce H_12
    H_12_reduced = H_12.copy()
    nodes = list(H_12_reduced.nodes())
    for node in tqdm(nodes, total=len(nodes), desc="Removing nodes"):
        if node < 0:
            H_12_reduced = utility.remove_node(
                H_12_reduced,
                node,
                path_weight_threshold,
            )
    H_12_reduced = H_12_reduced.to_undirected()
    H_12_reduced_relabeled = utility.convert_node_ids(H_12_reduced, bts)

    return (H_12_reduced_relabeled, bts)


def prune_sites(ts, threshold):
    return utility.prune_sites(ts, threshold)


def make_snplist(ts, site_metadata_id=None, population_dict=None):
    """
    Returns information on variant sites in the input
        :class:`tskit.TreeSequence`.

    :param TreeSequence bricked_ts: The input :class`tskit.TreeSequence`,
        for which the site list will be determined. This tree sequence
        must have been "bricked", i.e. ldgm.brick_ts() must have been
        run on the tree sequence.
    :param string site_metadata_id: The site ID field in the JSON-encoded
        site metadata. If None, does not return site IDs. Default: None.
    :param list population_dict: A dictionary where keys are population
        ids and values are lists of nodes in each population. This defines
        populations from which allele frequencies are computed. If None,
        does not return site frequencies. Default: None.
    :return: A dictionary containing the following key, value pairs:
        key: "index", value: a numpy.ndarray of the IDs of the
        node in the LDGM to which each site belongs. When multiple
        mutations appear on the same brick, they will refer to the same node
        in the LDGM.
        key: "site_ids", value: a numpy.ndarray containing the site_id of
        each site, only returned if `site_metadata_id` is specified.
        key: "anc_alleles", value: a numpy.ndarray of the ancestral allele
        of each site.
        key: "deriv_alleles", value: a numpy.ndarray of the derived allele
        of each site.
        key: "site_frequencies", value: a mxp numpy.ndarray, where :math:`m`
        is the number of sites in the input tree sequence and :math:`s` is
        the number of sample sets passed to the `population_dict` parameter.
        Only returned if `population_dict` is specified.
    :rtype: dict
    """
    return utility.make_snplist(
        ts, site_metadata_id=site_metadata_id, population_dict=population_dict
    )


def return_edgelist(ldgm_graph):
    return utility.return_edgelist(ldgm_graph)
