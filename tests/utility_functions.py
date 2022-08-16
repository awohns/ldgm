"""
A collection of utilities to edit and construct tree sequences for testing purposes
"""
import io

import tskit


def single_tree_ts_n2():
    r"""
    Simple case where we have n = 2 and one tree. [] marks a sample
         2
        / \
      [0] [1]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       0           1
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       2       0,1
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def single_tree_ts_n2_2_mutations():
    r"""
    Simple case where we have n = 2, one tree, and two mutations.
         2
        x x
      [0] [1]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       0           1
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       2       0,1
    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.5         0
    0.6         0
    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       0       1
    0       1       1
    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def single_tree_ts_n3():
    r"""
    Simple case where we have n = 3 and one tree.
            4
           / \
          3   \
         / \   \
       [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    4       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       4       2,3
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def single_tree_ts_n4():
    r"""
    Simple case where we have n = 4 and one tree.
              6
             / \
            5   \
           / \   \
          4   \   \
         / \   \   \
       [0] [1] [2] [3]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       1           0
    4       0           1
    5       0           2
    6       0           3
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       4       0,1
    0       1       5       2,4
    0       1       6       3,5
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def single_tree_ts_mutation_n3():
    r"""
    Simple case where we have n = 3 and one tree.
            4
           / \
          3   x
         / \   \
       [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    4       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       4       2,3
    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.5         0
    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       2       1
    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def site_no_mutations():
    r"""
    Simple case where we have n = 3 and one tree.
    The single site has no derived alleles.
            4
           / \
          3   x
         / \   \
       [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    4       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       4       2,3
    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.5         0
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, sites=sites, strict=False)


def single_tree_all_samples_one_mutation_n3():
    r"""
    Simple case where we have n = 3 and one tree.
            4
           / \
          3   x
         / \   \
       [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       1           1
    4       1           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       4       2,3
    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.5         0
    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       2       1
    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def gils_example_tree():
    r"""
    Simple case where we have n = 3 and one tree.
    Mutations marked on each branch by *.
             4
            / \
           /   \
          /     *
         3       *
        / \       *
       *   *       *
      *     \       \
    [0]     [1]     [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    4       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       4       2,3
    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.1         0
    0.2         0
    0.3         0
    0.4         0
    0.5         0
    0.6         0
    0.7         0
    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       0       1
    1       0       1
    2       1       1
    3       2       1
    4       2       1
    5       2       1
    6       2       1
    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def polytomy_tree_ts():
    r"""
    Simple case where we have n = 3 and a polytomy.
          3
         /|\
        / | \
      [0][1][2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1,2
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def single_tree_ts_internal_n3():
    r"""
    Simple case where we have n = 3 and one tree.
    Node 3 is an internal sample.
            4
           / \
          3   \
         / \   \
       [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       1           1
    4       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       4       2,3
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def two_tree_ts():
    r"""
    Simple case where we have n = 3 and 2 trees.
                   .    5
                   .   / \
            4      .  |   4
           / \     .  |   |\
          3   \    .  |   | \
         / \   \   .  |   |  \
       [0] [1] [2] . [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    4       0           2
    5       0           3
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       0.2     3       0,1
    0       1       4       2
    0       0.2     4       3
    0.2     1       4       1
    0.2     1       5       0,4
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def two_tree_ts_extra_length():
    r"""
    Simple case where we have n = 3 and 2 trees, but with extra length
    for testing keep_intervals() and delete_intervals().
                   .    5
                   .   / \
            4      .  |   4
           / \     .  |   |\
          3   \    .  |   | \
         / \   \   .  |   |  \
       [0] [1] [2] . [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    4       0           2
    5       0           3
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       0.2     3       0,1
    0       1.5     4       2
    0       0.2     4       3
    0.2     1.5     4       1
    0.2     1.5     5       0,4
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def two_tree_ts_n3_non_contemporaneous():
    r"""
    Simple case where we have n = 3 and two trees with node 2 ancient.
                   .    5
                   .   / \
            4      .  |   4
           / \     .  |   |\
          3  [2]   .  |   |[2]
         / \       .  |   |
       [0] [1]     . [0] [1]
    """
    ts = two_tree_ts()
    tables = ts.dump_tables()
    time = tables.nodes.time
    time[2] = time[3]
    tables.nodes.time = time
    return tables.tree_sequence()


def single_tree_ts_with_unary():
    r"""
    Simple case where we have n = 3 and some unary nodes.
            7
           / \
          5   \
          |    \
          4     6
          |     |
          3     |
         / \    |
       [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    4       0           2
    5       0           3
    6       0           2
    7       0           4
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       6       2
    0       1       4       3
    0       1       5       4
    0       1       7       5,6
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def two_tree_ts_with_unary_n3():
    r"""
    Simple case where we have n = 3 and node 5 is an internal, unary node in the first
    tree. In the second tree, node t is the root, but still unary.
             6        .      5
           /   \      .      |
          4     5     .      4
          |     |     .     /  \
          3     |     .    3    \
         / \    |     .   / \    \
       [0] [1] [2]    . [0]  [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    4       0           2
    5       0           3
    6       0           4
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       2       3       0,1
    0       1       5       2
    0       2       4       3
    0       1       6       4,5
    1       2       4       2
    1       2       5       4
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def two_tree_mutation_ts():
    r"""
    Simple case where we have n = 3, 2 trees, three mutations.
                   .     5
                   .    / \
            4      .   |   4
           / \     .   |   |\
          x   \    .   |   | \
         x     \   .   x   |  \
        /      |   .   |   |   |
       3       |   .   |   |   |
      / \      |   .   |   |   |
    [0] [1]   [2]  .  [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    4       0           2
    5       0           3
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       0.2     3       0,1
    0       1       4       2
    0       0.2     4       3
    0.2     1       4       1
    0.2     1       5       0,4
    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.1         0
    0.15         0
    0.8         0
    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       3       1
    1       3       1
    2       0       1
    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def two_tree_two_mrcas():
    r"""
    Simple case where we have n = 4, 2 trees, one mutation.
             6             |
            / \            |            7
           /   \           |           / \
          /     \          |          /   x
         /       \         |         /     \
        /         \        |        /       \
       4           5       |       4         5
      / \         / \      |      / \       / \
     /   \       /   \     |     /   \     /   \
   [0]   [1]   [2]   [3]   |   [0]   [1] [2]   [3]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       1           0
    4       0           1
    5       0           1
    6       0           3
    7       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       4       0,1
    0       1       5       2,3
    0       0.3     6       4
    0       0.3     6       5
    0.3     1       7       4
    0.3     1       7       5
    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.5         0
    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       5       1
    """
    )

    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def loopy_tree():
    r"""
    Testing a tree with a loop.
                   .          7
                   .         / \
                   .        /   |
                   .       /    |
         6         .      /     6
        / \        .     /     / \
       /   5       .    /     /   |
      /   / \      .   /     /    |
     /   |   \     .  |     |     |
    /    |    \    .  |     |     |
   |     4     |   .  |     4     |
   |    / \    |   .  |    / \    |
  [0] [1] [2] [3]  . [0] [1] [2] [3]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       1           0
    4       0           1
    5       0           2
    6       0           3
    7       0           4
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       4       0,1
    0       0.2     5       2,4
    0       0.2     6       5
    0       1       6       3
    0.2     1       6       4
    0.2     1       7       2
    0.2     1       7       6
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def single_tree_ts_n3_sample_as_parent():
    r"""
    Simple case where we have n = 3 and one tree. Node 3 is a sample.
            4
           / \
          3   \
         / \   \
       [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       1           1
    4       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       4       2,3
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def single_tree_ts_n2_dangling():
    r"""
    Simple case where we have n = 2 and one tree. Node 0 is dangling.
            4
           / \
          3   \
         / \   \
       [0] [1]  2
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       0           0
    3       0           1
    4       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       4       2,3
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def two_tree_ts_n2_part_dangling():
    r"""
    Simple case where we have n = 2 and two trees. Node 0 is dangling in the first tree.
            4                 4
           / \               / \
          3   \             3   \
         / \   \             \   \
      [0]   \   \            [0]  \
             \   \             \   \
             [1]  2           [1]   2
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0.5
    1       1           0
    2       0           0
    3       0           1
    4       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0
    0       0.5     3       1
    0.5     1       0       1
    0       1       4       2,3
    """
    )
    return tskit.load_text(nodes=nodes, edges=edges, strict=False)


def single_tree_ts_2mutations_multiallelic_n3():
    r"""
    Simple case where we have n = 3 and one tree.
    Site is multiallelic.
            4
           x \
          3   x
         / \   \
       [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    4       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       4       2,3
    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.5         0
    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       2       1
    0       3       2
    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def single_tree_ts_2mutations_singletons_n3():
    r"""
    Simple case where we have n = 3 and one tree.
    Site has two singleton mutations.
            4
           / \
          3   x
         / x   \
       [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           1
    4       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       4       2,3
    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.5         0
    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       1       1
    0       2       1
    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def single_tree_ts_2mutations_n3():
    r"""
    Simple case where we have n = 3 and one tree.
    Site has two mutations with different times.
            4
           x \
          3   \
         / x   \
       [0] [1] [2]
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       0           10
    4       0           20
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       3       0,1
    0       1       4       2,3
    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.5         0
    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       3       1
    0       1       0
    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def figure_one_example():
    r"""
    Simple case where we have n = 5 and 2 trees.
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       1           0
    4       1           0
    5       0           1
    6       0           1
    7       0           2
    8       0           2
    9       0           3

    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       7       0
    0       1       5       1
    0       0.5     5       2
    0.5     1       6       2
    0       1       7       5
    0       1       9       7,8
    0       1       8       6
    0       1       6       3
    0       1       8       4

    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.1         0
    0.15         0
    0.8         0
    0.9         0

    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       7       1
    1       5       1
    2       8       1
    3       6       1

    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def supplementary_example():
    r"""
    Simple case where we have n = 5 and 2 trees.
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       1           0
    4       1           0
    5       0           1
    6       0           1
    7       0           2
    8       0           2
    9       0           3

    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       7       0
    0       1       5       1
    0       0.5     5       2
    0.5     1       6       2
    0       1       7       5
    0       1       9       7,8
    0       1       8       6
    0       1       6       3
    0       1       8       4

    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.1         0
    0.15         0
    0.8         0
    0.9         0

    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       7       1
    1       6       1
    2       5       1
    3       8       1

    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def triangle_example():
    r"""
    Reduced graph should be a triangle
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       0           1
    3       1           0
    4       0           2
    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       2       0
    0       1       2       1
    0       1       4       2
    0       1       4       3
    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.8         0
    0.91        0
    0.95        0
    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       0       1
    1       1       1
    2       2       1
    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )


def multiple_snps_branch():
    r"""
    n = 5 and 2 trees, with two mutations on a branch
    """
    nodes = io.StringIO(
        """\
    id      is_sample   time
    0       1           0
    1       1           0
    2       1           0
    3       1           0
    4       1           0
    5       0           1
    6       0           1
    7       0           2
    8       0           2
    9       0           3

    """
    )
    edges = io.StringIO(
        """\
    left    right   parent  child
    0       1       7       0
    0       1       5       1
    0       0.5     5       2
    0.5     1       6       2
    0       1       7       5
    0       1       9       7,8
    0       1       8       6
    0       1       6       3
    0       1       8       4

    """
    )
    sites = io.StringIO(
        """\
    position    ancestral_state
    0.1         0
    0.15        0
    0.8         0
    0.9         0
    0.95        0
    """
    )
    mutations = io.StringIO(
        """\
    site    node    derived_state
    0       7       1
    1       6       1
    2       5       1
    3       5       1
    4       8       1

    """
    )
    return tskit.load_text(
        nodes=nodes, edges=edges, sites=sites, mutations=mutations, strict=False
    )
