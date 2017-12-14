# -*- coding: utf-8 -*-
r"""
Partial injections

EXAMPLES::

    sage: from slabbe import number_of_partial_injection
    sage: number_of_partial_injection(10)
    [1,
    100,
    4050,
    86400,
    1058400,
    7620480,
    31752000,
    72576000,
    81648000,
    36288000,
    3628800]

::

    sage: from slabbe import random_partial_injection
    sage: random_partial_injection(7)
    [None, None, 1, 3, None, 0, None]
    sage: random_partial_injection(7)
    [5, 1, 0, 3, None, 4, None]

REFERENCES:

- Bassino, Frédérique; Nicaud, Cyril; Weil, Pascal Random generation of
  finitely generated subgroups of a free group.  Internat. J. Algebra
  Comput. 18 (2008), no. 2, 375–405. 

AUTHORS:

- Sébastien Labbé and Vincent Delecroix, Nov 30, 2017, Sage Thursdays at LaBRI

"""
from sage.rings.integer_ring import ZZ
from sage.probability.probability_distribution import GeneralDiscreteDistribution
from sage.misc.prandom import shuffle, sample
from sage.graphs.digraph import DiGraph

def number_of_partial_injection(n):
    r"""
    Return the number of partial injections on an set of `n` elements
    defined on a subset of `k` elements for each `k` in `0, 1, ..., n`.

    INPUT:

    - ``n`` -- integer

    OUTPUT:

        list

    .. NOTE::

        The code of this function was written by Vincent Delecroix (Nov 30,
        2017) the day after a discussion with Pascal Weil and me at LaBRI.

    EXAMPLES::

        sage: from slabbe import number_of_partial_injection
        sage: number_of_partial_injection(0)
        [1]
        sage: number_of_partial_injection(1)
        [1, 1]
        sage: number_of_partial_injection(2)
        [1, 4, 2]
        sage: number_of_partial_injection(3)
        [1, 9, 18, 6]
        sage: number_of_partial_injection(4)
        [1, 16, 72, 96, 24]
        sage: number_of_partial_injection(5)
        [1, 25, 200, 600, 600, 120]
        sage: number_of_partial_injection(6)
        [1, 36, 450, 2400, 5400, 4320, 720]
        sage: number_of_partial_injection(7)
        [1, 49, 882, 7350, 29400, 52920, 35280, 5040]
        sage: number_of_partial_injection(8)
        [1, 64, 1568, 18816, 117600, 376320, 564480, 322560, 40320]
    """
    L = [ZZ(1)]
    for t in range(1, n+1):
        L.append(t * L[-1])
        for k in range(t-1, 0, -1):
            L[k] = (t * L[k]) // (t-k) + t * L[k-1]
    return L

def random_partial_injection(n):
    r"""
    Return a uniformly chosen random partial injection on 0, 1, ..., n-1.

    INPUT:

    - ``n`` -- integer

    OUTPUT:

        list

    EXAMPLES::

        sage: from slabbe import random_partial_injection
        sage: random_partial_injection(10)
        [3, 5, 2, None, 1, None, 0, 8, 7, 6]
        sage: random_partial_injection(10)
        [1, 7, 4, 8, 3, 5, 9, None, 6, None]
        sage: random_partial_injection(10)
        [5, 6, 8, None, 7, 4, 0, 9, None, None]
    """
    L = number_of_partial_injection(n)
    X = GeneralDiscreteDistribution(L)
    k = X.get_random_element()
    codomain = sample(range(n), k)
    missing = [None]*(n-k)
    codomain.extend(missing)
    shuffle(codomain)
    return codomain

def random_stallings_graph(n, r=2, verbose=False, merge=False):
    r"""
    Return a uniformly chosen Stallings graph of n vertices over r letters.

    INPUT:

    - ``n`` -- integer, size of graph
    - ``r`` -- integer (default: ``2``), number of generators of the free
      group
    - ``verbose`` -- bool (default: ``False``)
    - ``merge`` -- bool (default: ``False``), whether to merge the
      multiedges

    .. NOTE::

        The probability that G is connected is 1 - 2^r / n^(r-1) +
        o(1/n^(r-1)) which is approx. 1

    OUTPUT:

        digraph

    EXAMPLES::

        sage: from slabbe import random_stallings_graph
        sage: G = random_stallings_graph(20, 2)
        sage: G
        Looped multi-digraph on 20 vertices

    ::

        sage: random_stallings_graph(20, 5)
        Looped multi-digraph on 20 vertices

    With verbose output::

        sage: G = random_stallings_graph(20, 2, verbose=True)   # random
        rejecting because graph is not connected
        rejecting because graph has a vertex of degree <=1
        rejecting because graph has a vertex of degree <=1
        rejecting because graph has a vertex of degree <=1

    For displaying purpose, one can merge the multiedges::

        sage: G = random_stallings_graph(20, 2, merge=True)
        sage: from slabbe import TikzPicture
        sage: _ = TikzPicture.from_graph(G).pdf(view=False)

    AUTHORS:

    - Sébastien Labbé and Pascal Weil, Dec 14, 2017, Sage Thursdays at LaBRI

    """
    from sage.misc.latex import LatexExpr
    from slabbe.graph import merge_multiedges

    while True:
        injections = [random_partial_injection(n) for _ in range(r)]

        edges = []
        for i,injection in enumerate(injections):
            label = LatexExpr('a_{}'.format(i))
            edges.extend([(j,image_j,label) for (j,image_j) in enumerate(injection)
                                        if not image_j is None])

        G = DiGraph([range(n), edges], format='vertices_and_edges',
                loops=True, multiedges=True)

        if not G.is_connected():
            if verbose:
                print "rejecting because graph is not connected"
            continue

        if not all(d>=2 for d in G.degree_iterator()):
            if verbose:
                print "rejecting because graph has a vertex of degree <=1"
            continue

        if merge:
            G = merge_multiedges(G)
            edges = [(a,b,LatexExpr(label)) for (a,b,label) in G.edges()]
            return DiGraph(edges, format='list_of_edges', loops=True)
        else:
            return G

