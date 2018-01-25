# -*- coding: utf-8 -*-
r"""
Random partial injections and Stallings graphs

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

Random partial injections on ``[0, 1, ..., 6]``::

    sage: from slabbe import random_partial_injection
    sage: random_partial_injection(7)
    [None, None, 1, 3, None, 0, None]
    sage: random_partial_injection(7)
    [5, 1, 0, 3, None, 4, None]

Random Stallings graph on ``[0, 1, ..., 19]`` over 2 letters::

    sage: from slabbe import random_cyclically_reduced_stallings_graph
    sage: G,_,_ = random_cyclically_reduced_stallings_graph(20, 2)
    sage: G
    Looped multi-digraph on 20 vertices

Visualisation of the graph::

    sage: from slabbe import TikzPicture
    sage: tikz = TikzPicture.from_graph(G)
    sage: path_to_file = tikz.pdf()    # not tested

REFERENCES:

- Bassino, Frédérique; Nicaud, Cyril; Weil, Pascal Random generation of
  finitely generated subgroups of a free group.  Internat. J. Algebra
  Comput. 18 (2008), no. 2, 375–405. 

"""
#*****************************************************************************
#       Copyright (C) 2017 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
from sage.rings.integer_ring import ZZ
from sage.probability.probability_distribution import GeneralDiscreteDistribution
from sage.misc.prandom import shuffle, sample
from sage.graphs.digraph import DiGraph
from sage.functions.other import binomial, factorial

def number_of_partial_injection(n, algorithm='binomial'):
    r"""
    Return the number of partial injections on an set of `n` elements
    defined on a subset of `k` elements for each `k` in `0, 1, ..., n`.

    INPUT:

    - ``n`` -- integer
    - ``algorithm`` -- string (default: ``'binomial'``), ``'binomial'``
      or ``'recursive'``. When n>50, the binomial coefficient approach is
      faster (linear time vs quadratic time).

    OUTPUT:

        list

    .. NOTE::

        The recursive code of this function was originally written by
        Vincent Delecroix (Nov 30, 2017) the day after a discussion with
        Pascal Weil and me at LaBRI.

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

    TESTS::

        sage: number_of_partial_injection(8, algorithm='recursive')
        [1, 64, 1568, 18816, 117600, 376320, 564480, 322560, 40320]

    REFERENCE:

        https://oeis.org/A144084
    """
    if algorithm == 'binomial':
        return [binomial(n,k)**2*factorial(k) for k in range(n+1)]
    elif algorithm == 'recursive':
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

    TODO::

        Adapt the code once this is merged:

        https://trac.sagemath.org/ticket/24416

    AUTHORS:

    - Sébastien Labbé and Vincent Delecroix, Nov 30, 2017, Sage Thursdays
      at LaBRI
    """
    L = number_of_partial_injection(n)
    s = sum(L).n()
    L = [a/s for a in L]   # because of the bug #24416
    X = GeneralDiscreteDistribution(L)
    k = X.get_random_element()
    codomain = sample(range(n), k)
    missing = [None]*(n-k)
    codomain.extend(missing)
    shuffle(codomain)
    return codomain

def random_cyclically_reduced_stallings_graph(n, r=2, verbose=False, merge=False):
    r"""
    Return a uniformly chosen Stallings graph of n vertices over r letters.

    INPUT:

    - ``n`` -- integer, size of graph
    - ``r`` -- integer (default: ``2``), number of generators of the free
      group
    - ``verbose`` -- bool (default: ``False``)

    .. NOTE::

        The probability that G is connected is 1 - 2^r / n^(r-1) +
        o(1/n^(r-1)) which is approx. 1

    OUTPUT:

        digraph, integer, integer

    EXAMPLES::

        sage: from slabbe import random_cyclically_reduced_stallings_graph
        sage: G,_,_ = random_cyclically_reduced_stallings_graph(20, 2)
        sage: G
        Looped multi-digraph on 20 vertices

    ::

        sage: random_cyclically_reduced_stallings_graph(20, 5)[0]
        Looped multi-digraph on 20 vertices

    With verbose output::

        sage: G = random_cyclically_reduced_stallings_graph(20, 2, verbose=True)   # random
        rejecting because graph is not connected
        rejecting because graph has a vertex of degree <=1
        rejecting because graph has a vertex of degree <=1
        rejecting because graph has a vertex of degree <=1

    For displaying purposes, the following merges the multiedges
    automatically::

        sage: G,_,_ = random_cyclically_reduced_stallings_graph(20, 2)
        sage: from slabbe import TikzPicture
        sage: tikz = TikzPicture.from_graph(G)
        sage: _ = tikz.pdf(view=False)

    AUTHORS:

    - Sébastien Labbé and Pascal Weil, Dec 14, 2017, Sage Thursdays at LaBRI
    """
    from sage.misc.latex import LatexExpr

    # reject statistics
    not_connected_count = 0
    has_degree_1_count = 0

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
            not_connected_count += 1
            if verbose:
                print("rejecting because graph is not connected")
            continue

        if not all(d>=2 for d in G.degree_iterator()):
            has_degree_1_count += 1
            if verbose:
                print("rejecting because graph has a vertex of degree <=1")
            continue

        return G, not_connected_count, has_degree_1_count

def reject_statistics(n, r=2, sample_size=50, verbose=False):
    r"""
    Return return reject statistics when randomly chosing Stallings graph
    of n vertices over r letters.

    INPUT:

    - ``n`` -- integer, size of graph
    - ``r`` -- integer (default: ``2``), number of generators of the free
      group
    - ``n`` -- integer (default: ``50``), size of sample
    - ``verbose`` -- bool (default: ``False``)

    OUTPUT:

        histogram

    EXAMPLES::

        sage: from slabbe.partial_injection import reject_statistics
        sage: h = reject_statistics(50, verbose=True)   # random
        not connected: Counter({0: 48, 1: 2})
        has degree 1: Counter({0: 27, 1: 18, 2: 3, 3: 1, 4: 1})
        sage: h.save('h_50.png', title='size of graph=50') # not tested

    ::

        sage: h = reject_statistics(100, verbose=True)  # random
        not connected: Counter({0: 48, 1: 2})
        has degree 1: Counter({0: 41, 1: 8, 2: 1})
        sage: h.save('h_100.png', title='size of graph=100') # not tested

    ::

        sage: h = reject_statistics(500, verbose=True)       # not tested (30s)
        not connected: Counter({2: 5, 4: 5, 5: 5, 0: 4, 1: 4, 8: 4, 3: 3, 16:
            3, 6: 2, 11: 2, 15: 2, 18: 2, 23: 2, 7: 1, 10: 1, 44: 1, 13: 1, 49:
            1, 19: 1, 21: 1})
        has degree 1: Counter({0: 14, 1: 9, 3: 8, 2: 5, 4: 3, 5: 3, 6: 2, 9: 2,
            7: 1, 8: 1, 13: 1, 15: 1})
        sage: h.save('h_500.png', title='size of graph=500') # not tested

    ::

        sage: h = reject_statistics(1000, verbose=True)  # not tested (2min30s)
        not connected: Counter({8: 4, 26: 3, 3: 2, 7: 2, 9: 2, 10: 2, 14: 2,
        15: 2, 17: 2, 18: 2, 27: 2, 40: 2, 59: 2, 0: 1, 1: 1, 4: 1, 5: 1, 11:
        1, 13: 1, 19: 1, 20: 1, 21: 1, 22: 1, 28: 1, 44: 1, 48: 1, 51: 1, 52:
        1, 53: 1, 58: 1, 63: 1, 66: 1, 75: 1, 121: 1})
        has degree 1: Counter({2: 9, 0: 7, 1: 6, 4: 6, 3: 4, 5: 4, 7: 4, 8: 2,
            6: 1, 9: 1, 11: 1, 12: 1, 13: 1, 15: 1, 17: 1, 26: 1})
        sage: h.save('h_1000.png', title='size of graph=1000') # not tested

    """
    from sage.plot.histogram import histogram
    A = [random_cyclically_reduced_stallings_graph(n, r) for _ in range(sample_size)]
    _,s,t = zip(*A)
    if verbose:
        from collections import Counter
        print("not connected:", Counter(s))
        print("has degree 1:", Counter(t))
    h = histogram([s,t])
    return h

