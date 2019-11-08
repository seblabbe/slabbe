# -*- coding: utf-8 -*-
r"""
Functions on graphs

.. TODO::

    - Make the doctests more simple
"""
#*****************************************************************************
#       Copyright (C) 2016-2017 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
from collections import Counter, defaultdict
import itertools
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph

def projection_graph(G, proj_fn, filename=None, verbose=False):
    r"""
    Return the image of a graph under a function on vertices.

    INPUT:

    - ``G`` -- graph
    - ``proj_fn`` -- function
    - ``filename`` -- integer (default:``None``), save the graph to this pdf
      filename if filename is not None
    - ``verbose`` -- bool (default:``False``), print a table of data about the
      projection

    EXAMPLES::

        sage: from slabbe.graph import projection_graph
        sage: g = graphs.PetersenGraph()
        sage: g.vertices()
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: f = lambda i: i % 5
        sage: projection_graph(g, f)
        Looped multi-digraph on 5 vertices

    With verbose information::

        sage: projection_graph(g, lambda i:i%4, verbose=True)
          Number of vertices   Projected vertices
        +--------------------+--------------------+
          2                    3
          2                    2
          3                    1
          3                    0
        Looped multi-digraph on 4 vertices
    """
    edges = set((proj_fn(A),proj_fn(B)) for A,B,_ in G.edges())
    G_proj = DiGraph(edges, format='list_of_edges', loops=True, multiedges=True)
    if verbose:
        d = dict(Counter(proj_fn(s) for s in G.vertices()))
        rows = [(value, key) for key,value in d.items()]
        rows.sort(reverse=True,key=lambda row:row[1])
        header_row = ['Number of vertices', 'Projected vertices']
        from sage.misc.table import table
        print(table(rows=rows, header_row=header_row))
    if filename:
        from slabbe import TikzPicture
        print(TikzPicture.from_graph(G_proj, prog='dot').pdf(filename))
    return G_proj

def digraph_move_label_to_edge(G, label_function=None, loops=True,
        multiedges=False):
    r"""
    Return a digraph with labels moved from the arrival vertices to
    corresponding edges.

    INPUT:

    - ``G`` -- graph, whose vertices are tuples of the form (vertex, label)
    - ``label_function`` -- function or None, a function to apply to each label
    - ``loops`` -- bool (default: True)
    - ``multiedges`` -- bool (default: False)

    EXAMPLES::

        sage: G = DiGraph()
        sage: G.add_edges([((i, None), ((i+1)%10, 'plusone')) for i in range(10)])
        sage: G.add_edges([((i, None), ((i+2)%10, 'plustwo')) for i in range(10)])
        sage: G
        Digraph on 30 vertices
        sage: from slabbe.graph import digraph_move_label_to_edge
        sage: digraph_move_label_to_edge(G)
        Looped digraph on 10 vertices

    Using a function to modify the labels::

        sage: f = lambda label:"A"+label
        sage: GG = digraph_move_label_to_edge(G, label_function=f)
        sage: GG
        Looped digraph on 10 vertices
        sage: GG.edges()[0]
        (0, 1, 'Aplusone')
    """
    if label_function:
        edges = [(u,v,label_function(label)) for ((u,_), (v,label), _) in G.edges()]
    else:
        edges = [(u,v,label) for ((u,_), (v,label), _) in G.edges()]
    return DiGraph(edges, format='list_of_edges', loops=loops,
            multiedges=multiedges)

def induced_subgraph(G, filter):
    r"""
    Return the induced subdigraph of a digraph keeping only vertices that are
    map to ``True`` by the filter.

    INPUT:

    - ``G`` -- graph
    - ``filter`` -- function, a function from vertices to boolean

    EXAMPLES::

        sage: from slabbe.graph import induced_subgraph
        sage: G = DiGraph()
        sage: G.add_edges([((i, ''), ((i+1)%10, 'plusone')) for i in range(10)])
        sage: G.add_edges([((i, ''), ((i+2)%10, 'plustwo')) for i in range(10)])
        sage: G
        Digraph on 30 vertices
        sage: GG = induced_subgraph(G, lambda v: v[0]%2 == 0)
        sage: GG
        Digraph on 15 vertices
        sage: GG.edges()[0]
        ((0, ''), (2, 'plustwo'), None)

    .. TODO::

        simplify the edges aaaaa*aaa*aaa* to a* only
    """
    GG = G.copy()
    loops = dict((u, label) for (u,v,label) in GG.loop_edges())
    for v in GG.vertices():
        if filter(v):
            continue
        incoming = [(x,y,l) for (x,y,l) in GG.incoming_edges(v) if x != y]
        outgoing = [(x,y,l) for (x,y,l) in GG.outgoing_edges(v) if x != y]
        GG.delete_vertex(v)
        if v in loops:
            it = itertools.product(outgoing, [(loops[v][0]*10,)], incoming) 
        else:
            it = itertools.product(outgoing, [tuple()], incoming) 
        for c,b,a in it:
            _,y,labelout = c
            x,_,labelin = a
            GG.add_edge(x,y, labelout + b + labelin)
    return GG

def merge_multiedges(G, label_function=tuple):
    r"""
    Return the (di)graph where multiedges are merged into one.

    INPUT:

    - ``G`` -- graph
    - ``label_function`` -- function (default:``tuple``), a function to
      apply to each list of labels

    OUTPUT:

        (looped) (di)graph

    EXAMPLES:

    A digraph::

        sage: from slabbe.graph import merge_multiedges
        sage: G = DiGraph(multiedges=True)
        sage: G.add_edge(0,1,'one')
        sage: G.add_edge(0,1,'two')
        sage: G.add_edge(0,1,'alpha')
        sage: GG = merge_multiedges(G)
        sage: GG
        Digraph on 2 vertices
        sage: GG.edges()
        [(0, 1, ('alpha', 'one', 'two'))]

    A graph::

        sage: G = Graph(multiedges=True)
        sage: G.add_edge(0,1,'one')
        sage: G.add_edge(0,1,'two')
        sage: G.add_edge(0,1,'alpha')
        sage: GG = merge_multiedges(G)
        sage: GG
        Graph on 2 vertices
        sage: GG.edges()
        [(0, 1, ('alpha', 'one', 'two'))]

    Using ``label_function``::

        sage: fn = lambda L: LatexExpr(','.join(map(str, L)))
        sage: GG = merge_multiedges(G, label_function=fn)
        sage: GG.edges()
        [(0, 1, alpha,one,two)]

    """
    d = defaultdict(list)
    for (u,v,label) in G.edges():
        d[(u,v)].append(label)

    edges = [(u,v,label_function(label_list)) for (u,v),label_list in d.items()]

    loops = G.has_loops()
    if G.is_directed():
        return DiGraph(edges, format='list_of_edges', loops=loops)
    else:
        return Graph(edges, format='list_of_edges', loops=loops)

def clean_sources_and_sinks(G):
    r"""
    Return a copy of the graph where every vertices of the graph that have
    in or out degree 0 is removed (recursively).

    EXAMPLES::

        sage: from slabbe.graph import clean_sources_and_sinks
        sage: L = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,3)]
        sage: G = DiGraph(L,format='list_of_edges')
        sage: H = clean_sources_and_sinks(G)
        sage: H
        Digraph on 3 vertices
        sage: H.vertices()
        [3, 4, 5]

    ::

        sage: L = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,3),(1,0)]
        sage: G = DiGraph(L, format='list_of_edges')
        sage: H = clean_sources_and_sinks(G)
        sage: H
        Digraph on 6 vertices
        sage: H.vertices()
        [0, 1, 2, 3, 4, 5]

    """
    H = G.copy()
    done = False
    while not done:
        done = True
        for v in H.vertices():
            if H.in_degree(v) == 0 or H.out_degree(v) == 0:
                done = False
                H.delete_vertex(v)
    return H

def get_funnel(G):
    r"""
    Return an edge (u,v) such that u and v are distinct, G.out_degree(u) is
    1 and G.in_degree(v) is 1. Return ``None`` if no such funnel is found.

    INPUT:

    - ``G`` -- digraph

    EXAMPLES::

        sage: from slabbe.graph import get_funnel
        sage: G = DiGraph([(str(a),str(a+1)) for a in range(5)], format='list_of_edges')
        sage: get_funnel(G)
        ('0', '1')
    """
    for (u,v,_) in G.edges(): 
        if (u != v and G.in_degree(v) == 1 and G.out_degree(u) == 1):
            return (u,v)
    else:
        return None

def reduce_funnel_edges(G, merge_function):
    r"""
    EXAMPLES::

        sage: from slabbe.graph import reduce_funnel_edges
        sage: G = DiGraph([(str(a),str(a+1)) for a in range(5)], format='list_of_edges')
        sage: merge_function = lambda a,b:a+b
        sage: GG = reduce_funnel_edges(G, merge_function)
        sage: GG.vertices()
        ['012345']

    ::

        sage: G = DiGraph([(str(a),str((a+1)%5)) for a in range(5)], format='list_of_edges')
        sage: merge_function = lambda a,b:a+b
        sage: GG = reduce_funnel_edges(G, merge_function)
        sage: GG.vertices()
        ['01234']

    The following result does not seem right::

        sage: w = words.FibonacciWord()[:100]
        sage: G = w.rauzy_graph(11)
        sage: merge_function = lambda a,b:a+b[-1:]
        sage: GG = reduce_funnel_edges(G, merge_function)
        sage: GG.vertices()
        [word: 01001010010, word: 100101001001, word: 100101001011]

    """
    from copy import copy
    GG = copy(G)
    GG.allow_loops(True)
    #GG.allow_multiple_edges(True)

    funnel = get_funnel(GG)
    while funnel:
        (u,v) = funnel
        u_v = merge_function(u, v)
        GG.add_vertex(u_v)
        for s in GG.neighbors_out(v):
            GG.add_edge(u_v, s)
        for s in GG.neighbors_in(u):
            GG.add_edge(s, u_v)
        GG.delete_vertex(u)
        GG.delete_vertex(v)
        funnel = get_funnel(GG)
    return GG

