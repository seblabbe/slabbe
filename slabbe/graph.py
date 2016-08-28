# -*- coding: utf-8 -*-
r"""
.. TODO::

    - Make the doctests more simple
"""
#*****************************************************************************
#       Copyright (C) 2016 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from collections import Counter
import itertools
from sage.graphs.digraph import DiGraph

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
        rows = [(value, key) for key,value in d.iteritems()]
        rows.sort(reverse=True,key=lambda row:row[1])
        header_row = ['Number of vertices', 'Projected vertices']
        from sage.misc.table import table
        print table(rows=rows, header_row=header_row)
    if filename:
        from slabbe import TikzPicture
        print TikzPicture.from_graph(G_proj, prog='dot').pdf(filename)
    return G_proj

def digraph_move_label_to_edge(G, label_function=None):
    r"""
    Return a digraph with labels moved from the arrival vertices to
    corresponding edges.

    INPUT:

    - ``G`` -- graph, whose vertices are tuples of the form (vertex, label)
    - ``label_function`` -- function or None, a function to apply to each label

    EXAMPLES::

        sage: G = DiGraph()
        sage: G.add_edges([((i, None), ((i+1)%10, 'plusone')) for i in range(10)])
        sage: G.add_edges([((i, None), ((i+2)%10, 'plustwo')) for i in range(10)])
        sage: G
        Digraph on 30 vertices
        sage: from slabbe.graph import digraph_move_label_to_edge
        sage: digraph_move_label_to_edge(G)
        Looped multi-digraph on 10 vertices

    Using a function to modify the labels::

        sage: f = lambda label:"A"+label
        sage: GG = digraph_move_label_to_edge(G, label_function=f)
        sage: GG
        Looped multi-digraph on 10 vertices
        sage: GG.edges()[0]
        (0, 1, 'Aplusone')
    """
    if label_function:
        edges = [(u,v,label_function(label)) for ((u,_), (v,label), _) in G.edges()]
    else:
        edges = [(u,v,label) for ((u,_), (v,label), _) in G.edges()]
    return DiGraph(edges, format='list_of_edges', loops=True, multiedges=True)

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
        sage: G.add_edges([((i, None), ((i+1)%10, 'plusone')) for i in range(10)])
        sage: G.add_edges([((i, None), ((i+2)%10, 'plustwo')) for i in range(10)])
        sage: G
        Digraph on 30 vertices
        sage: GG = induced_subgraph(G, lambda v: v[0]%2 == 0)
        sage: GG
        Digraph on 15 vertices
        sage: GG.edges()[0]
        ((0, None), (2, 'plustwo'), None)

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

