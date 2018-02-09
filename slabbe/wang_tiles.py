# -*- coding: utf-8 -*-
r"""
Wang tile solver

This uses MILP solvers like Coin or Gurobi. Coin can be installed with::

    sage -i cbc sagelib

EXAMPLES::

    sage: from slabbe import WangTileSolver
    sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: tiling = W.solve()
    sage: _ = tiling.tikz().pdf(view=False)

::

    sage: tiles = [(1/2,1/2,1/2,1/2), (1,1,1,1), (2,2,2,2)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: tiling = W.solve()
    sage: _ = tiling.tikz().pdf(view=False)

::

    sage: tiles = [(0,3,1,4), (1,4,0,3)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: tiling = W.solve()
    sage: tiling._table
    [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]

Kari-Culik (Here 0' is replaced by 10)::

    sage: divide_by_2 = [(10,0,10,10), (10,1,10,2), (1/2,0,10,1), (1/2,10,10,1),
    ....:     (1/2,0,1/2,10), (1/2,1,1/2,2), (10,1,1/2,1)]
    sage: times_3 = [(1,2,0,1), (2,1,0,1), (2,2,1,1), (0,1,1,0), (0,2,2,0),
    ....:     (1,1,2,0)]
    sage: tiles = divide_by_2 + times_3
    sage: W = WangTileSolver(tiles,3,4)
    sage: tiling = W.solve()
    sage: _ = tiling.tikz().pdf(view=False)

Rao-Jeandel::

    sage: tiles = [(2,4,2,1), (2,2,2,0), (1,1,3,1), (1,2,3,2), (3,1,3,3),
    ....: (0,1,3,1), (0,0,0,1), (3,1,0,2), (0,2,1,2), (1,2,1,4), (3,3,1,2)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: tiling = W.solve()
    sage: _ = tiling.tikz().pdf(view=False)

.. TODO::

 - Reduce the problem to SAT
                                         
"""
#*****************************************************************************
#       Copyright (C) 2017-2018 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
import itertools
from collections import Counter, defaultdict
from sage.misc.cachefunc import cached_method
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.rings.integer_ring import ZZ

from sage.misc.decorators import rename_keyword

def tile_to_tikz(tile, position, color=None, size=1,
        rotate=None, label_shift=.2, top_right_edges=True):
    r"""

    INPUT:

    - ``tile`` -- tuple of length 4
    - ``position`` -- tuple of two numbers
    - ``color`` -- dict (default: ``None``) from tile values -> tikz colors
    - ``size`` -- number (default: ``1``), size of the tile
    - ``rotate`` -- list or ``None`` (default:``None``) list of four angles
      in degrees like ``(0,0,0,0)``, the rotation angle to apply to each
      label of Wang tiles. If ``None``, it performs a 90 degres rotation
      for left and right labels taking more than one character.
    - ``label_shift`` -- number (default: ``.2``) translation distance of the
      label from the edge
    - ``top_right_edges`` -- bool (default: ``True``) 

    OUTPUT:

    - list of strings

    EXAMPLES::

        sage: from slabbe.wang_tiles import tile_to_tikz
        sage: color = {0:'white',1:'red',2:'cyan',3:'green',4:'white'}
        sage: tile_to_tikz((1,2,3,4), (10,100), color)
        ['% tile at position (x,y)=(10, 100)',
         '\\fill[red] (11, 100) -- (10.5, 100.5) -- (11, 101);',
         '\\fill[cyan] (10, 101) -- (10.5, 100.5) -- (11, 101);',
         '\\fill[green] (10, 100) -- (10.5, 100.5) -- (10, 101);',
         '\\fill[white] (10, 100) -- (10.5, 100.5) -- (11, 100);',
         '\\draw (10, 100) -- (11, 100);',
         '\\draw (10, 100) -- (10, 101);',
         '\\draw (11, 101) -- (11, 100);',
         '\\draw (11, 101) -- (10, 101);',
         '\\node[rotate=0] at (10.8, 100.5) {1};',
         '\\node[rotate=0] at (10.5, 100.8) {2};',
         '\\node[rotate=0] at (10.2, 100.5) {3};',
         '\\node[rotate=0] at (10.5, 100.2) {4};']
        sage: tile_to_tikz((1,2,3,4), (10,100), color=None)
        ['% tile at position (x,y)=(10, 100)',
         '\\draw (10, 100) -- (11, 100);',
         '\\draw (10, 100) -- (10, 101);',
         '\\draw (11, 101) -- (11, 100);',
         '\\draw (11, 101) -- (10, 101);',
         '\\node[rotate=0] at (10.8, 100.5) {1};',
         '\\node[rotate=0] at (10.5, 100.8) {2};',
         '\\node[rotate=0] at (10.2, 100.5) {3};',
         '\\node[rotate=0] at (10.5, 100.2) {4};']
        sage: tile_to_tikz((1,2,3,4), (10,100), color=None, rotate=(0,90,0,0))
        ['% tile at position (x,y)=(10, 100)',
         '\\draw (10, 100) -- (11, 100);',
         '\\draw (10, 100) -- (10, 101);',
         '\\draw (11, 101) -- (11, 100);',
         '\\draw (11, 101) -- (10, 101);',
         '\\node[rotate=0] at (10.8, 100.5) {1};',
         '\\node[rotate=90] at (10.5, 100.8) {2};',
         '\\node[rotate=0] at (10.2, 100.5) {3};',
         '\\node[rotate=0] at (10.5, 100.2) {4};']
        sage: tile_to_tikz((1,2,3,4), (10,100), color=None, label_shift=.1)
        ['% tile at position (x,y)=(10, 100)',
         '\\draw (10, 100) -- (11, 100);',
         '\\draw (10, 100) -- (10, 101);',
         '\\draw (11, 101) -- (11, 100);',
         '\\draw (11, 101) -- (10, 101);',
         '\\node[rotate=0] at (10.9000000000000, 100.5) {1};',
         '\\node[rotate=0] at (10.5, 100.900000000000) {2};',
         '\\node[rotate=0] at (10.1000000000000, 100.5) {3};',
         '\\node[rotate=0] at (10.5, 100.100000000000) {4};']

    ::

        sage: tile_to_tikz((10,20,30,40), (10,100), color=None)
        ['% tile at position (x,y)=(10, 100)',
         '\\draw (10, 100) -- (11, 100);',
         '\\draw (10, 100) -- (10, 101);',
         '\\draw (11, 101) -- (11, 100);',
         '\\draw (11, 101) -- (10, 101);',
         '\\node[rotate=90] at (10.8, 100.5) {10};',
         '\\node[rotate=0] at (10.5, 100.8) {20};',
         '\\node[rotate=90] at (10.2, 100.5) {30};',
         '\\node[rotate=0] at (10.5, 100.2) {40};']
    """
    if rotate is None:
        rotate = []
        for i,a in enumerate(tile):
            if (i == 0 or i == 2) and len(str(a)) > 1:
                rotate.append(90)
            else:
                rotate.append(0)
    s = size        # because it is shorter to write below
    t = label_shift # because it is shorter to write below
    x,y = position
    lines = []
    #lines.append(r'\begin{tikzpicture}')
    lines.append('% tile at position (x,y)={}'.format((x,y)))
    if color:
        triangle = r'\fill[{}] {} -- {} -- {};'
        c = (x+.5*s,y+.5*s)
        lines.append(triangle.format(color[tile[0]],(x+s,y),c,(x+s,y+s)))
        lines.append(triangle.format(color[tile[1]],(x,y+s),c,(x+s,y+s)))
        lines.append(triangle.format(color[tile[2]],(x,y),c,(x,y+s)))
        lines.append(triangle.format(color[tile[3]],(x,y),c,(x+s,y)))
    lines.append(r'\draw {} -- {};'.format((x,y), (x+s,y)))
    lines.append(r'\draw {} -- {};'.format((x,y), (x,y+s)))
    if top_right_edges:
        lines.append(r'\draw {} -- {};'.format((x+s,y+s), (x+s,y)))
        lines.append(r'\draw {} -- {};'.format((x+s,y+s), (x,y+s)))
    node_str = r'\node[rotate={}] at {} {{{}}};'
    lines.append(node_str.format(rotate[0],(x+s-t,y+.5),  tile[0]))
    lines.append(node_str.format(rotate[1],(x+.5, y+s-t), tile[1]))
    lines.append(node_str.format(rotate[2],(x+t,  y+.5),  tile[2]))
    lines.append(node_str.format(rotate[3],(x+.5, y+t),   tile[3]))
    #lines.append(r'\end{tikzpicture}')
    return lines
    #return TikzPicture('\n'.join(lines))


def hexagonal_tile_to_tikz(tile, position, color=None, radius=1,
        rotate=(0,0,0,0,0,0), label_shift=.2, digits=5):
    r"""

    INPUT:

    - ``tile`` -- tuple of length 6
    - ``position`` -- tuple of two numbers
    - ``color`` -- dict (default: ``None``) from tile values -> tikz colors
    - ``radius`` -- number (default: ``1``), radius of the tile
    - ``rotate`` -- list (default:``(0,0,0,0,0,0)``) of four angles in
      degrees, the rotation angle to apply to each label of Wang tiles
    - ``label_shift`` -- number (default: ``.2``) translation distance of the
      label from the edge
    - ``digits`` -- number (default: ``5``), number of digits of coordinates

    OUTPUT:

    - list of strings

    EXAMPLES::

        sage: from slabbe.wang_tiles import hexagonal_tile_to_tikz
        sage: color = {0:'white',1:'red',2:'cyan',3:'green',4:'white'}
        sage: hexagonal_tile_to_tikz((1,2,3,4,1,2), (10,10), color)
        ['% hexagonal tile at position (x,y)=(10, 10)',
         '\\fill[red] (10, 10) -- (10.866, 9.5000) -- (10.866, 10.500) -- cycle;',
         '\\fill[cyan] (10, 10) -- (10.866, 10.500) -- (10.000, 11.000) -- cycle;',
         '\\fill[green] (10, 10) -- (10.000, 11.000) -- (9.1340, 10.500) -- cycle;',
         '\\fill[white] (10, 10) -- (9.1340, 10.500) -- (9.1340, 9.5000) -- cycle;',
         '\\fill[red] (10, 10) -- (9.1340, 9.5000) -- (10.000, 9.0000) -- cycle;',
         '\\fill[cyan] (10, 10) -- (10.000, 9.0000) -- (10.866, 9.5000) -- cycle;',
         '\\draw (10.866, 9.5000) -- (10.866, 10.500) -- (10.000, 11.000) -- (9.1340, 10.500) -- (9.1340, 9.5000) -- (10.000, 9.0000) -- (10.866, 9.5000);',
         '\\node[rotate=0] at (10.666, 10.000) {1};',
         '\\node[rotate=0] at (10.333, 10.577) {2};',
         '\\node[rotate=0] at (9.6670, 10.577) {3};',
         '\\node[rotate=0] at (9.3340, 10.000) {4};',
         '\\node[rotate=0] at (9.6670, 9.4232) {1};',
         '\\node[rotate=0] at (10.333, 9.4232) {2};']
    """
    from sage.modules.free_module_element import vector
    from sage.symbolic.constants import pi
    from sage.functions.trig import sin,cos
    from slabbe import TikzPicture
    lines = []
    #lines.append(r'\begin{tikzpicture}')
    x,y = position = vector(position)
    lines.append('% hexagonal tile at position (x,y)={}'.format(position))
    angles = [-pi/6+i*pi/3 for i in range(7)]
    vertices = [position+vector((cos(angle), sin(angle))) for angle in angles]
    if color:
        triangle = r'\fill[{}] {} -- {} -- {} -- cycle;'
        c = (x,y)
        for i in range(6):
            p = vertices[i]
            q = vertices[i+1]
            lines.append(triangle.format(color[tile[i]],c,p.n(digits=digits),q.n(digits=digits)))
    contour = ' -- '.join(['{}'.format(v.n(digits=digits)) for v in vertices])
    lines.append(r'\draw {};'.format(contour))
    node_str = r'\node[rotate={}] at {} {{{}}};'
    for i in range(6):
        angle = i*pi/3
        coeff = radius*cos(pi/6) - label_shift
        p = position + coeff*vector((cos(angle), sin(angle)))
        lines.append(node_str.format(rotate[i],p.n(digits=digits),tile[i]))
    return lines
    #lines.append(r'\end{tikzpicture}')
    #return TikzPicture('\n'.join(lines))


class WangTileSet_generic(object):
    r"""
    Construct a Wang tile set.

    INPUT:

    - ``tiles`` -- list of tiles, a tile is a 4-tuple (right color, top
        color, left color, bottom color)

    EXAMPLES::

        sage: from slabbe import WangTileSet
        sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
        sage: T = WangTileSet(tiles)

    ::

        sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
        ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
        sage: T = WangTileSet(tiles)
    """
    def __init__(self, tiles):
        r"""
        See documentation of the class.

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: T = WangTileSet(tiles)
        """
        self._tiles = tiles

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: T = WangTileSet(tiles)
            sage: next(iter(T))
            (0, 0, 0, 0)
        """
        return iter(self._tiles)

    def __len__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: T = WangTileSet(tiles)
            sage: len(T)
            3
        """
        return len(self._tiles)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: T = WangTileSet(tiles)
            sage: T
            Wang tile set of cardinality 3
        """
        return r"Wang tile set of cardinality {}".format(len(self))

    def tiles(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: T = WangTileSet(tiles)
            sage: T.tiles()
            [(0, 0, 0, 0), (1, 1, 1, 1), (2, 2, 2, 2)]
        """
        return self._tiles

class WangTileSet(WangTileSet_generic):
    r"""
    Construct a Wang tile set.

    INPUT:

    - ``tiles`` -- list of tiles, a tile is a 4-tuple (right color, top
        color, left color, bottom color)

    EXAMPLES::

        sage: from slabbe import WangTileSet
        sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
        sage: T = WangTileSet(tiles)

    ::

        sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
        ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
        sage: T = WangTileSet(tiles)
    """
    def vertical_alphabet(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,1,2,3), (4,5,6,7), (8,9,10,11)]
            sage: T = WangTileSet(tiles)
            sage: T.vertical_alphabet()
            {0, 2, 4, 6, 8, 10}
        """
        right, top, left, bottom = zip(*self._tiles)
        return set(left) | set(right)

    def horizontal_alphabet(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,1,2,3), (4,5,6,7), (8,9,10,11)]
            sage: T = WangTileSet(tiles)
            sage: T.horizontal_alphabet()
            {1, 3, 5, 7, 9, 11}
        """
        right, top, left, bottom = zip(*self._tiles)
        return set(top) | set(bottom)

    def dual(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
            ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
            sage: T = WangTileSet(tiles)
            sage: dual = T.dual()
            sage: dual
            Wang tile set of cardinality 7
            sage: dual.tiles()
            [(0, 0, 2, 0),
             (0, 1, 1, 0),
             (1, 2, 0, 0),
             (0, 0, 0, 1),
             (2, 1, 1, 1),
             (1, 1, 0, 2),
             (0, 2, 1, 2)]
        """
        tiles = []
        for (right, top, left, bottom) in self:
            tile = (top, right, bottom, left)
            tiles.append(tile)
        return WangTileSet(tiles)

    def to_transducer(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
            ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
            sage: T = WangTileSet(tiles)
            sage: T.to_transducer()
            Transducer with 3 states
        """
        from sage.combinat.finite_state_machine import Transducer
        transitions = []
        for (right, top, left, bottom) in self:
            # transition = (stateA, stateB, word_in, word_out)
            transition = (left, right, bottom, top)
            transitions.append(transition)
        return Transducer(transitions)

    def to_transducer_graph(self, label_function=tuple):
        r"""
        Return the graph of the transducer.

        Labels are cleaned. Label of multiedges are merged with commas.

        INPUT:

        - ``label_function`` -- function (default:``tuple``), a function to
          apply to each list of labels when merging multiedges into one

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = ['ABCD', 'EFGH', 'AXCY']
            sage: tiles = map(tuple, tiles)
            sage: T = WangTileSet(tiles)
            sage: G = T.to_transducer_graph()
            sage: G
            Digraph on 4 vertices

        The edge labels are clean::

            sage: G.edges()
            [('C', 'A', ('D|B', 'Y|X')), ('G', 'E', ('H|F',))]

        Using ``label_function``::

            sage: fn = lambda L: ','.join(map(str, L))
            sage: G = T.to_transducer_graph(label_function=fn)
            sage: G.edges()
            [('C', 'A', 'D|B,Y|X'), ('G', 'E', 'H|F')]

        Using ``label_function`` with latex expressions::

            sage: fn = lambda L: LatexExpr(','.join(map(str, L)))
            sage: G = T.to_transducer_graph(label_function=fn)
            sage: G.edges()
            [('C', 'A', D|B,Y|X), ('G', 'E', H|F)]

        This is to compared to::

            sage: T.to_transducer().graph().edges()
            [('C', 'A', "'D'|'B'"), ('C', 'A', "'Y'|'X'"), ('G', 'E', "'H'|'F'")]

        It works for integers entries::

            sage: tiles = [(0,1,2,3), (0,5,2,3)]
            sage: T = WangTileSet(tiles)
            sage: G = T.to_transducer_graph()
            sage: G
            Digraph on 2 vertices
            sage: G.edges()
            [(2, 0, ('3|1', '3|5'))]
        """
        from slabbe.graph import merge_multiedges
        def edge_labels(t):
            assert len(t.word_in) == 1, "we assume word_in is of length 1"
            assert len(t.word_out) == 1, "we assume word_out is of length 1"
            return "{}|{}".format(t.word_in[0], t.word_out[0])
        G = self.to_transducer().graph(edge_labels)
        return merge_multiedges(G, label_function)

    def system_of_density_equations(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
            ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
            sage: T = WangTileSet(tiles)
            sage: M = T.system_of_density_equations()
            sage: M
            [ 0  1  1 -1  0  0  0  0]
            [ 0 -1  0  1  0 -1  0  0]
            [ 0  0 -1  0  0  1  0  0]
            [ 1  1 -1  0  0 -1  1  0]
            [ 0 -1  1  0 -1  1 -1  0]
            [-1  0  0  0  1  0  0  0]
            [ 1  1  1  1  1  1  1  1]
            sage: M.rank()
            5
        """
        from sage.modules.free_module import FreeModule
        from sage.rings.integer_ring import ZZ
        from sage.matrix.constructor import matrix
        M = FreeModule(ZZ, len(self)+1)
        vertical = defaultdict(M)
        horizontal = defaultdict(M)
        for i, (right, top, left, bottom) in enumerate(self):
            vertical[left][i] += 1
            vertical[right][i] -= 1
            horizontal[top][i] += 1
            horizontal[bottom][i] -= 1
        rows = vertical.values()
        rows.extend(horizontal.values())
        rows.append(M([1 for _ in range(len(self)+1)]))
        return matrix(rows)

    def polyhedron_of_densities(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
            ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
            sage: T = WangTileSet(tiles)
            sage: P = T.polyhedron_of_densities()
            sage: P
            A 2-dimensional polyhedron in QQ^7 defined as the convex hull of 3 vertices
            sage: P.vertices()
            (A vertex at (0, 2/7, 1/7, 3/7, 0, 1/7, 0),
             A vertex at (0, 0, 1/5, 1/5, 0, 1/5, 2/5),
             A vertex at (2/7, 0, 1/7, 1/7, 2/7, 1/7, 0))
            sage: (0, 0, 1/5, 1/5, 0, 1/5, 2/5) in P
            True

        Jeandel-Rao tiles::

            sage: tiles = [(2,4,2,1), (2,2,2,0), (1,1,3,1), (1,2,3,2), (3,1,3,3),
            ....: (0,1,3,1), (0,0,0,1), (3,1,0,2), (0,2,1,2), (1,2,1,4), (3,3,1,2)]
            sage: T = WangTileSet(tiles)
            sage: P = T.polyhedron_of_densities()
            sage: P
            A 4-dimensional polyhedron in QQ^11 defined as the convex hull of 10 vertices
            sage: P.vertices()
            (A vertex at (0, 1/5, 1/5, 0, 1/5, 0, 1/5, 0, 0, 0, 1/5),
             A vertex at (0, 1/5, 0, 1/5, 1/5, 0, 1/5, 0, 0, 0, 1/5),
             A vertex at (0, 1/5, 1/5, 0, 0, 0, 1/5, 1/5, 1/5, 0, 0),
             A vertex at (0, 1/5, 0, 1/5, 0, 0, 1/5, 1/5, 1/5, 0, 0),
             A vertex at (0, 1/4, 0, 0, 0, 1/4, 1/4, 1/4, 0, 0, 0),
             A vertex at (1/4, 0, 0, 0, 0, 1/4, 0, 1/4, 0, 1/4, 0),
             A vertex at (1/5, 0, 1/5, 0, 1/5, 0, 0, 0, 0, 1/5, 1/5),
             A vertex at (1/5, 0, 1/5, 0, 0, 0, 0, 1/5, 1/5, 1/5, 0),
             A vertex at (1/5, 0, 0, 1/5, 1/5, 0, 0, 0, 0, 1/5, 1/5),
             A vertex at (1/5, 0, 0, 1/5, 0, 0, 0, 1/5, 1/5, 1/5, 0))
        """
        from sage.modules.free_module import FreeModule
        from sage.rings.integer_ring import ZZ
        from sage.geometry.polyhedron.constructor import Polyhedron
        M = FreeModule(ZZ, len(self)+1)
        vertical = defaultdict(M)
        horizontal = defaultdict(M)
        for i, (right, top, left, bottom) in enumerate(self):
            vertical[left][i+1] += 1
            vertical[right][i+1] -= 1
            horizontal[top][i+1] += 1
            horizontal[bottom][i+1] -= 1
        eqns = vertical.values()
        eqns.extend(horizontal.values())
        v = M([1 for _ in range(len(self)+1)])
        v[0] = -1
        eqns.append(v)
        # An entry equal to "[-1,7,3,4]" represents the equality
        # 7x_1+3x_2+4x_3= 1.
        ieqs = []
        for i in range(len(self)):
            v = M()
            v[i+1] = 1
            ieqs.append(v)
        return Polyhedron(ieqs=ieqs, eqns=eqns)

    @rename_keyword(fontsize='font')
    def tikz(self, ncolumns=10, color=None, size=1, space=.1, scale=1,
             font=r'\normalsize', rotate=None, label_shift=.2):
        r"""
        INPUT:

        - ``ncolumns`` -- integer (default: ``10``)
        - ``color`` -- dict (default: None)
        - ``size`` -- number (default: ``1``)
        - ``space`` -- number (default: ``.1``)
        - ``scale`` -- number (default: ``1``)
        - ``font`` -- string (default: ``r'\normalsize'``)
        - ``rotate`` -- list or ``None`` (default:``None``) list of four angles
          in degrees like ``(0,0,0,0)``, the rotation angle to apply to each
          label of Wang tiles. If ``None``, it performs a 90 degres rotation
          for left and right labels taking more than one character.
        - ``label_shift`` -- number (default: ``.2``) translation distance
          of the label from the edge

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
            ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
            sage: T = WangTileSet(tiles)
            sage: color = {0:'white',1:'red',2:'cyan',3:'green',4:'white'}
            sage: _ = T.tikz(color=color).pdf(view=False)
        """
        from slabbe import TikzPicture
        lines = []
        lines.append(r'\begin{tikzpicture}')
        lines.append('[scale={}]'.format(scale))
        lines.append(r'\tikzstyle{{every node}}=[font={}]'.format(font))
        for i,tile in enumerate(self):
            x = i % ncolumns
            y = - (i // ncolumns)
            position = (x * (size + space), y * (size + space))
            new_lines = tile_to_tikz(tile, position, color=color,
                    size=size, rotate=rotate, label_shift=label_shift,
                    top_right_edges=True)
            lines.extend(new_lines)
        lines.append(r'\end{tikzpicture}')
        return TikzPicture('\n'.join(lines))

    def create_tikz_pdf_files(self, prefix='tile', color=None):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
            ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
            sage: T = WangTileSet(tiles)
            sage: T.create_tikz_pdf_files() # not tested

        This creates tile0.pdf, tile1.pdf, etc. in the repository.
        """
        if color is None:
            color = {0:'white',1:'red',2:'cyan',3:'green',4:'white'}
        for i in range(11):
            tiling = WangTiling([[i]], tiles)
            tiling.tikz().pdf('{}-{}.pdf'.format(prefix,i))
            tiling.tikz(color).pdf('{}-{}_colored.pdf'.format(prefix,i))

    @rename_keyword(fontsize='font')
    def create_macro_file(self, filename='macro.tex', command_name='Tile',
            color=None, size=1, scale=1, font=r'\normalsize',
            rotate=None, label_shift=.2):
        r"""
        INPUT:

        - ``filename`` -- string (default: ``r'macro.tex'``)
        - ``comand_name`` -- string (default: ``r'Tile'``)
        - ``color`` -- dict (default: None)
        - ``size`` -- number (default: ``1``)
        - ``scale`` -- number (default: ``1``)
        - ``font`` -- string (default: ``r'\normalsize'``)
        - ``rotate`` -- list or ``None`` (default:``None``) list of four angles
          in degrees like ``(0,0,0,0)``, the rotation angle to apply to each
          label of Wang tiles. If ``None``, it performs a 90 degres rotation
          for left and right labels taking more than one character.
        - ``label_shift`` -- number (default: ``.2``) translation distance
          of the label from the edge

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
            ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
            sage: T = WangTileSet(tiles)
            sage: T.create_macro_file() # not tested
            creation of file macro.tex

        ::

            sage: color = {0:'white',1:'red',2:'cyan',3:'green',4:'white'}
            sage: T.create_macro_file(color=color) # not tested
            creation of file macro.tex
        """
        from roman import toRoman
        lines = []
        for i,tile in enumerate(self):
            roman_i = toRoman(i) if i!=0 else 'O'
            lines.append(r'\newcommand\{}{}{{'.format(command_name, roman_i))
            lines.append(r'\begin{tikzpicture}')
            lines.append('[scale={}]'.format(scale))
            lines.append(r'\tikzstyle{{every node}}=[font={}]'.format(font))
            new_lines = tile_to_tikz(tile, position=(0,0), color=color,
                    size=size, rotate=rotate, label_shift=label_shift,
                    top_right_edges=True)
            lines.extend(new_lines)
            lines.append(r'\end{tikzpicture}')
            lines.append(r'} % end of newcommand')
        s = '\n'.join(lines)
        with open(filename, 'w') as f:
            f.write(s)
            print("creation of file {}".format(filename))

    @rename_keyword(fontsize='font')
    def substitution_tikz(self, substitution, function=None, color=None,
            size=1, scale=1, font=r'\normalsize', rotate=None,
            label_shift=.2, transformation_matrix=None, ncolumns=4,
            tabular='tabular', align='l'):
        r"""
        Return the tikz code showing what the substitution A->B* does on
        Wang tiles.

        Note: we assume that the tiles in self are the elements of B.

        INPUT:

        - ``substitution`` -- substitution 2d
        - ``fn`` -- a function (default: ``None``) to apply to the
          new colors which are tuple of previous colors
        - ``color`` -- dict (default: ``None``) from tile values -> tikz colors
        - ``size`` -- number (default: ``1``), size of the tile
        - ``scale`` -- number (default: ``1``), scale of tikzpicture
        - ``font`` -- string (default: ``r'\normalsize'``
        - ``rotate`` -- list or ``None`` (default:``None``) list of four angles
          in degrees like ``(0,0,0,0)``, the rotation angle to apply to each
          label of Wang tiles. If ``None``, it performs a 90 degres rotation
          for left and right labels taking more than one character.
        - ``label_shift`` -- number (default: ``.2``) translation distance
          of the label from the edge
        - ``transformation_matrix`` -- matrix (default: ``None``), a matrix
          to apply to the coordinate before drawing, it can be in
          ``SL(2,ZZ)`` or not.
        - ``ncolumns`` -- integer (default: ``4``)
        - ``tabular`` -- string (default: ``'tabular'``) or ``'longtable'``
        - ``align`` -- character (default:``'l'``), latex alignment symbol
          ``'l'``, ``'r'`` or ``'c'``.

        OUTPUT:

            dict, key -> tile
        """
        raise NotImplementedError("use method of the substitution instead")

    def desubstitute(self, substitution, function=None):
        r"""
        Return the Wang tile set obtained from the desubstitution of the
        given Wang tile set.

        INPUT:

        - ``substitution`` -- substitution 2d
        - ``fn`` -- a function (default: ``None``) to apply to the
          new colors which are tuple of previous colors

        OUTPUT:

            dict, key -> tile

        EXAMPLES::

            sage: from slabbe import Substitution2d, WangTileSet
            sage: A = [[0,1,2],[1,0,0]]
            sage: B = [[0,1,2]]
            sage: d = {4:A, 5:B}
            sage: s = Substitution2d(d)
            sage: tiles = [(0,3,1,4), (1,4,0,3), (5,6,7,8)]
            sage: W = WangTileSet(tiles)
            sage: W.desubstitute(s)
            {4: ((1, 0, 0), (6, 3), (1, 0, 7), (4, 3)),
             5: ((0, 1, 5), (6,), (1, 0, 7), (4,))}

        Providing a function which gets back to integers::

            sage: fn = lambda colors:int(''.join(map(str, colors)))
            sage: W.desubstitute(s, fn)
            {4: (100, 63, 107, 43), 5: (15, 6, 107, 4)}

        Providing a function which concatenate label as strings::

            sage: fn = lambda colors:''.join(map(str, colors))
            sage: W.desubstitute(s, fn)
            {4: ('100', '63', '107', '43'), 5: ('015', '6', '107', '4')}
        """
        if function is None:
            function = lambda x:x
        d = {}
        for a,image_a in substitution._d.items():
            # get the border tiles
            west = image_a[0]
            east = image_a[-1]
            north = [column[-1] for column in image_a]
            south = [column[0] for column in image_a]
            # get the good color for each
            west = tuple(self._tiles[b][2] for b in west)
            east = tuple(self._tiles[b][0] for b in east)
            north = tuple(self._tiles[b][1] for b in north)
            south = tuple(self._tiles[b][3] for b in south)
            # create the tile and save
            tile = tuple(function(color) for color in (east, north, west, south))
            d[a] = tile
        return d


    def admissible_horizontal_words(self, length, width, height):
        r"""
        Return the horizontal word of given length appearing in every
        position inside a rectangle of given width and height.

        INPUT:

        - ``length`` -- integer
        - ``width`` -- integer
        - ``height`` -- integer

        OUTPUT:

            set of tuples

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: T = WangTileSet(tiles)
            sage: T.admissible_horizontal_words(2,2,2)
            {(0, 0), (1, 1), (2, 2)}

        The horizontal word 22 is impossible after looking at large enough
        boxes::

            sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
            ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
            sage: T = WangTileSet(tiles)
            sage: T.admissible_horizontal_words(2,2,2)
            {(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0), (2, 2)}
            sage: T.admissible_horizontal_words(2,3,3)
            {(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0)}
            sage: T.admissible_horizontal_words(2,4,4)
            {(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0)}
            sage: T.admissible_horizontal_words(2,5,5)
            {(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0)}

        """
        wang_tile_solver = WangTileSolver(self._tiles, width, height)
        seen = defaultdict(set)
        for tiling in wang_tile_solver.solutions_iterator():
            for pos,word in tiling.horizontal_words_dict(length).items():
                seen[pos].add(word)
        return set.intersection(*seen.values())

    def admissible_vertical_words(self, length, width, height):
        r"""
        Return the vertical word of given length appearing in every
        position inside a rectangle of given width and height.

        INPUT:

        - ``length`` -- integer
        - ``width`` -- integer
        - ``height`` -- integer

        OUTPUT:

            set of tuples

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: T = WangTileSet(tiles)
            sage: T.admissible_vertical_words(2,2,2)
            {(0, 0), (1, 1), (2, 2)}

        Every word of length 2 appear as a vertical word in every position
        of a `5\times 5` box::

            sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
            ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
            sage: T = WangTileSet(tiles)
            sage: T.admissible_vertical_words(2,2,2)
            {(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)}
            sage: T.admissible_vertical_words(2,5,5)
            {(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)}

        """
        wang_tile_solver = WangTileSolver(self._tiles, width, height)
        seen = defaultdict(set)
        for tiling in wang_tile_solver.solutions_iterator():
            for pos,word in tiling.vertical_words_dict(length).items():
                seen[pos].add(word)
        return set.intersection(*seen.values())

    def fusion(self, other, direction, function=str.__add__, initial=''):
        r"""
        Return the fusion of wang tile sets in the given direction.

        We keep only the strongly connected components.

        INPUT:

        - ``other`` -- WangTileSet
        - ``direction`` -- integer (1 or 2)
        - ``function`` -- function (default:``str.__add__``), monoid
          operation
        - ``initial`` -- object (default:``''``), monoid neutral

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = ['ABCD', 'EFGH', 'AXCY', 'ABAB']
            sage: tiles = map(tuple, tiles)
            sage: T = WangTileSet(tiles)
            sage: T1T = T.fusion(T, 1)
            sage: T1T.tiles()
            [('A', 'BB', 'A', 'BB')]
            sage: T2T = T.fusion(T, 2)
            sage: T2T.tiles()
            [('AA', 'B', 'AA', 'B')]

        To keep integers, one way is to wrap them into a tuple and do tuple
        operations::

            sage: tiles = [(0,1,0,1)]
            sage: tiles = [tuple((a,) for a in tile) for tile in tiles]
            sage: T = WangTileSet(tiles)
            sage: T2T = T.fusion(T, 2, function=tuple.__add__, initial=tuple())
            sage: T2T2T = T2T.fusion(T, 2, function=tuple.__add__, initial=tuple())
            sage: T2T2T.tiles()
            [((0, 0, 0), (1,), (0, 0, 0), (1,))]

        """
        if not isinstance(other, WangTileSet):
            raise TypeError('other(={}) must be a'
                    ' WangTileSet'.format(other))
        if direction == 1:
            return self.dual().fusion(other.dual(), direction=2,
                    function=function, initial=initial).dual()
        if not direction == 2:
            raise ValueError('direction(={}) must be 1 or 2'.format(direction))

        T = self.to_transducer()
        U = other.to_transducer()
        TU = T.composition(U)

        #edge_labels = lambda t:"{}|{}".format(t.word_in, t.word_out)
        edge_labels = lambda t:(t.word_in, t.word_out)
        G = TU.graph(edge_labels)
        from slabbe.graph import clean_sources_and_sinks
        H = clean_sources_and_sinks(G)

        tiles = []
        for (u,v,(word_in, word_out)) in H.edges():
            left = reduce(function, (a.label() for a in u), initial)
            right = reduce(function, (a.label() for a in v), initial)
            bottom = reduce(function, word_in, initial)
            top = reduce(function, word_out, initial)
            tiles.append((right, top, left, bottom))
        return WangTileSet(tiles)

    def solver(self, width, height, preassigned_color=None,
            preassigned_tiles=None, color=None):
        r"""
        Return the Wang tile solver of this Wang tile set inside a
        rectangle of given width and height.

        INPUT:

        - ``width`` -- integer
        - ``height`` -- integer
        - ``preassigned_color`` -- None or list of 4 dict or the form ``[{},
          {}, {}, {}]`` right, top, left, bottom colors preassigned to some
          positions (on the border or inside)
        - ``preassigned_tiles`` -- None or dict of tiles preassigned to some
          positions
        - ``color`` -- None or dict

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
            sage: T = WangTileSet(tiles)
            sage: W = T.solver(3,3)
            sage: W.solve()
            A wang tiling of a 3 x 3 rectangle

        ::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: T = WangTileSet(tiles)
            sage: W = T.solver(3,3, preassigned_tiles={(1,1):0})
            sage: W.solve().table()
            [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
            sage: W = T.solver(3,3, preassigned_tiles={(1,1):1})
            sage: W.solve().table()
            [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
            sage: W = T.solver(3,3, preassigned_tiles={(1,1):2})
            sage: W.solve().table()
            [[2, 2, 2], [2, 2, 2], [2, 2, 2]]

        TESTS::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
            sage: T = WangTileSet(tiles)
            sage: W = T.solver(3,3, preassigned_tiles={(1,1):3})
            sage: W.solve().table()
            Traceback (most recent call last):
            ...
            MIPSolverException: ...

        """
        return WangTileSolver(self._tiles, width, height,
                preassigned_color=preassigned_color,
                preassigned_tiles=preassigned_tiles,
                color=color)

    def tiles_allowing_surrounding(self, radius, solver=None, verbose=False):
        r"""
        Return the subset of tiles allowing a surrounding of given radius.

        INPUT:

        - ``radius`` - integer
        - ``solver`` - string or None
        - ``verbose`` - boolean

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
            sage: T = WangTileSet(tiles)
            sage: U = T.tiles_allowing_surrounding(1)
            sage: U
            Wang tile set of cardinality 3
            sage: U.tiles()
            [(0, 0, 0, 0), (1, 1, 1, 1), (2, 2, 2, 2)]

        ::

            sage: T.tiles_allowing_surrounding(1, verbose=True)
            Solution found for tile 0:
            [[0, 0, 3], [0, 0, 0], [0, 0, 0]]
            Solution found for tile 1:
            [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
            Solution found for tile 2:
            [[2, 2, 2], [2, 2, 2], [2, 2, 2]]
            Wang tile set of cardinality 3

        """
        diameter = 2*radius+1
        L = []
        for i,t in enumerate(self):
            d = {(radius,radius):i}
            s = self.solver(diameter, diameter, preassigned_tiles=d)
            if s.has_solution(solver=solver):
                if verbose:
                    print("Solution found for tile {}:\n{}".format(i,
                                s.solve(solver)._table))
                L.append(t)
        return WangTileSet(L)

    def not_forbidden_tilings(self, width, height, radius=1,
            function=str.__add__, initial='', verbose=False):
        r"""
        Return the set of valid of a rectangle of given width and height
        allowing a surrounding of itself of given radius.

        INPUT:

        - ``width`` - integer
        - ``height`` - integer
        - ``radius`` - integer
        - ``function`` -- function (default:``str.__add__``), monoid
          operation for fusion of colors
        - ``initial`` -- object (default:``''``), monoid neutral, for
          fusion of colors
        - ``verbose`` - boolean

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
            sage: tiles = [map(str, tile) for tile in tiles]
            sage: T = WangTileSet(tiles)
            sage: S = T.not_forbidden_tilings(2,2)
            sage: S
            [A wang tiling of a 2 x 2 rectangle,
             A wang tiling of a 2 x 2 rectangle,
             A wang tiling of a 2 x 2 rectangle]
            sage: [a.table() for a in S]
            [[[0, 0], [0, 0]], [[1, 1], [1, 1]], [[2, 2], [2, 2]]]

        ::

            sage: S = T.not_forbidden_tilings(3,3)
            sage: S
            [A wang tiling of a 3 x 3 rectangle,
             A wang tiling of a 3 x 3 rectangle,
             A wang tiling of a 3 x 3 rectangle]
            sage: [a.table() for a in S]
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
             [[2, 2, 2], [2, 2, 2], [2, 2, 2]]]

        """
        T = base = self
        for _ in range(width-1):
            T = T.fusion(base, 1, function=function, initial=initial)
        if verbose:
            print("After fusion in the direction e1: ", T)
        base = T
        for _ in range(height-1):
            T = T.fusion(base, 2, function=function, initial=initial)
        if verbose:
            print("After fusion in the direction e2: ", T)
        T = T.tiles_allowing_surrounding(1)
        if verbose:
            print("After filtering tiles without surrounding of "
                  "radius {} : {}".format(radius, T))
        L = []
        for t in T:
            right = {(width-1,i):a for (i,a) in enumerate(t[0])}
            top = {(i,height-1):a for (i,a) in enumerate(t[1])}
            left = {(0,i):a for (i,a) in enumerate(t[2])}
            bottom = {(i,0):a for (i,a) in enumerate(t[3])}
            preassigned_color=[right,top,left,bottom]
            W = self.solver(width, height,
                    preassigned_color=preassigned_color)
            L.extend(W.solutions_iterator())
        return L

class HexagonalWangTileSet(WangTileSet_generic):
    r"""
    Construct an hexagonal Wang tile set.

    INPUT:

    - ``tiles`` -- list of tiles, a tile is a 6-tuple of colors (right, top
      right, top left, left, bottom left, bottom right)

    EXAMPLES::

        sage: from slabbe import HexagonalWangTileSet
        sage: tiles = [(0,0,0,0,0,0), (1,1,1,1,1,1), (2,2,2,2,2,2)]
        sage: T = HexagonalWangTileSet(tiles)
    """
    @rename_keyword(fontsize='font')
    def tikz(self, ncolumns=10, color=None, radius=1, space=.1, scale=1,
             font=r'\normalsize', rotate=(0,0,0,0,0,0), label_shift=.2):
        r"""
        INPUT:

        - ``ncolumns`` -- integer (default: ``10``)
        - ``color`` -- dict (default: None)
        - ``radius`` -- number (default: ``1``)
        - ``space`` -- number (default: ``.1``)
        - ``scale`` -- number (default: ``1``)
        - ``font`` -- string (default: ``r'\normalsize'``
        - ``rotate`` -- list (default:``(0,0,0,0,0,0)``) of four angles in
          degrees, the rotation angle to apply to each label of Wang tiles
        - ``label_shift`` -- number (default: ``.2``) translation distance
          of the label from the edge

        EXAMPLES::

            sage: from slabbe import HexagonalWangTileSet
            sage: tiles = [(0,0,0,0,0,2), (1,0,1,0,0,1), (1,0,2,1,0,0),
            ....:          (1,2,2,1,1,1), (1,1,1,1,2,0), (0,2,2,0,2,1)]
            sage: T = HexagonalWangTileSet(tiles)
            sage: color = {0:'white',1:'red',2:'cyan',3:'green',4:'white'}
            sage: _ = T.tikz(color=color).pdf(view=False)

        TESTS::

            sage: _ = T.tikz(color=color,label_shift=.5).pdf(view=False)
        """
        from slabbe import TikzPicture
        lines = []
        lines.append(r'\begin{tikzpicture}')
        lines.append('[scale={}]'.format(scale))
        lines.append(r'\tikzstyle{{every node}}=[font={}]'.format(font))
        for i,tile in enumerate(self):
            x = i % ncolumns
            y = - (i // ncolumns)
            position = (x * (2*radius + space), y * (2*radius + space))
            new_lines = hexagonal_tile_to_tikz(tile, position, color=color,
                    radius=radius, rotate=rotate, label_shift=label_shift)
            lines.extend(new_lines)
        lines.append(r'\end{tikzpicture}')
        return TikzPicture('\n'.join(lines))

class WangTileSolver(object):
    r"""
    Wang tile solver inside a rectangle of given width and height.

    INPUT:

    - ``tiles`` -- list of tiles, a tile is a 4-tuple (right color, top
      color, left color, bottom color)
    - ``width`` -- integer
    - ``height`` -- integer
    - ``preassigned_color`` -- None or list of 4 dict or the form ``[{},
      {}, {}, {}]`` right, top, left, bottom colors preassigned to some
      positions (on the border or inside)
    - ``preassigned_tiles`` -- None or dict of tiles preassigned to some
      positions
    - ``color`` -- None or dict

    EXAMPLES::

        sage: from slabbe import WangTileSolver
        sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
        sage: W = WangTileSolver(tiles, 3, 3)
        sage: tiling = W.solve()
        sage: tiling._table
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    With color 2 preassigned to the right part of tile at position (1,1)::

        sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
        sage: right = {(1,1):2}
        sage: W = WangTileSolver(tiles,3,3,preassigned_color=[right,{},{},{}])
        sage: tiling = W.solve()
        sage: tiling._table
        [[2, 2, 2], [2, 2, 2], [2, 2, 2]]

    With tile 2 preassigned at position (0,1)::

        sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
        sage: preassigned = {(0,1):1}
        sage: W = WangTileSolver(tiles,3,3,preassigned_tiles=preassigned)
        sage: tiling = W.solve()
        sage: tiling._table
        [[1, 1, 1], [1, 1, 1], [1, 1, 1]]

    When constraints are inconsistent::

        sage: right = {(1,1):1, (2,2):0}
        sage: W = WangTileSolver(tiles,3,3,preassigned_color=[right,{},{},{}])
        sage: W.solve(solver='GLPK')
        Traceback (most recent call last):
        ...
        MIPSolverException: GLPK: Problem has no feasible solution

    TESTS::

        sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
        sage: preassigned = {(0,1):1}
        sage: W = WangTileSolver(tiles,3,3,preassigned_tiles=preassigned)
        sage: tiling = W.solve()
        sage: tiling._table
        [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
    """
    def __init__(self, tiles, width, height, preassigned_color=None,
            preassigned_tiles=None, color=None):
        r"""
        See class for documentation.

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles, 3, 4)
        """
        self._tiles = tiles
        self._width = width
        self._height = height
        if preassigned_color is None:
            preassigned_color = [{}, {}, {}, {}]
        self._preassigned_color = preassigned_color
        if preassigned_tiles is None:
            preassigned_tiles = {}
        self._preassigned_tiles = preassigned_tiles
        self._color = color

    def vertical_alphabet(self):
        right, top, left, bottom = zip(*self._tiles)
        return set(left) | set(right)

    def horizontal_alphabet(self):
        right, top, left, bottom = zip(*self._tiles)
        return set(top) | set(bottom)

    @cached_method
    def milp(self, solver=None):
        r"""
        Return the Mixed integer linear program.

        INPUT:

        - ``solver`` -- string or None (default: ``None``), other possible
          values are ``'Coin'`` or ``'Gurobi'``

        OUTPUT:

        a tuple (p,x) where p is the MILP and x is the variable

        .. NOTE::

            In some cases, calling this method takes much more time (few
            minutes) than calling the method ``solve`` which takes few
            seconds.

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: p,x = W.milp(solver='GLPK')
            sage: p
            Boolean Program (maximization, 36 variables, 29 constraints)
            sage: x
            MIPVariable of dimension 1

        Then you can solve it and get the solutions::

            sage: p.solve()
            1.0
            sage: soln = p.get_values(x)
            sage: support = [key for key in soln if soln[key]]
            sage: support
            [(0, 1, 1), (0, 1, 3), (0, 2, 1), (0, 2, 0), (0, 2, 3), (0, 2, 2), 
             (0, 1, 2), (0, 0, 3), (0, 0, 2), (0, 0, 1), (0, 0, 0), (0, 1, 0)]

        Other solver can be used::

            sage: p,x = W.milp(solver='Gurobi')   # optional gurobi

        TESTS:

        Colors do not have to be integers::

            sage: tiles = [('a','a','a','a'), ('b','b','b','b')]
            sage: W = WangTileSolver(tiles,3,4)
            sage: p,x = W.milp()
            sage: tiling = W.solve()
        """
        tiles = self._tiles
        indices = range(len(tiles))

        # x[i,j,k] == 1 iff tile i is at position (j,k)
        p = MixedIntegerLinearProgram(solver=solver)
        x = p.new_variable(binary=True)

        # exactly one tile at each position (j,k)
        for j in range(self._width):
            for k in range(self._height):
                S = p.sum(x[i,j,k] for i in indices)
                name = "one tile at {}".format((j,k))
                p.add_constraint(S==1, name=name)

        # preassigned tiles at position (j,k)
        for j,k in self._preassigned_tiles:
            i = self._preassigned_tiles[(j,k)]
            name = "preassigned tile {} at {}".format(i, (j,k))
            p.add_constraint(x[i,j,k]==1, name=name)

        # matching vertical colors
        Va = sorted(self.vertical_alphabet())
        Va_to_int = {a:i for i,a in enumerate(Va)}
        for j in range(self._width-1):
            for k in range(self._height):
                A = p.sum(Va_to_int[tiles[i][0]]*x[i,j,k] for i in indices)
                B = p.sum(Va_to_int[tiles[i][2]]*x[i,j+1,k] for i in indices)
                name = "matching right of {}".format((j,k))
                p.add_constraint(A==B, name=name)

        # matching horizontal colors
        Ha = sorted(self.horizontal_alphabet())
        Ha_to_int = {a:i for i,a in enumerate(Ha)}
        for j in range(self._width):
            for k in range(self._height-1):
                A = p.sum(Ha_to_int[tiles[i][1]]*x[i,j,k] for i in indices)
                B = p.sum(Ha_to_int[tiles[i][3]]*x[i,j,k+1] for i in indices)
                name = "matching top of {}".format((j,k))
                p.add_constraint(A==B, name=name)

        # matching preassigned color constraints
        legend = {0:'right',1:'top',2:'left',3:'bottom'}
        for angle, D in enumerate(self._preassigned_color):
            if angle == 0 or angle == 2:
                to_int = Va_to_int
            else:
                to_int = Ha_to_int
            for j,k in D:
                A = p.sum(to_int[tiles[i][angle]]*x[i,j,k] for i in indices)
                name = "preassigned color {} of {}".format(legend[angle], (j,k))
                p.add_constraint(A==to_int[D[(j,k)]], name=name)

        p.set_objective(x[0,0,0])
        return p, x

    def solve(self, solver=None, solver_parameters=None, ncpus=1):
        r"""
        Return a dictionary associating to each tile a list of positions
        where to find this tile.

        INPUT:

        - ``solver`` -- string or None (default: ``None``), 
          ``'dancing_links'`` or the name of a MILP solver in Sage like
          ``'GLPK'``, ``'Coin'`` or ``'Gurobi'``.
        - ``solver_parameters`` -- dict (default: ``{}``), parameters given
          to the MILP solver using method ``solver_parameter``. For a list
          of available parameters for example for the Gurobi backend, see
          dictionary ``parameters_type`` in the file
          ``sage/numerical/backends/gurobi_backend.pyx``
        - ``ncpus`` -- integer (default: ``1``), maximal number of
          subprocesses to use at the same time, used only if ``solver`` is
          ``'dancing_links'``.

        OUTPUT:

            a wang tiling object

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: table = tiling._table
            sage: table
            [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]

        The tile at position (1,3) is::

            sage: table[1][3]
            0

        Allowing more threads while using Gurobi::

            sage: W = WangTileSolver(tiles,3,4)
            sage: kwds = dict(Threads=4)
            sage: tiling = W.solve(solver='Gurobi', kwds) # optional Gurobi
            sage: tiling._table                           # optional Gurobi
            [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]

        Using dancing links::

            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve(solver='dancing_links', ncpus=8)
            sage: tiling
            A wang tiling of a 3 x 4 rectangle

        Using dancing links with tile 2 preassigned at position (0,1)::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: preassigned = {(0,1):1}
            sage: W = WangTileSolver(tiles,3,3,preassigned_tiles=preassigned)
            sage: tiling = W.solve(solver='dancing_links')
            sage: tiling._table
            [[1, 1, 1], [1, 1, 1], [1, 1, 1]]

        Using dancing links when constraints are inconsistent::

            sage: right = {(1,1):1, (2,2):0}
            sage: W = WangTileSolver(tiles,3,3,preassigned_color=[right,{},{},{}])
            sage: W.solve(solver='dancing_links')
            Traceback (most recent call last):
            ...
            ValueError: no solution found using dancing links

        REFERENCES:

            How do I set solver_parameter to make Gurobi use more than one
            processor?, https://ask.sagemath.org/question/37726/
        """
        if solver == 'dancing_links':
            from sage.combinat.matrices.dancing_links import dlx_solver
            rows,row_info = self.rows_and_information()
            dlx = dlx_solver(rows)
            solution = dlx.one_solution(ncpus=ncpus)
            if solution is None:
                raise ValueError('no solution found using dancing links')
            table = [[None]*self._height for _ in range(self._width)]
            for a in solution:
                j,k,i = row_info[a]
                assert table[j][k] is None, "table[{}][{}](={}) is not None".format(j,k,table[j][k])
                table[j][k] = i
            return WangTiling(table, self._tiles, color=self._color)
        else:
            p,x = self.milp(solver=solver)
            if solver_parameters is None:
                solver_parameters = {}
            for key, value in solver_parameters.items():
                p.solver_parameter(key, value)
            p.solve()
            soln = p.get_values(x)
            support = [key for key in soln if soln[key]]
            assert len(support) == self._width * self._height, ("len(support)={} "
                    "!= width*height".format(len(support)))
            table = [[None]*self._height for _ in range(self._width)]
            for i,j,k in support:
                table[j][k] = i
            return WangTiling(table, self._tiles, color=self._color)

    def has_solution(self, solver=None, solver_parameters=None, ncpus=1):
        r"""
        Return whether there is a solution.

        .. TODO::

            - Currently, the dancing links reduction ignores the
              preassigned parameters.

        INPUT:

        - ``solver`` -- string or None (default: ``None``),
          ``'dancing_links'`` or the name of a MILP solver in Sage like
          ``'GLPK'``, ``'Coin'`` or ``'Gurobi'``.
        - ``solver_parameters`` -- dict (default: ``{}``), parameters given
          to the MILP solver using method ``solver_parameter``. For a list
          of available parameters for example for the Gurobi backend, see
          dictionary ``parameters_type`` in the file
          ``sage/numerical/backends/gurobi_backend.pyx``
        - ``ncpus`` -- integer (default: ``1``), maximal number of
          subprocesses to use at the same time, used only if ``solver`` is
          ``'dancing_links'``.

        OUTPUT:

            a wang tiling object

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: W.has_solution()
            True

        Allowing more threads while using Gurobi::

            sage: W = WangTileSolver(tiles,3,4)
            sage: kwds = dict(Threads=4)
            sage: tiling = W.has_solution(solver='Gurobi', kwds) # optional Gurobi
            True

        Using dancing links::

            sage: W = WangTileSolver(tiles,3,4)
            sage: W.has_solution(solver='dancing_links', ncpus=8)
            True
        """
        if solver == 'dancing_links':
            from sage.combinat.matrices.dancing_links import dlx_solver
            rows,row_info = self.rows_and_information()
            dlx = dlx_solver(rows)
            solution = dlx.one_solution(ncpus=ncpus)
            return solution is not None
        else:
            p,x = self.milp(solver=solver)
            if solver_parameters is None:
                solver_parameters = {}
            for key, value in solver_parameters.items():
                p.solver_parameter(key, value)
            from sage.numerical.mip import MIPSolverException
            try:
                p.solve()
            except MIPSolverException:
                return False
            return True

    def rows_and_information(self, verbose=False):
        r"""
        Return the rows to give to the dancing links solver.

        INPUT:

        - ``verbose`` -- bool (default: ``False``)

        OUTPUT:

        Two lists:

            - the rows
            - row information (j,k,i) meaning tile i is at position (j,k)

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles, 4, 1)
            sage: rows,row_info = W.rows_and_information()
            sage: rows
            [[1, 2, 9],
             [0, 2, 9],
             [2, 9],
             [0, 4, 5, 10],
             [1, 3, 5, 10],
             [0, 1, 5, 10],
             [3, 7, 8, 11],
             [4, 6, 8, 11],
             [3, 4, 8, 11],
             [6, 12],
             [7, 12],
             [6, 7, 12]]
            sage: row_info
            [(0, 0, 0),
             (0, 0, 1),
             (0, 0, 2),
             (1, 0, 0),
             (1, 0, 1),
             (1, 0, 2),
             (2, 0, 0),
             (2, 0, 1),
             (2, 0, 2),
             (3, 0, 0),
             (3, 0, 1),
             (3, 0, 2)]
            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: dlx = dlx_solver(rows)
            sage: dlx
            Dancing links solver for 13 columns and 12 rows
            sage: dlx.search()
            1
            sage: dlx.get_solution()
            [1, 4, 7, 10]
            sage: row_info[1]
            (0, 0, 1)
            sage: row_info[4]
            (1, 0, 1)
            sage: row_info[7]
            (2, 0, 1)
            sage: row_info[10]
            (3, 0, 1)

        ... which means tile 1 is at position (0,0), (1,0), (2,0) and (3,0)

        TESTS::

            sage: tiles = [(0,0,0,0), (1,1,1,1)]
            sage: W = WangTileSolver(tiles, 4, 1)
            sage: W.rows_and_information(verbose=True)
            Vertical colors (coded using 3 bits):
            color 0 represented by bits [0] when on left
            color 0 represented by bits [1, 2] when on right
            color 1 represented by bits [1] when on left
            color 1 represented by bits [0, 2] when on right
            Horizontal colors (coded using 3 bits):
            color 0 represented by bits [0] when on bottom
            color 0 represented by bits [1, 2] when on top
            color 1 represented by bits [1] when on bottom
            color 1 represented by bits [0, 2] when on top
            ([[1, 2, 9],
              [0, 2, 9],
              [0, 4, 5, 10],
              [1, 3, 5, 10],
              [3, 7, 8, 11],
              [4, 6, 8, 11],
              [6, 12],
              [7, 12]],
             [(0, 0, 0),
              (0, 0, 1),
              (1, 0, 0),
              (1, 0, 1),
              (2, 0, 0),
              (2, 0, 1),
              (3, 0, 0),
              (3, 0, 1)])

        ::

            sage: tiles = [(0,0,0,0)]
            sage: W = WangTileSolver(tiles, 4, 1)
            sage: W.rows_and_information(verbose=True)
            Vertical colors (coded using 2 bits):
            color 0 represented by bits [0] when on left
            color 0 represented by bits [1] when on right
            Horizontal colors (coded using 2 bits):
            color 0 represented by bits [0] when on bottom
            color 0 represented by bits [1] when on top
            ([[1, 6], [0, 3, 7], [2, 5, 8], [4, 9]],
             [(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)])

        """
        from math import log, ceil

        # mapping the vertical colors to complementary binary strings
        vertical_alphabet = sorted(self.vertical_alphabet())
        padtoV = int(ceil(log(len(vertical_alphabet)+1, 2)))+1
        left_color_to_digits = {c:[k for (k,d) in 
                                enumerate(ZZ(i+1).digits(2,padto=padtoV)) if d!=0]
                                for i,c in enumerate(vertical_alphabet)}
        right_color_to_digits = {c:[k for (k,d) in 
                                enumerate(ZZ(2**padtoV-1-i-1).digits(2,padto=padtoV)) if d!=0]
                                for i,c in enumerate(vertical_alphabet)}

        # mapping the horizontal colors to complementary binary strings
        horizontal_alphabet = sorted(self.horizontal_alphabet())
        padtoH = int(ceil(log(len(horizontal_alphabet)+1, 2)))+1
        bottom_color_to_digits = {c:[k for (k,d) in 
                                enumerate(ZZ(i+1).digits(2,padto=padtoH)) if d!=0]
                                for i,c in enumerate(horizontal_alphabet)}
        top_color_to_digits = {c:[k for (k,d) in 
                                enumerate(ZZ(2**padtoH-1-i-1).digits(2,padto=padtoH)) if d!=0]
                                for i,c in enumerate(horizontal_alphabet)}

        # Note: above, we want to avoid the binary string 0000 which may
        # create empty rows below for the dlx solver which is not good we
        # want each color to be map to a binary string containing both 0's
        # and 1's

        if verbose:
            phrase = "color {} represented by bits {} when on {}"
            print("Vertical colors (coded using {} bits):".format(padtoV))
            for color in left_color_to_digits:
                digits_left = left_color_to_digits[color]
                digits_right = right_color_to_digits[color]
                print(phrase.format(color, digits_left, 'left'))
                print(phrase.format(color, digits_right, 'right'))
            #
            print("Horizontal colors (coded using {} bits):".format(padtoH))
            for color in bottom_color_to_digits:
                digits_bottom = bottom_color_to_digits[color]
                digits_top = top_color_to_digits[color]
                print(phrase.format(color, digits_bottom, 'bottom'))
                print(phrase.format(color, digits_top, 'top'))

        W = self._width
        H = self._height
        dict_of_rows = defaultdict(list)

        # Preassigned colors
        rightPRE,topPRE,leftPRE,bottomPRE = self._preassigned_color

        # matching vertical colors
        for j in range(W):
            for k in range(H):
                position = (k*(W-1)+j)*padtoV
                for i,tile in enumerate(self._tiles):
                    # the tile i at position (j,k)
                    right,top,left,bottom = tile
                    if (j,k) in rightPRE and rightPRE[(j,k)] != right:
                        continue
                    if (j,k) in leftPRE and leftPRE[(j,k)] != left:
                        continue
                    A = left_color_to_digits[left]
                    B = right_color_to_digits[right]
                    row = []
                    if j > 0:
                        row.extend([position+a-padtoV for a in A])
                    if j < W-1:
                        row.extend([position+b for b in B])
                    dict_of_rows[(j,k,i)] = row
        column_shift = H*(W-1)*padtoV

        # matching horizontal colors
        for j in range(W):
            for k in range(H):
                position = (j*(H-1)+k)*padtoH
                for i,tile in enumerate(self._tiles):
                    # the tile i at position (j,k)
                    right,top,left,bottom = tile
                    if (j,k) in bottomPRE and bottomPRE[(j,k)] != bottom:
                        continue
                    if (j,k) in topPRE and topPRE[(j,k)] != top:
                        continue
                    A = bottom_color_to_digits[bottom]
                    B = top_color_to_digits[top]
                    row = []
                    if k > 0:
                        row.extend([column_shift+position-padtoH+a for a in A])
                    if k < H-1:
                        row.extend([column_shift+position+b for b in B])
                    dict_of_rows[(j,k,i)].extend(row)
        column_shift += W*(H-1)*padtoH

        # exactly one tile at each position
        for j in range(W):
            for k in range(H):
                position = j*H + k
                for i,tile in enumerate(self._tiles):
                    # the tile i at position (j,k)
                    dict_of_rows[(j,k,i)].append(column_shift+position)
        column_shift += W*H

        # Preassigned tiles
        col = itertools.count(column_shift)
        for (j,k),i in self._preassigned_tiles.items():
            dict_of_rows[(j,k,i)].append(next(col))

        # Creation of the rows and row information
        dict_of_rows = dict(dict_of_rows)
        sorted_keys = sorted(dict_of_rows)
        rows = [dict_of_rows[key] for key in sorted_keys]

        return rows, sorted_keys

    def dlx_solver(self):
        r"""
        Return the sage DLX solver of that Wang tiling problem.

        OUTPUT:

            DLX Solver

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: dlx = W.dlx_solver()
            sage: dlx
            Dancing links solver for 63 columns and 24 rows
            sage: dlx.number_of_solutions()
            2

        TESTS::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles,2,2)
            sage: dlx = W.dlx_solver()
            sage: list(dlx.solutions_iterator())
            [[1, 7, 4, 10], [6, 0, 9, 3], [8, 2, 5, 11]]
        """
        from sage.combinat.matrices.dancing_links import dlx_solver
        rows,row_info = self.rows_and_information()
        return dlx_solver(rows)

    def number_of_solutions(self, ncpus=8):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(2,4,2,1), (2,2,2,0), (1,1,3,1), (1,2,3,2), (3,1,3,3),
            ....: (0,1,3,1), (0,0,0,1), (3,1,0,2), (0,2,1,2), (1,2,1,4), (3,3,1,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: W.number_of_solutions()
            908

        ::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles,2,2)
            sage: W.number_of_solutions()
            3
        """
        return self.dlx_solver().number_of_solutions(ncpus=ncpus)

    def solutions_iterator(self):
        r"""
        Iterator over all solutions

        .. NOTE::

            This uses the reduction to dancing links.

        OUTPUT:

            iterator of wang tilings

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(2,4,2,1), (2,2,2,0), (1,1,3,1), (1,2,3,2), (3,1,3,3),
            ....: (0,1,3,1), (0,0,0,1), (3,1,0,2), (0,2,1,2), (1,2,1,4), (3,3,1,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: it = W.solutions_iterator()
            sage: next(it)
            A wang tiling of a 3 x 4 rectangle
            sage: next(it)
            A wang tiling of a 3 x 4 rectangle

        ::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles,2,2)
            sage: list(W.solutions_iterator())
            [A wang tiling of a 2 x 2 rectangle,
             A wang tiling of a 2 x 2 rectangle,
             A wang tiling of a 2 x 2 rectangle]

        With preassigned colors and tiles::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
            sage: t = {(0,1):0}
            sage: c = [{},{},{(1,1):0},{}]
            sage: W = WangTileSolver(tiles,3,3,preassigned_tiles=t,preassigned_color=c)
            sage: S = list(W.solutions_iterator())
            sage: [s._table for s in S]
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[0, 0, 3], [0, 0, 0], [0, 0, 0]]]

        """
        from sage.combinat.matrices.dancing_links import dlx_solver
        rows,row_info = self.rows_and_information()
        dlx = dlx_solver(rows)
        for solution in dlx.solutions_iterator():
            table = [[None]*self._height for _ in range(self._width)]
            for a in solution:
                j,k,i = row_info[a]
                assert table[j][k] is None, "table[{}][{}](={}) is not None".format(j,k,table[j][k])
                table[j][k] = i
            yield WangTiling(table, self._tiles, color=self._color)

class WangTiling(object):
    r"""
    INPUT:

    - ``table`` -- list of lists
    - ``tiles`` -- list of tiles, a tile is a 4-tuple (right color, top
        color, left color, bottom color)
    - ``color`` -- dict (default: None)

    .. NOTE::

        ``table[x][y]`` refers to the tile at position `(x,y)` using the
        cartesian coordinates. Thus, it is **not** using the matrix-like
        coordinates.

    EXAMPLES::

        sage: from slabbe import WangTiling
        sage: tiles = [(0,3,1,4), (1,4,0,3)]
        sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
        sage: tiling = WangTiling(table, tiles)
        sage: tiling
        A wang tiling of a 3 x 4 rectangle

    Using some blank tiles::

        sage: tiles = [(0,3,1,4), (1,4,0,3)]
        sage: table = [[0, 1, None, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
        sage: tiling = WangTiling(table, tiles)
        sage: tiling
        A wang tiling of a 3 x 4 rectangle
    """
    def __init__(self, table, tiles, color=None):
        r"""
        See class for documentation.

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles, 3, 4)
            sage: W.solve()
            A wang tiling of a 3 x 4 rectangle
        """
        self._tiles = tiles
        self._table = table
        self._color = color

    def __repr__(self):
        return "A wang tiling of a {} x {} rectangle".format(self.width(), self.height())

    def height(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.height()
            4
        """
        return len(self._table[0])

    def width(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.width()
            3
        """
        return len(self._table)

    def table(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.table()
            [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
        """
        return self._table

    def transpose(self):
        table = [[column[j] for column in self._table]
                            for j in range(self.height())]
        return WangTiling(table, self._tiles, self._color)

    def horizontal_words_list(self, side=3):
        r"""
        Return a list of horizontal words of colors appearing on a given
        side.

        INPUT

        - ``side`` -- integer in [0,1,2,3], 3 is for bottom

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: tiling.horizontal_words_list()
            [[4, 3, 4], [3, 4, 3], [4, 3, 4], [3, 4, 3]]
            sage: tiling.horizontal_words_list(0)
            [[0, 1, 0], [1, 0, 1], [0, 1, 0], [1, 0, 1]]
        """
        rep = []
        for i in range(self.height()):
            row = [self._table[j][i] for j in range(self.width())]
            L = [self._tiles[r][side] for r in row]
            rep.append(L)
        return rep

    def horizontal_words_dict(self, length):
        r"""
        Return a dict of horizontal words (left to right) of given length
        starting at each position (x,y).

        INPUT:

        - ``length`` -- integer

        OUTPUT:

            dict position -> word

        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.horizontal_words_dict(2)
            {(0, 0): (4, 3),
             (0, 1): (3, 4),
             (0, 2): (4, 3),
             (0, 3): (3, 4),
             (0, 4): (4, 3),
             (1, 0): (3, 4),
             (1, 1): (4, 3),
             (1, 2): (3, 4),
             (1, 3): (4, 3),
             (1, 4): (3, 4)}

        """
        d = {}
        for i in range(self.width()+1-length):
            for j in range(self.height()):
                tile_sequence = [self._table[i+k][j] for k in range(length)]
                color_sequence = tuple(self._tiles[r][3] for r in tile_sequence)
                d[(i,j)] = color_sequence
        # we do once more for the top color of the top row
        H = self.height()
        for i in range(self.width()+1-length):
            tile_sequence = [self._table[i+k][H-1] for k in range(length)]
            color_sequence = tuple(self._tiles[r][1] for r in tile_sequence)
            d[(i,H)] = color_sequence
        return d

    def vertical_words_dict(self, length):
        r"""
        Return a dict of vertical words (bottom to top) of given length
        starting at each position (x,y).

        INPUT:

        - ``length`` -- integer

        OUTPUT:

            dict position -> word

        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.vertical_words_dict(2)
            {(0, 0): (1, 0),
             (0, 1): (0, 1),
             (0, 2): (1, 0),
             (1, 0): (0, 1),
             (1, 1): (1, 0),
             (1, 2): (0, 1),
             (2, 0): (1, 0),
             (2, 1): (0, 1),
             (2, 2): (1, 0),
             (3, 0): (0, 1),
             (3, 1): (1, 0),
             (3, 2): (0, 1)}
            sage: tiling.vertical_words_dict(3)
            {(0, 0): (1, 0, 1),
             (0, 1): (0, 1, 0),
             (1, 0): (0, 1, 0),
             (1, 1): (1, 0, 1),
             (2, 0): (1, 0, 1),
             (2, 1): (0, 1, 0),
             (3, 0): (0, 1, 0),
             (3, 1): (1, 0, 1)}

        """
        d = {}
        for i in range(self.width()):
            for j in range(self.height()+1-length):
                tile_sequence = [self._table[i][j+k] for k in range(length)]
                color_sequence = tuple(self._tiles[r][2] for r in tile_sequence)
                d[(i,j)] = color_sequence
        # we do once more for the right color of the right-most column
        W = self.width()
        for j in range(self.height()+1-length):
            tile_sequence = [self._table[W-1][j+k] for k in range(length)]
            color_sequence = tuple(self._tiles[r][0] for r in tile_sequence)
            d[(W,j)] = color_sequence
        return d

    def plot_points_on_torus(self, M, pointsize=5, color_dict=None):
        r"""
        Plot points modulo some values in x and y.

        INPUT

        - ``M`` -- M is the matrix projection to `\mathbb{R}^2/\mathbb{Z}^2`
        - ``pointsize`` -- positive real number (default:``5``)
        - ``color_dict`` -- dict, tile index -> color or None
          (default:``None``)

        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: z = polygen(QQ, 'z')
            sage: K.<phi> = NumberField(z^2-z-1, 'phi', embedding=AA(golden_ratio))
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: T = WangTiling(table, tiles)
            sage: M = matrix(2, [phi, 0, 0, 0.01])
            sage: G = T.plot_points_on_torus(M)
        """
        number_of_tiles = len(self._tiles)
        alphabet = range(number_of_tiles)

        # Create the color dict
        if color_dict is None:
            from sage.plot.colors import hue
            from random import shuffle
            color = [hue(i/float(number_of_tiles)) for i in range(number_of_tiles)]
            shuffle(color)
            color_dict = dict(zip(alphabet, color))

        from sage.functions.other import floor
        def frac(x):
            return x-floor(x)

        # compute the points
        from sage.modules.free_module_element import vector
        PTS = {a:[] for a in alphabet}
        for i,column in enumerate(self._table):
            for j,a in enumerate(column):
                x,y = M*vector((i,j))
                PTS[a].append((frac(x),frac(y)))

        from sage.plot.graphics import Graphics
        from sage.plot.point import points
        from sage.plot.text import text
        from sage.misc.latex import latex
        from random import randrange

        G = Graphics()
        title = r"Drawing {} on $R^2/Z^2$ with ${}$".format(self, latex(M))
        G += text(title, (.5,1.), axis_coords=True)

        # draw the points
        for a in alphabet:
            PTSa = PTS[a]
            if not PTSa:
                continue
            color = color_dict[a]
            G += points(PTSa, color=color, size=pointsize)
            # add a label on one of them
            px,py = PTSa[randrange(len(PTSa))]
            G += text(str(a), (px,py), color=color)

        # color of points in legend
        for d,a in enumerate(alphabet):
            color = color_dict[a]
            G += points([(-.1,0)], color=color, size=pointsize, legend_label=str(d))

        return G

    def number_of_occurences(self, pattern, avoid_border=0):
        r"""
        Return the number of occurences of the given pattern in the tiling.

        INPUT

        - ``pattern`` -- dict
        - ``avoid_border`` -- integer (default: 0), the size of the border
          to avoid during the computation

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: tiling.number_of_occurences({(0,0):0})
            6
            sage: tiling.number_of_occurences({(0,0):1})
            6
            sage: tiling.number_of_occurences({(0,0):1, (1,0):1})
            0
            sage: tiling.number_of_occurences({(0,0):1, (1,0):1, (0,1):1})
            0
            sage: tiling.number_of_occurences({(0,0):1, (1,0):0, (0,1):0})
            3

        The pattern is translation invariant::

            sage: tiling.number_of_occurences({(0,-1):1})
            6
            sage: tiling.number_of_occurences({(-1,-1):1})
            6
            sage: tiling.number_of_occurences({(-100,-100):1})
            6

        The x coordinates of the pattern corresponds to the x coordinates
        when you plot it::

            sage: tiles = [(0,3,0,4), (1,4,1,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: tiling.number_of_occurences({(0,0):1})
            6
            sage: tiling.number_of_occurences({(0,0):1, (1,0):1})
            4
            sage: tiling.number_of_occurences({(0,0):1, (0,1):1})
            0
            sage: tiling.tikz().pdf(view=False)   # not tested

        When avoiding the border::

            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: tiling.number_of_occurences({(0,0):0}, avoid_border=1)
            1
        """
        xmin = min(x for (x,y) in pattern)
        xmax = max(x for (x,y) in pattern)
        ymin = min(y for (x,y) in pattern)
        ymax = max(y for (x,y) in pattern)
        a = 0
        for i in range(0-xmin+avoid_border, self.width()-xmax-avoid_border):
            for j in range(0-ymin+avoid_border, self.height()-ymax-avoid_border):
                if all(self._table[i+x][j+y] == pattern[(x,y)] for (x,y) in pattern):
                    a += 1
        return a

    def pattern_occurrences(self, shape, avoid_border=0):
        r"""
        Return the number of occurences of every pattern having a given
        shape.

        INPUT

        - ``shape`` -- list, list of coordinates
        - ``avoid_border`` -- integer (default: 0), the size of the border
          to avoid during the computation

        OUTPUT

        a dict where each key is a tuple giving the tiles at each
        coordinate of the shape (in the same order) and values are integers

        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: table = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
            sage: tiles = [(0, 0, 0, 0), (1, 1, 1, 1), (2, 2, 2, 2)]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.pattern_occurrences([(0,0)])
            Counter({(0,): 12})

        ::

            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.pattern_occurrences([(0,0)])
            Counter({(0,): 6, (1,): 6})
            sage: tiling.pattern_occurrences([(0,0), (1,0), (0,1)])
            Counter({(1, 0, 0): 3, (0, 1, 1): 3})

        When avoiding the border::

            sage: tiling.pattern_occurrences([(0,0)], avoid_border=1)
            Counter({(0,): 1, (1,): 1})
            sage: tiling.pattern_occurrences([(0,0)], avoid_border=2)
            Counter()
        """
        xmin = min(x for (x,y) in shape)
        xmax = max(x for (x,y) in shape)
        ymin = min(y for (x,y) in shape)
        ymax = max(y for (x,y) in shape)
        C = Counter()
        for i in range(0-xmin+avoid_border, self.width()-xmax-avoid_border):
            for j in range(0-ymin+avoid_border, self.height()-ymax-avoid_border):
                pattern = tuple(self._table[i+x][j+y] for (x,y) in shape)
                C[pattern] += 1
        return C

    def tile_frequency(self, avoid_border=1):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.tile_frequency()
            {(0,): 1/2,
             (1,): 1/2}
        """
        from sage.rings.rational_field import QQ
        C = self.pattern_occurrences([(0,0)], avoid_border=avoid_border)
        s = sum(C.values())
        return {k:QQ((v,s)) for k,v in C.items()}

    @rename_keyword(fontsize='font')
    def tikz(self, color=None, font=r'\normalsize', rotate=None,
            label_shift=.2, scale=1, transformation_matrix=None):
        r"""
        Return a tikzpicture showing one solution.

        INPUT:

        - ``color`` -- None or dict from tile values -> tikz colors
        - ``font`` -- string (default: ``r'\normalsize'``
        - ``rotate`` -- list or ``None`` (default:``None``) list of four angles
          in degrees like ``(0,0,0,0)``, the rotation angle to apply to each
          label of Wang tiles. If ``None``, it performs a 90 degres rotation
          for left and right labels taking more than one character.
        - ``label_shift`` -- number (default: ``.2``) translation distance
          of the label from the edge
        - ``scale`` -- number (default: ``1``), tikzpicture scale
        - ``transformation_matrix`` -- matrix (default: ``None``), a matrix
          to apply to the coordinate before drawing, it can be in
          ``SL(2,ZZ)`` or not.

        .. TODO::

            - Fix the top and right lines when using the transformation
              matrix.
            - Fix the double drawn edges

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: t = tiling.tikz()
            sage: t
            \documentclass[tikz]{standalone}
            \usepackage{amsmath}
            \begin{document}
            \begin{tikzpicture}[scale=1]
            \tikzstyle{every node}=[font=\normalsize]
            % tile at position (x,y)=(0, 0)
            \draw (0, 0) -- (1, 0);
            \draw (0, 0) -- (0, 1);
            ...
            ... 101 lines not printed (3304 characters in total) ...
            ...
            \node[rotate=0] at (2.8, 3.5) {0};
            \node[rotate=0] at (2.5, 3.8) {0};
            \node[rotate=0] at (2.2, 3.5) {0};
            \node[rotate=0] at (2.5, 3.2) {0};
            \end{tikzpicture}
            \end{document}

        With colors::

            sage: tiles = [(0,2,1,3), (1,3,0,2)]
            sage: color = {0:'white',1:'red',2:'blue',3:'green'}
            sage: W = WangTileSolver(tiles,3,4,color=color)
            sage: tiling = W.solve()
            sage: t = tiling.tikz()

        With colors, alternatively::

            sage: tiles = [(0,2,1,3), (1,3,0,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve('GLPK')
            sage: color = {0:'white',1:'red',2:'blue',3:'green'}
            sage: t = tiling.tikz(color=color)

        Using some blank tiles::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,2), (1,2,0,3)]
            sage: table = [[0, 1, None, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: color = {0:'white',1:'red',2:'blue',3:'green'}
            sage: tiling = WangTiling(table, tiles, color)
            sage: t = tiling.tikz()

        Testing the options::

            sage: tiles = [(0,3,1,2), (1,2,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: color = {0:'white',1:'red',2:'blue',3:'green'}
            sage: t = WangTiling(table, tiles, color).tikz(font=r'\Huge')
            sage: t = WangTiling(table, tiles, color).tikz(rotate=(0,90,0,0))
            sage: t = WangTiling(table, tiles, color).tikz(label_shift=.05)
            sage: t = WangTiling(table, tiles, color).tikz(scale=4)
            sage: m = matrix(2,[1,1,0,1])
            sage: t = WangTiling(table, tiles, color).tikz(transformation_matrix=m)
        """
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector
        if color is None:
            color = self._color
        if transformation_matrix is None:
            transformation_matrix = matrix.identity(2)
        lines = []
        lines.append(r'\begin{{tikzpicture}}[scale={}]'.format(scale))
        lines.append(r'\tikzstyle{{every node}}=[font={}]'.format(font))
        W = self.width()
        H = self.height()
        ## missing lines on the top
        #for j in range(W):
        #    lines.append(r'\draw {} -- {};'.format((j,H), (j+1,H)))
        ## missing lines on the right
        #for k in range(H):
        #    lines.append(r'\draw {} -- {};'.format((W,k), (W,k+1)))
        # the tiles with borders left and below
        for j in range(W):
            for k in range(H):
                i = self._table[j][k]
                if i is None:
                    # this is a blank tile
                    continue
                position = transformation_matrix*vector((j,k))
                tile = self._tiles[i]
                more_lines = tile_to_tikz(tile, position, color=color,
                        size=1, rotate=rotate, label_shift=label_shift,
                        top_right_edges=True)
                lines.extend(more_lines)
        lines.append(r'\end{tikzpicture}')
        from slabbe import TikzPicture
        return TikzPicture('\n'.join(lines))

