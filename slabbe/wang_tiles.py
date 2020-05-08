# -*- coding: utf-8 -*-
r"""
Wang tile solver

We solve the problem of tiling a rectangle by Wang tiles by reducing it to
other well-known problems like linear problem, exact cover problem and SAT.

We thus use MILP solvers like Coin or Gurobi, Sat solvers like
cryptominisat, picosat or glucose and dancing links solver which is already
in Sage.

Coin can be installed with::

    sage -i cbc sagelib

Cryptominisat can be installed with::

    sage -i cryptominisat sagelib

Glucose can be installed with::

    sage -i glucose

EXAMPLES::

    sage: from slabbe import WangTileSolver
    sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: tiling = W.solve()
    sage: _ = tiling.tikz().pdf(view=False)

Using different kind of solvers::

    sage: tiling = W.solve(solver='GLPK')
    sage: tiling = W.solve(solver='dancing_links')
    sage: tiling = W.solve(solver='Gurobi')         # optional Gurobi
    sage: tiling = W.solve(solver='cryptominisat')  # optional cryptominisat

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
from sage.graphs.graph import Graph

from sage.misc.decorators import rename_keyword
from sage.misc.superseded import deprecated_function_alias


def tile_to_tikz(tile, position, color=None, id=None, id_color='',
        id_format='{}', sizex=1, sizey=1, rotate=None, label=True,
        label_shift=.2, label_color='black', right_edges=True,
        top_edges=True, left_edges=True, bottom_edges=True, draw_H=None,
        draw_V=None):
    r"""

    INPUT:

    - ``tile`` -- tuple of length 4
    - ``position`` -- tuple of two numbers
    - ``color`` -- dict (default: ``None``) from tile values -> tikz colors
    - ``id`` -- id (default: ``None``) of the tile to be printed in the center
    - ``id_color`` -- string (default: ``''``) 
    - ``id_format`` -- string (default: ``r'{}'``) to be called with
      ``id_format.format(key)``
    - ``sizex`` -- number (default: ``1``), horizontal size of the tile
    - ``sizey`` -- number (default: ``1``), vertical size of the tile
    - ``rotate`` -- list or ``None`` (default:``None``) list of four angles
      in degrees like ``(0,0,0,0)``, the rotation angle to apply to each
      label of Wang tiles. If ``None``, it performs a 90 degres rotation
      for left and right labels taking more than one character.
    - ``label`` -- boolean (default: ``True``) or a tuple of four boolean
    - ``label_shift`` -- number (default: ``.2``) translation distance of the
      label from the edge
    - ``label_color`` -- string (default: ``'black'``)
    - ``right_edges`` -- bool (default: ``True``) 
    - ``top_edges`` -- bool (default: ``True``) 
    - ``left_edges`` -- bool (default: ``True``) 
    - ``bottom_edges`` -- bool (default: ``True``) 
    - ``draw_H`` -- dict (default: ``None``) from tile values -> tikz
      draw commands. If ``None`` the values of the dict get replaced by
      straight lines, more precisely by ``r'\draw {{}} -- ++ (1,0);'``.
      Dict values must be strings ``s`` such that ``s.format((x,y))``
      works.
    - ``draw_V`` -- dict (default: ``None``) from tile values -> tikz
      draw commands. If ``None`` the values of the dict get replaced by
      straight lines, more precisely by ``r'\draw {{}} -- ++ (0,1);'``.
      Dict values must be strings ``s`` such that ``s.format((x,y))``
      works.

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
         '\\draw (11, 100) -- ++ (0,1);',
         '\\draw (10, 101) -- ++ (1,0);',
         '\\draw (10, 100) -- ++ (0,1);',
         '\\draw (10, 100) -- ++ (1,0);',
         '\\node[rotate=0,black] at (10.8, 100.5) {1};',
         '\\node[rotate=0,black] at (10.5, 100.8) {2};',
         '\\node[rotate=0,black] at (10.2, 100.5) {3};',
         '\\node[rotate=0,black] at (10.5, 100.2) {4};']
        sage: tile_to_tikz((1,2,3,4), (10,100), color=None)
        ['% tile at position (x,y)=(10, 100)',
         '\\draw (11, 100) -- ++ (0,1);',
         '\\draw (10, 101) -- ++ (1,0);',
         '\\draw (10, 100) -- ++ (0,1);',
         '\\draw (10, 100) -- ++ (1,0);',
         '\\node[rotate=0,black] at (10.8, 100.5) {1};',
         '\\node[rotate=0,black] at (10.5, 100.8) {2};',
         '\\node[rotate=0,black] at (10.2, 100.5) {3};',
         '\\node[rotate=0,black] at (10.5, 100.2) {4};']
        sage: tile_to_tikz((1,2,3,4), (10,100), color=None, rotate=(0,90,0,0))
        ['% tile at position (x,y)=(10, 100)',
         '\\draw (11, 100) -- ++ (0,1);',
         '\\draw (10, 101) -- ++ (1,0);',
         '\\draw (10, 100) -- ++ (0,1);',
         '\\draw (10, 100) -- ++ (1,0);',
         '\\node[rotate=0,black] at (10.8, 100.5) {1};',
         '\\node[rotate=90,black] at (10.5, 100.8) {2};',
         '\\node[rotate=0,black] at (10.2, 100.5) {3};',
         '\\node[rotate=0,black] at (10.5, 100.2) {4};']
        sage: tile_to_tikz((1,2,3,4), (10,100), color=None, label_shift=.1)
        ['% tile at position (x,y)=(10, 100)',
         '\\draw (11, 100) -- ++ (0,1);',
         '\\draw (10, 101) -- ++ (1,0);',
         '\\draw (10, 100) -- ++ (0,1);',
         '\\draw (10, 100) -- ++ (1,0);',
         '\\node[rotate=0,black] at (10.9000000000000, 100.5) {1};',
         '\\node[rotate=0,black] at (10.5, 100.900000000000) {2};',
         '\\node[rotate=0,black] at (10.1000000000000, 100.5) {3};',
         '\\node[rotate=0,black] at (10.5, 100.100000000000) {4};']

    ::

        sage: tile_to_tikz((10,20,30,40), (10,100), color=None)
        ['% tile at position (x,y)=(10, 100)',
         '\\draw (11, 100) -- ++ (0,1);',
         '\\draw (10, 101) -- ++ (1,0);',
         '\\draw (10, 100) -- ++ (0,1);',
         '\\draw (10, 100) -- ++ (1,0);',
         '\\node[rotate=90,black] at (10.8, 100.5) {10};',
         '\\node[rotate=0,black] at (10.5, 100.8) {20};',
         '\\node[rotate=90,black] at (10.2, 100.5) {30};',
         '\\node[rotate=0,black] at (10.5, 100.2) {40};']
    """
    if rotate is None:
        rotate = []
        for i,a in enumerate(tile):
            if (i == 0 or i == 2) and len(str(a)) > 1:
                rotate.append(90)
            else:
                rotate.append(0)
    sx = sizex      # because it is shorter to write below
    sy = sizey      # because it is shorter to write below
    t = label_shift # because it is shorter to write below
    x,y = position
    lines = []
    #lines.append(r'\begin{tikzpicture}')
    lines.append('% tile at position (x,y)={}'.format((x,y)))
    if color:
        triangle = r'\fill[{}] {} -- {} -- {};'
        c = (x+.5*sx,y+.5*sy)
        lines.append(triangle.format(color[tile[0]],(x+sx,y),c,(x+sx,y+sy)))
        lines.append(triangle.format(color[tile[1]],(x,y+sy),c,(x+sx,y+sy)))
        lines.append(triangle.format(color[tile[2]],(x,y),c,(x,y+sy)))
        lines.append(triangle.format(color[tile[3]],(x,y),c,(x+sx,y)))

    if id is not None:
        c = (x+.5*sx,y+.5*sy)
        id = id_format.format(id)
        lines.append(r'\node[{}] at {} {{{}}};'.format(id_color, c, id))

    if draw_H is None:
        draw_H = {tile[1]:r'\draw {{}} -- ++ ({},0);'.format(sx),
                  tile[3]:r'\draw {{}} -- ++ ({},0);'.format(sx)}
    if draw_V is None:
        draw_V = {tile[0]:r'\draw {{}} -- ++ (0,{});'.format(sy),
                  tile[2]:r'\draw {{}} -- ++ (0,{});'.format(sy)}

    if right_edges:
        lines.append(draw_V[tile[0]].format((x+sx,y)))
    if top_edges:
        lines.append(draw_H[tile[1]].format((x,y+sy)))
    if left_edges:
        lines.append(draw_V[tile[2]].format((x,y)))
    if bottom_edges:
        lines.append(draw_H[tile[3]].format((x,y)))

    node_str = r'\node[rotate={},{}] at {} {{{}}};'
    if isinstance(label, bool):
        if label:
            label = (True, True, True, True)
        else:
            label = (False, False, False, False)
    if label[0]:
        lines.append(node_str.format(rotate[0],label_color,(x+sx-t,  y+.5*sy),  tile[0]))
    if label[1]:
        lines.append(node_str.format(rotate[1],label_color,(x+.5*sx, y+sy-t), tile[1]))
    if label[2]:
        lines.append(node_str.format(rotate[2],label_color,(x+t,     y+.5*sy),  tile[2]))
    if label[3]:
        lines.append(node_str.format(rotate[3],label_color,(x+.5*sx, y+t),   tile[3]))

    #lines.append(r'\end{tikzpicture}')
    return lines
    #return TikzPicture('\n'.join(lines))


def fusion(tile0, tile1, direction, function=str.__add__, initial=''):
    r"""
    Return the fusion of wang tile sets in the given direction.

    We keep only the strongly connected components.

    INPUT:

    - ``tile0`` -- 4-uple
    - ``tile1`` -- 4-uple
    - ``direction`` -- integer (1 or 2)
    - ``function`` -- function (default:``str.__add__``), monoid
        operation
    - ``initial`` -- object (default:``''``), monoid neutral

    EXAMPLES::

        sage: from slabbe.wang_tiles import fusion
        sage: t0 = 'abcd'
        sage: t1 = 'xyaz'
        sage: fusion(t0,t1,1)
        ('x', 'by', 'c', 'dz')

    ::

        sage: t0 = 'abcd'
        sage: t1 = 'xyzb'
        sage: fusion(t0,t1,2)
        ('ax', 'y', 'cz', 'd')

    TESTS::

        sage: t0 = 'abcd'
        sage: t1 = 'efgh'
        sage: fusion(t0,t1,1)
        Traceback (most recent call last):
        ...
        AssertionError: A must be equal to Y

    """
    A,B,C,D = tile0
    W,X,Y,Z = tile1
    if direction == 1:
        assert A == Y, "A must be equal to Y"
        t = ((W,), (B,X), (C,), (D,Z))
    elif direction == 2:
        assert B == Z, "B must be equal to Z"
        t = ((A,W), (X,), (C,Y), (D,))
    from functools import reduce
    return tuple(reduce(function, a, initial) for a in t)

class WangTileSet(object):
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
        self._tiles = list(tiles)

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

    def __getitem__(self, i):
        r"""
        INPUT:

        - ``i`` -- integer, index

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: T = WangTileSet(tiles)
            sage: T[1]
            (1, 1, 1, 1)
        """
        return self._tiles[i]

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

    def table(self):
        r"""
        Return a table representation of the tile set.

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
            ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
            sage: T = WangTileSet(tiles)
            sage: T.table()
              Id   Right   Top   Left   Bottom
            +----+-------+-----+------+--------+
              0    0       0     0      2
              1    1       0     0      1
              2    2       1     0      0
              3    0       0     1      0
              4    1       2     1      1
              5    1       1     2      0
              6    2       0     2      1

        """
        from sage.misc.table import table
        header_row = ['Id', 'Right', 'Top', 'Left', 'Bottom']
        rows = [(i, a, b, c, d) for i,(a,b,c,d) in enumerate(self.tiles())]
        return table(rows, header_row=header_row)

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

    def to_transducer_graph(self, label_function=tuple,
            merge_multiedges=True):
        r"""
        Return the graph of the transducer.

        Labels are cleaned. Label of multiedges are merged with commas.

        INPUT:

        - ``label_function`` -- function (default:``tuple``), a function to
          apply to each list of labels when merging multiedges into one
        - ``merge_multiedges`` -- boolean (default:``True``)

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
        def edge_labels(t):
            assert len(t.word_in) == 1, "we assume word_in is of length 1"
            assert len(t.word_out) == 1, "we assume word_out is of length 1"
            return "{}|{}".format(t.word_in[0], t.word_out[0])
        G = self.to_transducer().graph(edge_labels)
        if merge_multiedges:
            from slabbe.graph import merge_multiedges
            return merge_multiedges(G, label_function)
        else:
            return  G

    def system_of_density_equations(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,2), (1,0,0,1), (2,1,0,0), (0,0,1,0),
            ....:          (1,2,1,1), (1,1,2,0), (2,0,2,1)]
            sage: T = WangTileSet(tiles)
            sage: M = T.system_of_density_equations()
            sage: M
            [ 1  1  1  1  1  1  1  1]
            [ 1  1 -1  0  0 -1  1  0]
            [ 0  1  1 -1  0  0  0  0]
            [ 0  0 -1  0  0  1  0  0]
            [ 0 -1  1  0 -1  1 -1  0]
            [ 0 -1  0  1  0 -1  0  0]
            [-1  0  0  0  1  0  0  0]
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
        rows = list(vertical.values())
        rows.extend(horizontal.values())
        rows.append(M([1 for _ in range(len(self)+1)]))
        rows.sort(reverse=True)
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
        eqns = list(vertical.values())
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
             font=r'\normalsize', rotate=None, label=True, id=True,
             id_format='{}', id_color='', label_shift=.2,
             label_color='black', right_edges=True, top_edges=True,
             left_edges=True, bottom_edges=True, draw_H=None, draw_V=None):
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
        - ``label`` -- boolean (default: ``True``), presence of the colors
        - ``id`` -- boolean (default: ``True``), presence of the tile id
        - ``id_color`` -- string (default: ``''``) 
        - ``id_format`` -- string (default: ``r'{}'``) to be called with
          ``id_format.format(key)``
        - ``label_shift`` -- number (default: ``.2``) translation distance of the
          label from the edge
        - ``label_color`` -- string (default: ``'black'``)
        - ``right_edges`` -- bool (default: ``True``) 
        - ``top_edges`` -- bool (default: ``True``) 
        - ``left_edges`` -- bool (default: ``True``) 
        - ``bottom_edges`` -- bool (default: ``True``) 
        - ``draw_H`` -- dict (default: ``None``) from tile values -> tikz
          draw commands. If ``None`` the values of the dict get replaced by
          straight lines, more precisely by ``r'\draw {{}} -- ++ (1,0);'``.
          Dict values must be strings ``s`` such that ``s.format((x,y))``
          works.
        - ``draw_V`` -- dict (default: ``None``) from tile values -> tikz
          draw commands. If ``None`` the values of the dict get replaced by
          straight lines, more precisely by ``r'\draw {{}} -- ++ (0,1);'``.
          Dict values must be strings ``s`` such that ``s.format((x,y))``
          works.

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
            this_id = i if id else None
            new_lines = tile_to_tikz(tile, position, color=color,
                    id=this_id, id_color=id_color, id_format=id_format,
                    sizex=size, sizey=size, rotate=rotate, label=label,
                    label_shift=label_shift, label_color=label_color,
                    right_edges=right_edges, top_edges=top_edges,
                    left_edges=left_edges, bottom_edges=bottom_edges,
                    draw_H=draw_H, draw_V=draw_V)
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
            label_color='black', rotate=None, label=True, label_shift=.2,
            id=True, id_color='', id_format='{}', draw_H=None, draw_V=None):
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
        - ``label`` -- boolean (default: ``True``), presence of the colors
        - ``label_shift`` -- number (default: ``.2``) translation distance
          of the label from the edge
        - ``label_color`` -- string (default: ``'black'``)
        - ``id`` -- boolean (default: ``True``), presence of the tile id
        - ``id_color`` -- string (default: ``''``) 
        - ``id_format`` -- string (default: ``r'{}'``) to be called with
          ``id_format.format(key)``
        - ``draw_H`` -- dict (default: ``None``) from tile values -> tikz
          draw commands. If ``None`` the values of the dict get replaced by
          straight lines, more precisely by ``r'\draw {{}} -- ++ (1,0);'``.
          Dict values must be strings ``s`` such that ``s.format((x,y))``
          works.
        - ``draw_V`` -- dict (default: ``None``) from tile values -> tikz
          draw commands. If ``None`` the values of the dict get replaced by
          straight lines, more precisely by ``r'\draw {{}} -- ++ (0,1);'``.
          Dict values must be strings ``s`` such that ``s.format((x,y))``
          works.

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
            this_id = i if id else None
            new_lines = tile_to_tikz(tile, position=(0,0), color=color,
                    label=label, label_color=label_color, id=this_id,
                    id_color=id_color, id_format=id_format, sizex=size,
                    sizey=size, rotate=rotate, draw_H=draw_H,
                    draw_V=draw_V, label_shift=label_shift)
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
            label_shift=.2, ncolumns=4, tabular='tabular', align='l'):
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

    def clean_sources_and_sinks(self):
        r"""
        TODO: do it for the dual?

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (3,2,4,8), (0,5,0,7)]
            sage: T = WangTileSet(tiles)
            sage: T.clean_sources_and_sinks().tiles()
            [(0, 0, 0, 0), (0, 5, 0, 7), (1, 1, 1, 1)]
            sage: T.dual().clean_sources_and_sinks().tiles()
            [(0, 0, 0, 0), (1, 1, 1, 1)]

        """
        from slabbe.graph import clean_sources_and_sinks
        def edge_labels(t):
            return (t.word_in[0], t.word_out[0])
        G = self.to_transducer().graph(edge_labels)
        G = clean_sources_and_sinks(G)
        tiles = []
        for (u,v,(word_in, word_out)) in G.edges():
            left = u
            right = v
            bottom = word_in
            top = word_out
            tiles.append((right, top, left, bottom))
        return WangTileSet(tiles)

    def fusion(self, other, direction, function=str.__add__, initial='',
            clean_graph=True):
        r"""
        Return the fusion of wang tile sets in the given direction.

        TODO: check if and when to do the clean
        
        INPUT:

        - ``other`` -- WangTileSet
        - ``direction`` -- integer (1 or 2)
        - ``function`` -- function (default:``str.__add__``), monoid
          operation
        - ``initial`` -- object (default:``''``), monoid neutral
        - ``clean_graph`` -- boolean (default: ``False``), clean the graph
          by recursively removing sources and sinks transitions (or tiles).

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

        TESTS::

            sage: tiles = [('02', '2', '02', '2'), ('32', '2', '02', '2')]
            sage: T = WangTileSet(tiles)
            sage: T.fusion(T, 1)
            Wang tile set of cardinality 2
            sage: T.fusion(T, 2)
            Wang tile set of cardinality 1
        """
        if not isinstance(other, WangTileSet):
            raise TypeError('other(={}) must be a'
                    ' WangTileSet'.format(other))
        if direction == 1:
            return self.dual().fusion(other.dual(), direction=2,
                    function=function, initial=initial,
                    clean_graph=clean_graph).dual()
        if not direction == 2:
            raise ValueError('direction(={}) must be 1 or 2'.format(direction))

        T = self.to_transducer()
        U = other.to_transducer()
        TU = T.composition(U)

        #edge_labels = lambda t:"{}|{}".format(t.word_in, t.word_out)
        edge_labels = lambda t:(t.word_in, t.word_out)
        G = TU.graph(edge_labels)
        if clean_graph:
            # BUG? maybe this should be done after the reduce?
            from slabbe.graph import clean_sources_and_sinks
            G = clean_sources_and_sinks(G)

        from functools import reduce
        tiles = []
        for (u,v,(word_in, word_out)) in G.edges():
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

        When incompatible preassigned tiles::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: T = WangTileSet(tiles)
            sage: W = T.solver(3,3, preassigned_tiles={(0,0):0,(0,1):1})
            sage: W.has_solution()
            False

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

    def tiles_allowing_surrounding(self, radius, solver=None, ncpus=None, verbose=False):
        r"""
        Return the subset of tiles allowing a surrounding of given radius.

        INPUT:

        - ``radius`` - integer
        - ``solver`` - string or None
        - ``ncpus`` - integer
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
            [[0, 0, ...], [0, 0, 0], [0, 0, 0]]
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
            if s.has_solution(solver=solver, ncpus=ncpus):
                if verbose:
                    print("Solution found for tile {}:\n{}".format(i,
                                s.solve(solver)._table))
                L.append(t)
        return WangTileSet(L)

    def _tilings_with_surrounding(self, width, height, radius=1, solver=None,
            verbose=False):
        r"""
        Return the set of valid tiling of a rectangle of given width and
        height allowing a surrounding of itself of given radius.

        This algorithm is using fusion of Wang tiles horizontally and
        vertically and then using a solver for the new set of tiles and
        then tiling the obtained boundary solutions.

        See also the other method :meth:`tilings_with_surrounding`.

        INPUT:

        - ``width`` - integer
        - ``height`` - integer
        - ``radius`` - integer
        - ``solver`` - string or None
        - ``verbose`` - boolean

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
            sage: tiles = [[str(a) for a in tile] for tile in tiles]
            sage: T = WangTileSet(tiles)
            sage: S = T._tilings_with_surrounding(2,2)
            sage: S
            [A wang tiling of a 2 x 2 rectangle,
             A wang tiling of a 2 x 2 rectangle,
             A wang tiling of a 2 x 2 rectangle]
            sage: [a.table() for a in S]
            [[[0, 0], [0, 0]], [[1, 1], [1, 1]], [[2, 2], [2, 2]]]

        ::

            sage: S = T._tilings_with_surrounding(3,3)
            sage: S
            [A wang tiling of a 3 x 3 rectangle,
             A wang tiling of a 3 x 3 rectangle,
             A wang tiling of a 3 x 3 rectangle]
            sage: [a.table() for a in S]
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
             [[2, 2, 2], [2, 2, 2], [2, 2, 2]]]

        TESTS::

            sage: tiles = ['ABCD', 'EFGH', 'AXCY', 'ABAB', 'EBEB']
            sage: T = WangTileSet(tiles)
            sage: solutions = T._tilings_with_surrounding(1,2)
            sage: [t.table() for t in solutions]
            [[[3, 3]], [[3, 4]], [[4, 3]], [[4, 4]]]

        ::

            sage: tiles = [('02', '4', '02', '4'), ('32', '4', '02', '4')]
            sage: T = WangTileSet(tiles)
            sage: [t.table() for t in T._tilings_with_surrounding(1,2)]
            [[[0, 0]]]
            sage: [t.table() for t in T._tilings_with_surrounding(2,1)]
            [[[0], [0]]]

        """
        tiles_tuple = [tuple((a,) for a in t) for t in self.tiles()]
        T = base = WangTileSet(tiles_tuple)
        for _ in range(width-1):
            T = T.fusion(base, 1, function=tuple.__add__, initial=tuple(),
                    clean_graph=False)
        if verbose:
            print("After fusion in the direction e1:\n", T.table())
        base = T
        for _ in range(height-1):
            T = T.fusion(base, 2, function=tuple.__add__, initial=tuple(),
                    clean_graph=False)
        if verbose:
            print("After fusion in the direction e2:\n", T.table())
        T = T.tiles_allowing_surrounding(radius, solver=solver, verbose=verbose)
        if verbose:
            print("After filtering tiles without surrounding of "
                  "radius {} :\n {}".format(radius, T.table()))
        L = []
        for t in T:
            assert len(t[0]) == height
            assert len(t[1]) == width
            assert len(t[2]) == height
            assert len(t[3]) == width
            right = {(width-1,i):a for (i,a) in enumerate(t[0])}
            top = {(i,height-1):a for (i,a) in enumerate(t[1])}
            left = {(0,i):a for (i,a) in enumerate(t[2])}
            bottom = {(i,0):a for (i,a) in enumerate(t[3])}
            preassigned_color=[right,top,left,bottom]
            W = self.solver(width, height,
                    preassigned_color=preassigned_color)
            L.extend(W.solutions_iterator())
        return L

    @cached_method
    def tilings_with_surrounding(self, width, height, radius=1, solver=None,
            ncpus=1, use_previous_smaller_results=True, verbose=False):
        r"""
        Return the set of valid tiling of a rectangle of given width and
        height allowing a surrounding of itself of given radius.

        INPUT:

        - ``width`` - integer
        - ``height`` - integer
        - ``radius`` - integer or 2-tuple (default: ``1``), if 2-tuple is
          given, then it is interpreted as ``(xradius, yradius)``
        - ``solver`` - string or None (default: ``None``)
        - ``ncpus`` -- integer (default: ``1``), maximal number of
          subprocesses to use at the same time, used only if ``solver`` is
          ``'dancing_links'``.
        - ``use_previous_smaller_results`` - bool (default: ``True``)
        - ``verbose`` - bool (default: ``False``)

        .. NOTE::

            The ``solver='dancing_links'`` is fast for this question (I think)

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
            sage: tiles = [[str(a) for a in tile] for tile in tiles]
            sage: T = WangTileSet(tiles)
            sage: S = T.tilings_with_surrounding(2,2)
            sage: sorted(tiling.table() for tiling in S)
            [[[0, 0], [0, 0]], [[1, 1], [1, 1]], [[2, 2], [2, 2]]]

        ::

            sage: S = T.tilings_with_surrounding(3,3)
            sage: sorted(tiling.table() for tiling in S)
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
             [[2, 2, 2], [2, 2, 2], [2, 2, 2]]]

        TESTS::

            sage: tiles = ['ABCD', 'EFGH', 'AXCY', 'ABAB', 'EBEB']
            sage: T = WangTileSet(tiles)
            sage: solutions = T.tilings_with_surrounding(1,2)
            sage: sorted(tiling.table() for tiling in solutions)   # random
            [[[3, 3]], [[4, 4]], [[3, 4]], [[4, 3]]]

        ::

            sage: tiles = [('02', '4', '02', '4'), ('32', '4', '02', '4')]
            sage: T = WangTileSet(tiles)
            sage: [tiling.table() for tiling in T.tilings_with_surrounding(1,2)]
            [[[0, 0]]]
            sage: [tiling.table() for tiling in T.tilings_with_surrounding(2,1)]
            [[[0], [0]]]

        """
        if isinstance(radius, tuple):
            xradius,yradius = radius
        else:
            xradius = yradius = radius

        # Check if we cached the result of previous calls with smaller radius
        previous_calls = [args for (args, kwds) in self.tilings_with_surrounding.cache
                               if args[0] == width and args[1] == height]
        smaller_previous_calls = {}
        for args in previous_calls:
            radius_ = args[2]
            if isinstance(radius_, tuple):
                xradius_,yradius_ = radius_
            else:
                xradius_ = yradius_ = radius_
            if xradius_ <= xradius and yradius_ <= yradius:
                # the following may erase another previous call, we don't care
                smaller_previous_calls[xradius_+yradius_] = args 
        if use_previous_smaller_results and smaller_previous_calls:
            args = smaller_previous_calls[max(smaller_previous_calls)]
            previous_solutions = self.tilings_with_surrounding.cache[(args,())]
        else:
            W = WangTileSolver(self._tiles, width, height)
            previous_solutions = W.solutions_iterator()

        total_width = 2*xradius+width
        total_height = 2*yradius+height

        L = set()
        for tiling in previous_solutions:
            table = tiling._table
            d = {(x+xradius,y+yradius):table[x][y] for x in range(width) 
                                                   for y in range(height)}
            s = self.solver(total_width, total_height, preassigned_tiles=d)
            if s.has_solution(solver=solver, ncpus=ncpus):
                if verbose:
                    print("Solution found for pattern {}:\n {}".format(d, s.solve(solver)._table))
                L.add(tiling)
        return L

    @cached_method
    def dominoes_with_surrounding(self, i=2, radius=1, solver=None,
            ncpus=1, verbose=False):
        r"""

        INPUT:

        - ``i`` - integer (default: ``2``), 1 or 2
        - ``radius`` - integer or 2-tuple (default: ``1``), if 2-tuple is
          given, then it is interpreted as ``(xradius, yradius)``
        - ``solver`` - string or None (default: ``None``)
        - ``ncpus`` -- integer (default: ``1``), maximal number of
          subprocesses to use at the same time, used only if ``solver`` is
          ``'dancing_links'``.
        - ``verbose`` - bool

        .. NOTE::

            The ``solver='dancing_links'`` is fast for this question.

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = ['ABCD', 'EFGH', 'AXCY', 'ABAB', 'EBEB']
            sage: T = WangTileSet(tiles)
            sage: sorted(T.dominoes_with_surrounding(i=1))
            [(3, 3), (4, 4)]
            sage: sorted(T.dominoes_with_surrounding(i=2))
            [(3, 3), (3, 4), (4, 3), (4, 4)]
            sage: sorted(T.dominoes_with_surrounding(i=2, radius=2))
            [(3, 3), (3, 4), (4, 3), (4, 4)]
            sage: sorted(T.dominoes_with_surrounding(i=2, radius=(1,2)))
            [(3, 3), (3, 4), (4, 3), (4, 4)]

        TESTS::

            sage: tiles = [('02', '4', '02', '4'), ('32', '4', '02', '4')]
            sage: T = WangTileSet(tiles)
            sage: sorted(T.dominoes_with_surrounding(1))
            [(0, 0)]
            sage: sorted(T.dominoes_with_surrounding(2))
            [(0, 0)]

        """
        if isinstance(radius, tuple):
            xradius,yradius = radius
        else:
            xradius = yradius = radius

        if i == 2:
            width = 2*xradius+1
            height = 2*yradius+2
            p = (xradius,yradius)
            q = (xradius,yradius+1)
        elif i == 1:
            width = 2*xradius+2
            height = 2*yradius+1
            p = (xradius,yradius)
            q = (xradius+1,yradius)
        else:
            raise ValueError('i={} must be 1 or 2'.format(i))

        L = set()
        for i,j in itertools.product(range(len(self)), repeat=2):
            d = {p:i, q:j}
            s = self.solver(width, height, preassigned_tiles=d)
            if s.has_solution(solver=solver, ncpus=ncpus):
                if verbose:
                    print("Solution found for tiles {} and {}:\n{}".format(i,j,
                                s.solve(solver)._table))
                L.add((i,j))
        return L

    tiling_with_surrounding = deprecated_function_alias(123456,tilings_with_surrounding)
    not_forbidden_tilings = deprecated_function_alias(123456,tilings_with_surrounding)
    not_forbidden_dominoes = deprecated_function_alias(123456,dominoes_with_surrounding)
    def is_pattern_surrounded(self, pattern, radius=1, solver=None, ncpus=None):
        r"""
        Return whether the rectangular pattern allows a surrounding of
        given radius.

        INPUT:

        - ``pattern`` -- list of lists of tile indices
        - ``radius`` - integer or 2-tuple (default: ``1``), if 2-tuple is
          given, then it is interpreted as ``(xradius, yradius)``
        - ``solver`` -- string (default:``None``)
        - ``ncpus`` -- integer (default:``None``)

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: T = WangTileSet(tiles)
            sage: T.is_pattern_surrounded([[0,0]], solver='dancing_links')
            True
            sage: T.is_pattern_surrounded([[0,1]], solver='dancing_links')
            False
            sage: T.is_pattern_surrounded([[0],[1]], solver='dancing_links')
            False
            sage: T.is_pattern_surrounded([[0],[0]], solver='dancing_links')
            True

        """
        if isinstance(radius, tuple):
            xradius,yradius = radius
        else:
            xradius = yradius = radius

        width = len(pattern) + 2 * xradius
        height = len(pattern[0]) + 2 * yradius

        if isinstance(pattern, list):
            preassigned_tiles = {(xradius+a,yradius+b):tile 
                                for (a,col) in enumerate(pattern)
                                for (b,tile) in enumerate(col)}
        elif isinstance(pattern, dict):
            raise NotImplementedError('when pattern is a dict')
        else:
            raise TypeError('pattern(={}) must be a list or a dict'.format(pattern))

        S = self.solver(width=width, height=height,
                             preassigned_tiles=preassigned_tiles)

        return S.has_solution(solver=solver, ncpus=ncpus)

    def is_forbidden_product(self, A, B, i=2, radius=1, solver=None, ncpus=None):
        r"""
        Return whether A \odot^i B is forbidden using a given radius around
        the product and a given solver.

        INPUT:

        - ``A`` -- list of tile indices
        - ``B`` -- list of tile indices
        - ``i`` -- integer, 1 or 2
        - ``radius`` -- integer (default:``1``)
        - ``solver`` -- string (default:``None``)
        - ``ncpus`` -- integer (default:``None``)

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = ['ABCD', 'EFGH', 'AXCY', 'ABAB']
            sage: T = WangTileSet(tiles)
            sage: T.is_forbidden_product([3],[3])
            False
            sage: T.is_forbidden_product([0],[0])
            True
            sage: T.is_forbidden_product([0,1],[0,1])
            True
            sage: T.is_forbidden_product([0,1],[0,1,2])
            True
            sage: T.is_forbidden_product([0,1],[0,1,2,3])
            True
            sage: T.is_forbidden_product([0,1,3],[0,1,2,3])
            False

        """
        if i == 1:
            p = (radius,radius)
            q = (radius+1,radius)
            width = 2*radius + 2
            height = 2*radius + 1
        elif i == 2:
            p = (radius,radius)
            q = (radius,radius+1)
            width = 2*radius + 1
            height = 2*radius + 2
        for ta,tb in itertools.product(A, B):
            s = self.solver(width, height, preassigned_tiles={p:ta, q:tb})
            if s.has_solution(solver=solver, ncpus=ncpus):
                return False
        return True

    def find_markers_with_slope(self, i=2, slope=None, radius=1, solver=None,
            ncpus=1, verbose=False):
        r"""
        Return a list of lists of marker tiles.

        INPUT:

        - ``i`` -- integer (default:``2``), 1 or 2. 
        - ``slope`` -- -1, 0, 1 or Infinity (default:``None``)
        - ``radius`` - integer or 2-tuple (default: ``1``), if 2-tuple is
          given, then it is interpreted as ``(xradius, yradius)``
        - ``solver`` -- string (default:``None``)
        - ``ncpus`` -- integer (default:``1``)
        - ``verbose`` -- boolean (default:``False``)

        OUTPUT:

            list of lists

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = ['ABCD', 'EFGH', 'AXCY', 'ABAB']
            sage: T = WangTileSet(tiles)
            sage: T.find_markers_with_slope(i=1, slope=1)   # known bug
            [{0, 3}, {1}, {2}]
            sage: T.find_markers_with_slope(i=2, slope=1)   # known bug
            [{0, 2, 3}, {1}]

        """
        from sage.rings.infinity import Infinity
        if slope is None:
            if i == 1:
                slope = Infinity
            elif i == 2:
                slope = 0

        if not i in [1, 2]:
            raise ValueError("i(={}) should be 1 or 2".format(i))
        if i == 1 and slope not in [1, -1, Infinity]:
            raise ValueError("slope(={}) should be -1, 1 or +Infinity when i=1".format(slope))
        elif i == 2 and slope not in [1, -1, 0]:
            raise ValueError("slope(={}) should be 0, -1, 1 when i=2".format(slope))

        edges = []
        dominoes = self.dominoes_with_surrounding(i=i, radius=radius,
                                        solver=solver, ncpus=ncpus)
        if verbose:
            print('dominoes=', dominoes)
        for A,B in dominoes:
            rightA, topA, leftA, bottomA = self[A]
            rightB, topB, leftB, bottomB = self[B]
            if i == 1:
                if slope == Infinity:
                    edges.append((bottomB,topB,B))
                elif slope == 1:
                    edges.append((bottomA,topB,B))
                elif slope == -1:
                    edges.append((topA,bottomB,B))
            elif i == 2:
                if slope == 0:
                    edges.append((leftB,rightB,B))
                elif slope == 1:
                    edges.append((leftA,rightB,B))
                elif slope == -1:
                    edges.append((rightA,leftB,B))
        G = Graph(edges, format='list_of_edges', loops=True, multiedges=True)
        subgraphs = G.connected_components_subgraphs()
        candidates = sorted(set(k for (u,v,k) in subgraph.edges()) for subgraph in subgraphs)
        if verbose:
            print('candidates=', candidates)

        ans = []
        for candidate in candidates:
            I = dominoes.intersection(itertools.product(candidate, repeat=2))
            if len(I) == 0:
                ans.append(candidate)
            else:
                if verbose:
                    print('rejecting {} since {} allows surrounding of '
                          'radius {}'.format(candidate, I, radius))
        return ans

    def find_markers(self, i=2, radius=1, solver=None, ncpus=1,
            verbose=False):
        r"""
        Return a list of lists of marker tiles.

        INPUT:

        - ``i`` -- integer (default:``2``), 1 or 2. 
        - ``radius`` - integer or 2-tuple (default: ``1``), if 2-tuple is
          given, then it is interpreted as ``(xradius, yradius)``
        - ``solver`` -- string (default:``None``)
        - ``ncpus`` -- integer (default:``1``)
        - ``verbose`` -- boolean (default:``False``)

        .. NOTE::

            The ``solver='dancing_links'`` is fast for this question.

        OUTPUT:

            list of lists

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(2,4,2,1), (2,2,2,0), (1,1,3,1), (1,2,3,2), (3,1,3,3),
            ....: (0,1,3,1), (0,0,0,1), (3,1,0,2), (0,2,1,2), (1,2,1,4), (3,3,1,2)]
            sage: tiles = [[str(a) for a in t] for t in tiles]
            sage: T = WangTileSet(tiles)
            sage: T.find_markers(i=1)
            []
            sage: T.find_markers(i=2)
            [[0, 1]]

        ::

            sage: tiles = ['ABCD', 'EFGH', 'AXCY', 'ABAB']
            sage: T = WangTileSet(tiles)
            sage: T.find_markers(i=1)
            [[0], [1], [2]]
            sage: T.find_markers(i=2)
            [[0], [1], [2]]

        """
        from sage.sets.disjoint_set import DisjointSet

        if not i in [1, 2]:
            raise ValueError("i(={}) should be 1 or 2".format(i))

        dominoes_i = self.dominoes_with_surrounding(i=i, radius=radius,
                                        solver=solver, ncpus=ncpus)
        dominoes_j = self.dominoes_with_surrounding(i=3-i, radius=radius,
                                        solver=solver, ncpus=ncpus)

        union_find = DisjointSet(len(self))
        for A,B in dominoes_j:
            union_find.union(A,B)

        ans = []
        for subset in union_find:
            I = dominoes_i.intersection(itertools.product(subset, repeat=2))
            if len(I) == 0:
                ans.append(subset)
            else:
                if verbose:
                    print('rejecting {} since {} allows surrounding of '
                          'radius {}'.format(subset, I, radius))
        return ans

    def find_substitution(self, M=None, i=2, side='right', radius=1,
            solver=None, ncpus=1, function=str.__add__, initial='',
            verbose=False):
        r"""
        Return the derived Wang tile set obtained from desubstitution using
        a given set of marker tiles.

        INPUT:

        - ``M`` -- markers, set of tile indices
        - ``i`` -- integer 1 or 2
        - ``side`` -- ``'right'`` or ``'left'``
        - ``radius`` - integer or 2-tuple (default: ``1``), if 2-tuple is
          given, then it is interpreted as ``(xradius, yradius)``
        - ``solver`` -- string (default:``None``)
        - ``ncpus`` -- integer (default:``1``)
        - ``function`` -- function (default:``str.__add__``), monoid
            operation
        - ``initial`` -- object (default:``''``), monoid neutral
        - ``verbose`` -- boolean

        OUTPUT:

            a 3-tuple (Wang tile set, substitution2d, set of markers)

        .. NOTE::

            The ``solver='dancing_links'`` is fast for this question.

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(2,4,2,1), (2,2,2,0), (1,1,3,1), (1,2,3,2), (3,1,3,3),
            ....: (0,1,3,1), (0,0,0,1), (3,1,0,2), (0,2,1,2), (1,2,1,4), (3,3,1,2)]
            sage: tiles = [[str(a) for a in t] for t in tiles]
            sage: T = WangTileSet(tiles)
            sage: T.find_markers(i=2)
            [[0, 1]]
            sage: T.find_substitution(M=[0,1], i=2)
            (Wang tile set of cardinality 12,
             Substitution 2d: {0: [[2]], 1: [[3]], 2: [[4]], 3:
                [[5]], 4: [[7]], 5: [[8]], 6: [[9]], 7: [[10]], 8:
                [[4, 0]], 9: [[5, 0]], 10: [[6, 1]], 11: [[7, 0]]})
        """
        # find markers
        if M is None:
            raise ValueError("a set of markers must be given as input, "
                    "markers can be computed with ``find_markers`` method")

        # Make sure M is of type set
        if not isinstance(M, set):
            M = set(M)

        if verbose:
            print("markers M =",M)

        # Compute the dominoes
        dominoes = self.dominoes_with_surrounding(i=i, radius=radius, solver=solver, ncpus=ncpus)
        if verbose:
            print("dominoes =",dominoes)

        # Print a warning when there are some valid M x M dominoes
        MM = [(a,b) for (a,b) in dominoes if a in M and b in M]
        if MM:
            print("Warning: it is expected as hypothesis that M odot^i M is "
                  "forbidden but the following dominoes admit a radius "
                  "{} neighborhood: {}. The algorihm works if M are "
                  " a set of markers anyway.".format(radius, MM))

        # Compute K and dominoes ending in M
        dominoes_without_M = [(a,b) for (a,b) in dominoes if a not in M and b not in M]
        if side == 'right':
            K = [a for (a,b) in dominoes_without_M]
            dominoes_M = [(a,b) for (a,b) in dominoes if b in M]
        elif side == 'left':
            K = [b for (a,b) in dominoes_without_M]
            dominoes_M = [(a,b) for (a,b) in dominoes if a in M]
        else:
            raise ValueError("side(={}) must be 'left' or 'right'".format(side))
        K = sorted(set(K))

        # We keep tiles in K
        tiles_and_structure = []
        for k in K:
            tiles_and_structure.append((self[k], [k]))

        # We add fusion of dominoes starting/ending with a tile in M
        for (a,b) in dominoes_M:
            t = fusion(self[a],self[b],i,function=function,initial=initial)
            tiles_and_structure.append((t, [a,b]))

        from slabbe.finite_word import sort_word_by_length_lex_key
        tiles_and_structure.sort(key=lambda v:sort_word_by_length_lex_key(v[1]))
        new_tiles, images = zip(*tiles_and_structure)
        d = {i:image for i,image in enumerate(images)}

        from slabbe.substitution_2d import Substitution2d
        if i == 1:
            s = Substitution2d.from_1d_row_substitution(d)
        elif i == 2:
            s = Substitution2d.from_1d_column_substitution(d)
        else:
            raise ValueError('i(={}) must be 1 or 2'.format(i))

        return WangTileSet(new_tiles), s

    def shear(self, radius=0, solver=None, ncpus=1,
            function=str.__add__, verbose=False):
        r"""
        Shears the Wang Tile set by the ``matrix(2,(1,-1,0,1))``.

        It is currently not implemented for other matrices.

        INPUT:

        - ``radius`` - integer or 2-tuple (default: ``0``), if 2-tuple is
          given, then it is interpreted as ``(xradius, yradius)``
        - ``solver`` -- string (default:``None``)
        - ``ncpus`` -- integer (default:``1``)
        - ``function`` -- function (default:``str.__add__``), monoid
          operation
        - ``verbose`` -- boolean (default:``False``)

        OUTPUT:

        - (WangTileSet, Substitution2d)

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [('aa','bb','cc','bb'), ('cc','dd','aa','dd')]
            sage: T = WangTileSet(tiles)
            sage: U,s = T.shear()
            sage: s
            Substitution 2d: {0: [[0]], 1: [[1]]}
            sage: U.tiles()
            [('aadd', 'dd', 'ccbb', 'bb'), ('ccbb', 'bb', 'aadd', 'dd')]
            sage: T.shear()[0].shear()[0].tiles()
            [('aaddbb', 'bb', 'ccbbdd', 'bb'), ('ccbbdd', 'dd', 'aaddbb', 'dd')]

        ::

            sage: tiles = [('aa','bb','cc','bb'), ('aa','dd','cc','bb'), ('cc','dd','aa','dd')]
            sage: T = WangTileSet(tiles)
            sage: U,s = T.shear()
            sage: s
            Substitution 2d: {0: [[0]], 1: [[1]], 2: [[2]], 3: [[2]]}
            sage: sorted(U.tiles())
            [('aadd', 'dd', 'ccbb', 'bb'),
             ('aadd', 'dd', 'ccdd', 'bb'),
             ('ccbb', 'bb', 'aadd', 'dd'),
             ('ccdd', 'dd', 'aadd', 'dd')]
            sage: U,s = T.shear(radius=1)
            sage: s
            Substitution 2d: {0: [[0]], 1: [[2]]}
            sage: U.tiles()
            [('aadd', 'dd', 'ccbb', 'bb'), ('ccbb', 'bb', 'aadd', 'dd')]

        """
        # Compute the dominoes
        dominoes = self.dominoes_with_surrounding(i=1, radius=radius, solver=solver, ncpus=ncpus)
        if verbose:
            print("dominoes =",dominoes)

        G = defaultdict(set)
        for (u,v) in dominoes:
            U = self[u]
            (e,n,w,s) = self[v]
            G[u].add(n)
        if verbose:
            print("The map Tile -> North[0] is {}".format(dict(G)))

        tiles = []
        sub = {}
        for i,(e,n,w,s) in enumerate(self):
            for a in G[i]:
                tile = (function(e,a),a,function(w,n),s)
                if tile not in tiles:
                    sub[len(tiles)] = i
                    tiles.append(tile)
        from slabbe.substitution_2d import Substitution2d
        return WangTileSet(tiles), Substitution2d.from_permutation(sub)

    def is_equivalent(self, other, certificate=False, verbose=False):
        r"""
        Return whether self and other are equivalent.

        INPUT:

        - ``other`` -- wang tile set
        - ``certificate`` -- boolean (default:``False``)
        - ``verbose`` -- boolean (default:``False``)

        .. NOTE::

            This code depends on the following bug to be fixed:
            https://trac.sagemath.org/ticket/24964

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(1,6,1,8), (2,6,1,7), (3,7,1,6), (1,6,2,6),
            ....:          (2,8,2,7), (2,7,3,6), (3,6,3,7)]
            sage: T = WangTileSet(tiles)
            sage: d = {1:'a', 2:'b', 3:'c', 6:'x', 7:'y', 8:'z'}
            sage: L = [tuple(d[a] for a in t) for t in tiles]
            sage: U = WangTileSet(L)
            sage: T.is_equivalent(U)
            True
            sage: T.is_equivalent(U,certificate=True)
            (True,
             {1: 'a', 2: 'b', 3: 'c'},
             {6: 'x', 7: 'y', 8: 'z'},
             Substitution 2d: {0: [[0]], 1: [[1]], 2: [[2]], 3: [[3]], 4: [[4]], 5: [[5]], 6: [[6]]})

        Not equivalent example::

            sage: _ = L.pop()
            sage: U = WangTileSet(L)
            sage: T.is_equivalent(U)
            False
            sage: T.is_equivalent(U,certificate=True)
            (False, None, None, None)

        When graphs admits non trivial automorphisms::

            sage: T = WangTileSet([(1,3,0,2), (0,2,1,3)])
            sage: U = WangTileSet([(7,'c',6,'z'), (6,'z',7,'c')])
            sage: V = WangTileSet([(7,9,6,8), (6,8,7,9)])
            sage: W = WangTileSet([(7,8,6,9), (6,9,7,8)])
            sage: T.is_equivalent(T, certificate=True)
            (True, {0: 0, 1: 1}, {2: 2, 3: 3}, Substitution 2d: {0: [[0]], 1: [[1]]})
            sage: T.is_equivalent(U, certificate=True)
            (True, {0: 6, 1: 7}, {2: 'z', 3: 'c'}, Substitution 2d: {0: [[0]], 1: [[1]]})
            sage: T.is_equivalent(V, certificate=True)
            (True, {0: 6, 1: 7}, {2: 8, 3: 9}, Substitution 2d: {0: [[0]], 1: [[1]]})
            sage: T.is_equivalent(W, certificate=True)
            (True, {0: 6, 1: 7}, {2: 9, 3: 8}, Substitution 2d: {0: [[0]], 1: [[1]]})
            sage: T.is_equivalent(W, certificate=True, verbose=True)
            True V_perm= {0: 6, 1: 7}
            True H_perm= {2: 8, 3: 9}
            Found automorphisms p=() and q=(2,3)
            (True, {0: 6, 1: 7}, {2: 9, 3: 8}, Substitution 2d: {0: [[0]], 1: [[1]]})

        """
        # Compute an isomorphism on the Vertical colors
        G = self.to_transducer_graph(merge_multiedges=False)
        H = other.to_transducer_graph(merge_multiedges=False)
        is_iso, V_perm = G.is_isomorphic(H, certificate=True)
        if verbose:
            print(is_iso, "V_perm=", V_perm)
        if not is_iso:
            return (False, None, None, None) if certificate else False

        # Compute an isomorphism on the Horizontal colors
        Gd = self.dual().to_transducer_graph(merge_multiedges=False)
        Hd = other.dual().to_transducer_graph(merge_multiedges=False)
        is_iso, H_perm = Gd.is_isomorphic(Hd, certificate=True)
        if verbose:
            print(is_iso, "H_perm=", H_perm)
        if not is_iso:
            return (False, None, None, None) if certificate else False

        # Compute the automorphisms group on Vertical and Horizontal colors
        # of self
        automorphism_group_V = G.automorphism_group()
        automorphism_group_H = Gd.automorphism_group()

        # If we find a pair (p,q) of automorphisms which works, return the result
        for p,q in itertools.product(automorphism_group_V, automorphism_group_H):
            perm_self_tiles = sorted((V_perm[p(a)], H_perm[q(b)],
                                      V_perm[p(c)], H_perm[q(d)])
                                     for (a,b,c,d) in self)
            sorted_other_tiles = sorted((a,b,c,d) for (a,b,c,d) in other)
            if sorted_other_tiles == perm_self_tiles:
                if certificate:
                    if verbose:
                        print ("Found automorphisms p={} and q={}".format(p,q))
                    V_perm_p = {a:V_perm[p(a)] for a in V_perm}
                    H_perm_q = {a:H_perm[q(a)] for a in H_perm}

                    other_to_index = {t:i for i,t in enumerate(other)}
                    d = {i:other_to_index[(V_perm_p[a],H_perm_q[b],V_perm_p[c],H_perm_q[d])] 
                           for i,(a,b,c,d) in enumerate(self)}
                    from slabbe.substitution_2d import Substitution2d
                    tiles_perm = Substitution2d.from_permutation(d)
                    return True, V_perm_p, H_perm_q, tiles_perm
                else:
                    return True
        return (False, None, None, None) if certificate else False

    def is_equivalent_up_to_isometry(self, other, certificate=False):
        r"""
        Return whether self and other are equivalent up to isometry.

        INPUT:

        - ``other`` -- wang tile set
        - ``certificate`` -- boolean (default:``False``)

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [(1,6,1,8), (2,6,1,7), (3,7,1,6), (1,6,2,6),
            ....:          (2,8,2,7), (2,7,3,6), (3,6,3,7)]
            sage: T = WangTileSet(tiles)
            sage: d = {1:'a', 2:'b', 3:'c', 6:'x', 7:'y', 8:'z'}
            sage: L = [tuple(d[a] for a in t) for t in tiles]
            sage: p = DihedralGroup(4)[6]
            sage: pL = [p(tuple(d[a] for a in t)) for t in tiles]
            sage: U = WangTileSet(pL)
            sage: T.is_equivalent(U)
            False
            sage: T.is_equivalent_up_to_isometry(U)
            True
            sage: T.is_equivalent_up_to_isometry(U, certificate=True)
            (True,
             ((1,4)(2,3),
              {6: 'x', 7: 'y', 8: 'z'},
              {1: 'a', 2: 'b', 3: 'c'},
              Substitution 2d: {0: [[0]], 1: [[1]], 2: [[2]], 3: [[3]], 4: [[4]], 5: [[5]], 6: [[6]]}))
        """
        from sage.groups.perm_gps.permgroup_named import DihedralGroup
        G = DihedralGroup(4)
        for p in G:
            p_self = WangTileSet([p(tile) for tile in self._tiles])
            is_equiv, V_perm, H_perm, tiles_perm = p_self.is_equivalent(other, certificate=True)
            if is_equiv:
                return (True, (p,V_perm,H_perm,tiles_perm)) if certificate else True
        else:
            return (False, (None,None,None,None)) if certificate else False

    def unsynchronized_graph_size2(self, i=1):
        r"""
        INPUT:

        - ``i`` -- integer, 1 or 2

        Signification of the nodes (u,v,w,d)::

             d = 0     |w| = d > 0      -|w| = d < 0

               |               |           |
              v|              v|          v|
               |            w  |           |
               +         +-----+           +-----+
               |         |                   w   |
              u|        u|                      u|
               |         |                       |

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [('aa','bb','cc','bb'), ('cc','dd','aa','dd')]
            sage: T = WangTileSet(tiles)
            sage: G = T.unsynchronized_graph_size2()       # known bug (in Python 3)
            sage: sorted(G.vertices())                     # known bug (in Python 3)
            [('aa', 'aa', '', 0), ('cc', 'cc', '', 0)]

        """
        if i == 2:
            raise NotImplementedError
        elif not i == 1:
            raise ValueError

        # Prepare the seeds
        seeds = []
        for (a,b) in self.dominoes_with_surrounding(i=2):
            (a0,a1,a2,a3) = self._tiles[a]
            (b0,b1,b2,b3) = self._tiles[b]
            assert a1 == b3
            seed = (a0,b0,'',0)
            seeds.append((seed,None))

        # Classify the tiles by their left color
        tile_dict = defaultdict(list)
        for i,t in enumerate(self):
            right,top,left,bottom = t
            tile_dict[left].append((i,t))

        # Define the children function
        def children(node):
            (u,v,w,d),label = node
            L = []
            if d == 0:
                for i,(a0,a1,a2,a3) in tile_dict[u]:
                    for j,(b0,b1,b2,b3) in tile_dict[v]:
                        if a1.startswith(b3):
                            new_node = (a0,b0,a1[len(b3):],len(b3)-len(a1))
                            L.append((new_node,(i,j)))
                        elif b3.startswith(a1):
                            new_node = (a0,b0,b3[len(a1):],len(b3)-len(a1))
                            L.append((new_node,(i,j)))
            elif d > 0:
                for i,(right,top,left,bottom) in tile_dict[u]:
                    if top.startswith(w):
                        new_node = (right,v,top[len(w):],d-len(top))
                        L.append((new_node,(i,None)))
                    elif w.startswith(top):
                        new_node = (right,v,w[len(top):],d-len(top))
                        L.append((new_node,(i,None)))
            elif d < 0:
                for j,(right,top,left,bottom) in tile_dict[v]:
                    if bottom.startswith(w):
                        new_node = (u,right,bottom[len(w):],d+len(bottom))
                        L.append((new_node,(None,j)))
                    elif w.startswith(bottom):
                        new_node = (u,right,w[len(bottom):],d+len(bottom))
                        L.append((new_node,(None,j)))
            else:
                assert False
            return L

        from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
        from slabbe.graph import digraph_move_label_to_edge
        from slabbe.graph import clean_sources_and_sinks
        R = RecursivelyEnumeratedSet(seeds, children, structure=None)
        G = R.to_digraph()
        GG = digraph_move_label_to_edge(G)
        H = clean_sources_and_sinks(GG)
        return H


    @cached_method
    def unsynchronized_graph(self, i=1, size=2, verbose=False):
        r"""
        INPUT:

        - ``i`` -- integer, 1 or 2
        - ``size`` -- integer, 2 or more
        - ``verbose`` -- boolean (default:``False``)

        OUTPUT:

        - graph of vertices (delays, blocks)

        Signification of the nodes (d,b)::

                 +-----------+
                 |           |
                 |    b[1]   |
                 |           |
            +----+-----+-----+
            |          |     |
            |   b[0]   |    d[1]
            |          |
            +----------+
                       |
                      d[0]

        EXAMPLES::

            sage: from slabbe import WangTileSet
            sage: tiles = [('aa','bb','cc','bb'), ('cc','dd','aa','dd')]
            sage: T = WangTileSet(tiles)
            sage: G = T.unsynchronized_graph()           # known bug (in Python 3)
            sage: sorted(G.vertices())                            # known bug
            [(d=(0, 0), b=(0, 0)),
             (d=(0, 0), b=(1, 1)),
             (d=(2, 0), b=(0, 1)),
             (d=(2, 0), b=(1, 0))]
            sage: G.edges()                                       # known bug
            [((d=(0, 0), b=(1, 1)), (d=(2, 0), b=(0, 1)), (1, 0)),
             ((d=(0, 0), b=(0, 0)), (d=(2, 0), b=(1, 0)), (0, 1)),
             ((d=(2, 0), b=(1, 0)), (d=(0, 0), b=(1, 1)), (0, 1)),
             ((d=(2, 0), b=(0, 1)), (d=(0, 0), b=(0, 0)), (1, 0))]
            sage: [node.lengths_x() for node in G]                # known bug (in Python 3)
            [[2, 2], [2, 2], [2, 2], [2, 2]]
            sage: [node.is_synchronized() for node in G]          # known bug (in Python 3)
            [True, True, True, True]
            sage: from slabbe import TikzPicture
            sage: _ = TikzPicture.from_graph(G).pdf(view=False)   # known bug (in Python 3)

        """
        from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet

        if i == 2:
            raise NotImplementedError
        elif not i == 1:
            raise ValueError

        class Node(object):
            def __init__(_self, delays, blocks):
                _self._delays = delays
                _self._blocks = blocks
            def __repr__(_self):
                return "(d={}, b={})".format(_self._delays, _self._blocks)
            def __hash__(_self):
                return hash(_self._delays) + hash(_self._blocks)
            def __eq__(_self, other):
                return (_self._delays == other._delays and
                        _self._blocks == other._blocks)
            def _latex_(_self):
                lines = []
                lines.append(r'\begin{tikzpicture}')
                for y in range(len(_self._delays)):
                    d = _self._delays[y]
                    b = _self._blocks[y]
                    tile = self[b]
                    sizex = len(tile[1])
                    position = (d-sizex, y)
                    new_lines = tile_to_tikz(tile, position, color=None,
                            id=b, id_color='', id_format='{}', sizex=sizex,
                            sizey=1, rotate=None, label=True,
                            label_shift=.2, label_color='black',
                            draw_H=None, draw_V=None)
                    lines.extend(new_lines)
                lines.append(r'\end{tikzpicture}')
                return '\n'.join(lines)
            def lengths_x(_self):
                return [len(self[b][1]) for b in _self._blocks]
            def is_synchronized(_self):
                delays = _self._delays
                lengths_x = _self.lengths_x()
                for i in range(len(delays)-1):
                    da,db = delays[i:i+2]
                    la,lb = lengths_x[i:i+2]
                    if da == db:
                        if la == lb:
                            pass
                        else:
                            return  False
                    elif da < db:
                        if db - da == lb:
                            pass
                        else:
                            return False
                    else:
                        if da - db == la:
                            pass
                        else:
                            return False
                return True

        # Preparing the seeds
        if verbose:
            print("Computing the seeds of size {} ...".format(size))
        if size == 2:
            seeds = []
            for (a,b) in self.dominoes_with_surrounding(i=2):
                (a0,a1,a2,a3) = self._tiles[a]
                (b0,b1,b2,b3) = self._tiles[b]
                assert a1 == b3
                delays = (0,0)
                blocks = (a,b)
                node = Node(delays, blocks)
                seeds.append((node,None))
        elif size > 2:
            from slabbe.finite_word import are_overlapping_factors
            G = self.unsynchronized_graph(i=i, size=size-1, verbose=verbose)
            seeds = set()
            for node in G:
                delays = node._delays
                blocks = node._blocks
                assert len(delays) == len(blocks)
                u = blocks[-1]
                (u0,u1,u2,u3) = self[u]
                du = delays[-1]
                for z in range(len(self)):
                    (z0,z1,z2,z3) = self[z]
                    for dz in range(len(z3)+1):
                        if are_overlapping_factors(z3, u1, len(z3)-dz+du-len(u1)):
                            delays_copy = list(delays)
                            blocks_copy = list(blocks)
                            delays_copy.append(dz)
                            blocks_copy.append(z)
                            assert min(delays_copy) == 0
                            node = Node(tuple(delays_copy), tuple(blocks_copy))
                            ## TODO: add the node only if the suffix of
                            ## length - 1 is in G
                            seeds.add((node,None))
        else:
            raise ValueError("size(={}) is bad input".format(size))
        if verbose:
            print("Set of {} seeds of size {} found".format(len(seeds), size))

        # Classify the tiles by their left color
        # and the prefixes of their top and bottom colors
        if verbose:
            print("Classifying the tiles...")
        tile_dict = defaultdict(list)
        for i,t in enumerate(self):
            right,top,left,bottom = t
            for j,k in itertools.product(range(len(top)+1), repeat=2):
                top_prefix = top[:j]
                bottom_prefix = bottom[:k]
                tile_dict[(top_prefix, left, bottom_prefix)].append(i)
        tile_dict = dict(tile_dict)
        #print("tile_dict=",tile_dict)

        # Define the children function
        if verbose:
            print("Define the children function...")
        def children(node_label):
            node,label = node_label
            delays = node._delays
            blocks = node._blocks
            assert min(delays) == 0
            index = delays.index(0)
            u = blocks[index]
            (u0,u1,u2,u3) = self[u]
            if index > 0:
                # there is some block below
                v = blocks[index-1]
                (v0,v1,v2,v3) = self[v]
                length_below = delays[index-1]
                assert 0 < length_below <= len(v1)
                v1_suffix = v1[-length_below:]
            else:
                v1_suffix = ''
            if index + 1 < size:
                # there is some block above
                w = blocks[index+1]
                (w0,w1,w2,w3) = self[w]
                length_above = delays[index+1]
                assert 0 <= length_above <= len(w3)
                w3_suffix = w3[len(w3)-length_above:]
            else:
                w3_suffix = ''

            L = []
            for z in tile_dict.get((w3_suffix, u0, v1_suffix), []):
                (z0,z1,z2,z3) = self[z]
                assert len(z1) == len(z3)

                blocks_copy = list(blocks)
                blocks_copy[index] = z
                blocks_copy = tuple(blocks_copy)

                delays_copy = list(delays)
                delays_copy[index] += len(z1)
                m = min(delays_copy)
                delays_copy = tuple(a-m for a in delays_copy)

                node = Node(delays_copy,blocks_copy)
                label = (u,z)
                L.append((node,label))
            return L

        if verbose:
            print("Enumeration of the elements of size {} ...".format(size))
        R = RecursivelyEnumeratedSet(seeds, children, structure=None)
        G = R.to_digraph()

        from slabbe.graph import digraph_move_label_to_edge
        from slabbe.graph import clean_sources_and_sinks
        G = digraph_move_label_to_edge(G)
        H = clean_sources_and_sinks(G)

        if verbose:
            print("Digraph computed:",G)
            print("Digraph computed (after removing sources and sinks):",H)

        return H

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

        ::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
            sage: t = {(0,2):0}
            sage: c = [{},{},{(2,2):0},{}]
            sage: W = WangTileSolver(tiles,3,3,preassigned_tiles=t,preassigned_color=c)

        TESTS::

            sage: t = {(0,2):0}
            sage: c = [{},{},{(2,3):0},{}]
            sage: W = WangTileSolver(tiles,3,3,preassigned_tiles=t,preassigned_color=c)
            Traceback (most recent call last):
            ...
            AssertionError

        """
        self._tiles = list(tiles)
        self._width = width
        self._height = height
        if preassigned_color is None:
            preassigned_color = [{}, {}, {}, {}]
        self._preassigned_color = preassigned_color
        if preassigned_tiles is None:
            preassigned_tiles = {}
        self._preassigned_tiles = preassigned_tiles
        self._color = color

        assert all(0 <= j < self._width for (j,k) in self._preassigned_tiles)
        assert all(0 <= k < self._height for (j,k) in self._preassigned_tiles)
        assert all(0 <= j < self._width for d in self._preassigned_color for (j,k) in d)
        assert all(0 <= k < self._height for d in self._preassigned_color for (j,k) in d)

    def vertical_alphabet(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles, 3, 3)
            sage: W.vertical_alphabet()
            {0, 1, 2}
        r"""
        right, top, left, bottom = zip(*self._tiles)
        return set(left) | set(right)

    def horizontal_alphabet(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles, 3, 3)
            sage: W.horizontal_alphabet()
            {0, 1, 2}
        r"""
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
            sage: sorted(support)
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 1, 0), (0, 1, 1), 
             (0, 1, 2), (0, 1, 3), (0, 2, 0), (0, 2, 1), (0, 2, 2), (0, 2, 3)]

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
          ``'GLPK'``, ``'Coin'``, ``'cplex'`` or ``'Gurobi'``.
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
            ValueError: no solution found using dancing links, the return
            value from dancing links solver is None

        Using SatLP solver::

            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve('LP')
            sage: tiling
            A wang tiling of a 3 x 4 rectangle

        Using SatLP solver with preassigned tiles::

            sage: preassigned = {(0,0):0}
            sage: W = WangTileSolver(tiles,3,4,preassigned_tiles=preassigned)
            sage: tiling = W.solve(solver='LP')
            sage: tiling._table
            [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]

        Using cryptominisat solver::

            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve('cryptominisat')  # optional cryptominisat
            sage: tiling._table                      # optional cryptominisat
            [[1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]]

        REFERENCES:

            How do I set solver_parameter to make Gurobi use more than one
            processor?, https://ask.sagemath.org/question/37726/
        """
        if solver == 'dancing_links':
            from sage.combinat.matrices.dancing_links import dlx_solver
            rows,row_info = self.rows_and_information()
            dlx = dlx_solver(rows)
            if ncpus == 1:
                solution = dlx.get_solution() if dlx.search() else None
            else:
                solution = dlx.one_solution(ncpus=ncpus)
            if solution is None or isinstance(solution, str):
                raise ValueError('no solution found using dancing links,'
                                 ' the return value from dancing links'
                                 ' solver is {}'.format(solution))
            table = [[None]*self._height for _ in range(self._width)]
            for a in solution:
                j,k,i = row_info[a]
                assert table[j][k] is None, "table[{}][{}](={}) is not None".format(j,k,table[j][k])
                table[j][k] = i
            return WangTiling(table, self._tiles, color=self._color)
        elif solver in ['Gurobi', 'gurobi', 'GLPK', 'cplex', 'Coin', 'CVXOPT', 'PPL', None]:
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
        else: # we assume we use a sat solver
            (var_to_tile_pos,
             tile_pos_to_var) = self.sat_variable_to_tile_position_bijection()
            sat_solver = self.sat_solver(solver)
            solution = sat_solver()
            if not solution:
                raise ValueError('no solution found using SAT solver (={})'.format(solver))
            support = [key for (key,val) in enumerate(solution) if val]
            assert len(support) == self._width * self._height, ("len(support)={} "
                    "!= width*height".format(len(support)))
            table = [[None]*self._height for _ in range(self._width)]
            for val in support:
                i,j,k = var_to_tile_pos[val]
                table[j][k] = i
            return WangTiling(table, self._tiles, color=self._color)

    def has_solution(self, solver=None, solver_parameters=None, ncpus=1):
        r"""
        Return whether there is a solution.

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

        Using cryptominisat::

            sage: W = WangTileSolver(tiles,3,4)
            sage: W.has_solution(solver='cryptominisat') # optional cryptominisat
            True
        """
        if solver == 'dancing_links':
            from sage.combinat.matrices.dancing_links import dlx_solver
            rows,row_info = self.rows_and_information()
            dlx = dlx_solver(rows)
            if ncpus == 1:
                return dlx.search() == 1
            else:
                solution = dlx.one_solution(ncpus=ncpus)
                return solution is not None
        elif solver in ['Gurobi', 'gurobi', 'GLPK', 'Coin', 'CVXOPT', 'PPL', None]:
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
        else:
            solution = self.sat_solver(solver)()
            return bool(solution)

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

        With preassigned colors::

            sage: right = {(0, 1): 'A', (0, 0): 'A'}
            sage: top = {(0, 1): 'B'}
            sage: left = {(0, 1): 'A', (0, 0): 'A'}
            sage: bottom = {(0, 0): 'B'}
            sage: preassigned_color=[right,top,left,bottom]
            sage: tiles = ['ABCD', 'EFGH', 'AXCY', 'ABAB', 'EBEB']
            sage: W = WangTileSolver(tiles, 1, 2, preassigned_color=preassigned_color)
            sage: W.rows_and_information()
            ([[4], [4], [4], [1, 2, 3, 4], [4], [5], [5], [5], [0, 5], [5]],
             [(0, 0, 0),
              (0, 0, 1),
              (0, 0, 2),
              (0, 0, 3),
              (0, 0, 4),
              (0, 1, 0),
              (0, 1, 1),
              (0, 1, 2),
              (0, 1, 3),
              (0, 1, 4)])

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

        def is_forbidden_at_position(tile, j, k):
            right,top,left,bottom = tile
            return ((j,k) in rightPRE and rightPRE[(j,k)] != right
                 or (j,k) in leftPRE and leftPRE[(j,k)] != left
                 or (j,k) in bottomPRE and bottomPRE[(j,k)] != bottom
                 or (j,k) in topPRE and topPRE[(j,k)] != top)

        # matching vertical colors
        for j in range(W):
            for k in range(H):
                position = (k*(W-1)+j)*padtoV
                for i,tile in enumerate(self._tiles):
                    if is_forbidden_at_position(tile, j, k):
                        continue
                    # the tile i at position (j,k)
                    right,top,left,bottom = tile
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
                    if is_forbidden_at_position(tile, j, k):
                        continue
                    # the tile i at position (j,k)
                    right,top,left,bottom = tile
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

    def sat_variable_to_tile_position_bijection(self):
        r"""
        Return the dictionary giving the correspondence between variables
        and tiles indices i at position (j,k)

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: d1,d2 = W.sat_variable_to_tile_position_bijection()
            sage: d1
            {1: (0, 0, 0),
             2: (0, 0, 1),
             3: (0, 0, 2),
             4: (0, 0, 3),
             5: (0, 1, 0),
             6: (0, 1, 1),
             7: (0, 1, 2),
             8: (0, 1, 3),
             9: (0, 2, 0),
             10: (0, 2, 1),
             11: (0, 2, 2),
             12: (0, 2, 3),
             13: (1, 0, 0),
             14: (1, 0, 1),
             15: (1, 0, 2),
             16: (1, 0, 3),
             17: (1, 1, 0),
             18: (1, 1, 1),
             19: (1, 1, 2),
             20: (1, 1, 3),
             21: (1, 2, 0),
             22: (1, 2, 1),
             23: (1, 2, 2),
             24: (1, 2, 3)}

        """
        W = self._width
        H = self._height
        ntiles = len(self._tiles)
        L = list(itertools.product(range(ntiles), range(W), range(H)))
        var_to_tile_pos = dict(enumerate(L, start=1))
        tile_pos_to_var = dict((b,a) for (a,b) in enumerate(L, start=1))
        return var_to_tile_pos, tile_pos_to_var

    def sat_solver(self, solver=None):
        r"""
        Return the SAT solver.

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: s = W.sat_solver()
            sage: s           # random
            an ILP-based SAT Solver
            CryptoMiniSat solver: 24 variables, 58 clauses.
            sage: L = s()
            sage: list(L)
            [None, ...]

        """
        tiles = self._tiles
        indices = range(len(tiles))

        from sage.sat.solvers.satsolver import SAT
        s = SAT(solver)

        (var_to_tile_pos,
         tile_pos_to_var) = self.sat_variable_to_tile_position_bijection()

        # at least one tile at each position (j,k)
        # (exactly one if one could use a xor clause)
        for j in range(self._width):
            for k in range(self._height):
                constraint = [tile_pos_to_var[(i,j,k)] for i in indices]
                s.add_clause(constraint)

        # no two tiles at the same position (j,k)
        for j in range(self._width):
            for k in range(self._height):
                for i1,i2 in itertools.combinations(indices, 2):
                    constraint = [-tile_pos_to_var[(i1,j,k)], 
                                  -tile_pos_to_var[(i2,j,k)]]
                    s.add_clause(constraint)

        # preassigned tiles at position (j,k)
        for j,k in self._preassigned_tiles:
            i = self._preassigned_tiles[(j,k)]
            constraint = [tile_pos_to_var[(i,j,k)]]
            s.add_clause(constraint)

        # matching vertical colors
        for j in range(self._width-1):
            for k in range(self._height):
                for i1,i2 in itertools.product(indices, repeat=2):
                    if tiles[i1][0] != tiles[i2][2]:
                        constraint = [-tile_pos_to_var[(i1,j,k)], 
                                      -tile_pos_to_var[(i2,j+1,k)]]
                        s.add_clause(constraint)

        # matching horizontal colors
        for j in range(self._width):
            for k in range(self._height-1):
                for i1,i2 in itertools.product(indices, repeat=2):
                    if tiles[i1][1] != tiles[i2][3]:
                        constraint = [-tile_pos_to_var[(i1,j,k)], 
                                      -tile_pos_to_var[(i2,j,k+1)]]
                        s.add_clause(constraint)

        # matching preassigned color constraints
        legend = {0:'right',1:'top',2:'left',3:'bottom'}
        for angle, D in enumerate(self._preassigned_color):
            for j,k in D:
                for i in indices:
                    if tiles[i][angle] != D[(j,k)]:
                        constraint = [-tile_pos_to_var[(i,j,k)]]
                        s.add_clause(constraint)

        return s

    def number_of_solutions(self, ncpus=8):
        r"""
        Return the number of solutions

        INPUT:

        - ``ncpus`` -- integer (default: ``8``), maximal number of
          subprocesses to use at the same time

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

        With preassigned colors and tiles::

            sage: right = {(0, 1): 'A', (0, 0): 'A'}
            sage: top = {(0, 1): 'B'}
            sage: left = {(0, 1): 'A', (0, 0): 'A'}
            sage: bottom = {(0, 0): 'B'}
            sage: preassigned_color=[right,top,left,bottom]
            sage: tiles = ['ABCD', 'EFGH', 'AXCY', 'ABAB', 'EBEB']
            sage: W = WangTileSolver(tiles, 1, 2, preassigned_color=preassigned_color)
            sage: solutions = list(W.solutions_iterator())
            sage: [t.table() for t in solutions]
            [[[3, 3]]]

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

    def all_solutions(self, ncpus=8):
        r"""
        Return the list of all solutions.

        .. NOTE::

            This uses the reduction to dancing links.

        INPUT:

        - ``ncpus`` -- integer (default: ``8``), maximal number of
          subprocesses to use at the same time

        OUTPUT:

            list of wang tilings

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(2,4,2,1), (2,2,2,0), (1,1,3,1), (1,2,3,2), (3,1,3,3),
            ....: (0,1,3,1), (0,0,0,1), (3,1,0,2), (0,2,1,2), (1,2,1,4), (3,3,1,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: W.number_of_solutions()
            908
            sage: L = W.all_solutions()
            sage: len(L)
            908

        ::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles,2,2)
            sage: W.all_solutions()
            [A wang tiling of a 2 x 2 rectangle,
             A wang tiling of a 2 x 2 rectangle,
             A wang tiling of a 2 x 2 rectangle]

        With preassigned colors and tiles::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
            sage: t = {(0,1):0}
            sage: c = [{},{},{(1,1):0},{}]
            sage: W = WangTileSolver(tiles,3,3,preassigned_tiles=t,preassigned_color=c)
            sage: S = W.all_solutions()
            sage: sorted([s._table for s in S])
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[0, 0, 3], [0, 0, 0], [0, 0, 0]]]

        With preassigned colors and tiles::

            sage: right = {(0, 1): 'A', (0, 0): 'A'}
            sage: top = {(0, 1): 'B'}
            sage: left = {(0, 1): 'A', (0, 0): 'A'}
            sage: bottom = {(0, 0): 'B'}
            sage: preassigned_color=[right,top,left,bottom]
            sage: tiles = ['ABCD', 'EFGH', 'AXCY', 'ABAB', 'EBEB']
            sage: W = WangTileSolver(tiles, 1, 2, preassigned_color=preassigned_color)
            sage: [t.table() for t in W.all_solutions()]
            [[[3, 3]]]

        """
        from sage.combinat.matrices.dancing_links import dlx_solver
        rows,row_info = self.rows_and_information()
        dlx = dlx_solver(rows)
        L = []
        for solution in dlx.all_solutions(ncpus=ncpus):
            table = [[None]*self._height for _ in range(self._width)]
            for a in solution:
                j,k,i = row_info[a]
                assert table[j][k] is None, "table[{}][{}](={}) is not None".format(j,k,table[j][k])
                table[j][k] = i
            tiling = WangTiling(table, self._tiles, color=self._color)
            L.append(tiling)
        return L

    def all_solutions_tikz(self, ncpus=8):
        r"""
        INPUT:

        - ``ncpus`` -- integer (default: ``8``), maximal number of
          subprocesses to use at the same time

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(2,4,2,1), (2,2,2,0), (1,1,3,1), (1,2,3,2), (3,1,3,3),
            ....: (0,1,3,1), (0,0,0,1), (3,1,0,2), (0,2,1,2), (1,2,1,4), (3,3,1,2)]
            sage: W = WangTileSolver(tiles,2,2)
            sage: t = W.all_solutions_tikz()
            sage: view(t)    # long # not tested
        """
        solutions = self.all_solutions(ncpus=ncpus)
        L = [sol.tikz(scale=.7,fontsize=r'\tiny').tikz_picture_code() for sol in solutions]
        bigtikz = '\n'.join(L)
        from sage.misc.latex import LatexExpr
        return LatexExpr(bigtikz)

    def meet_of_all_solutions(self, ncpus=8):
        r"""
        Return the tiling of the rectangle with tiles that are imposed at
        each position (this is the meet of the partially ordered set of
        all partial solutions inside the rectangle).

        INPUT:

        - ``ncpus`` -- integer (default: ``8``), maximal number of
          subprocesses to use at the same time

        OUTPUT:

            A Wang tiling (with ``None`` at positions where more than one
            tile can occur)

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
            sage: t = {(0,1):0}
            sage: W = WangTileSolver(tiles,3,3,preassigned_tiles=t)
            sage: tiling = W.meet_of_all_solutions()
            sage: tiling
            A wang tiling of a 3 x 3 rectangle
            sage: tiling.table()
            [[0, 0, None], [0, 0, 0], [0, 0, 0]]

        ::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2), (0,1,2,0)]
            sage: W = WangTileSolver(tiles,3,3)
            sage: tiling = W.meet_of_all_solutions()
            sage: tiling.table()
            [[None, None, None], [None, None, None], [None, None, None]]

        """
        W = self._width
        H = self._height
        d = defaultdict(set)
        for tiling in self.all_solutions(ncpus=ncpus):
            for j in range(W):
                for k in range(H):
                    table = tiling.table()
                    tile = table[j][k]
                    d[(j,k)].add(tile)
        table = []
        for j in range(W):
            row = []
            for k in range(H):
                tiles = d[(j,k)]
                if len(tiles) == 1:
                    row.append(tiles.pop())
                else:
                    row.append(None)
            table.append(row)
        return WangTiling(table, self._tiles, color=self._color)

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

    def _matrix_(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: matrix(tiling)
            [1 0 1]
            [0 1 0]
            [1 0 1]
            [0 1 0]

        """
        from sage.matrix.constructor import matrix
        return matrix.column([col[::-1] for col in self._table])

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

    def rows(self):
        r"""
        Return the rows from the top to the bottom.

        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4)] * 10
            sage: table = [[1,2,3,4], [5,6,7,8]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.rows()
            [(4, 8), (3, 7), (2, 6), (1, 5)]

        """
        L = []
        for row in zip(*self.table()):
            L.insert(0, row)
        return L

    def to_image(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: img = tiling.to_image()
            sage: img
            <PIL.Image.Image image mode=RGB size=3x4 at ...>
            sage: img.resize((100,100)).show()     # not tested

        """
        # Chose a random color for each tiles
        from random import randrange
        from collections import defaultdict
        def random_color():
                return (randrange(255),randrange(255),randrange(255))
        color_dict = defaultdict(random_color)
        color_dict[None] = [0,0,0] # black (=0) as the color for None

        # Create the image
        import numpy as np
        from PIL import Image
        data = [tuple(color_dict[a] for a in row) for row in self.rows()]
        data = np.array(data, dtype=np.uint8)
        return Image.fromarray(data)#.resize((100,100))

    def transpose(self):
        r"""
        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: M = matrix(2, (1,1,0,1))
            sage: tiling_T = tiling.transpose()
            sage: tiling_T.table()
            [[0, 1, 0], [1, 0, 1], [0, 1, 0], [1, 0, 1]]

        """
        table = [[column[j] for column in self._table]
                            for j in range(self.height())]
        return WangTiling(table, self._tiles, self._color)

    def concatenate_right(self, other):
        r"""
        Return a concatenation of self with other on its right.

        INPUT:

        - ``other`` -- a Wang tiling

        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,2), (1,2,0,3)]
            sage: tableA = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
            sage: tableB = [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]
            sage: tilingA = WangTiling(tableA, tiles)
            sage: tilingB = WangTiling(tableB, tiles)
            sage: t = tilingA.concatenate_right(tilingB)
            sage: t
            A wang tiling of a 6 x 4 rectangle
            sage: t.table()
            [[0, 0, 0, 0],
             [0, 0, 0, 0],
             [0, 0, 0, 0],
             [1, 1, 1, 1],
             [1, 1, 1, 1],
             [1, 1, 1, 1]]

        """
        from copy import copy
        table = copy(self.table())
        table.extend(other.table())
        return WangTiling(table, self._tiles, color=self._color)

    def concatenate_top(self, other):
        r"""
        Return a concatenation of self with other on its top.

        INPUT:

        - ``other`` -- a Wang tiling

        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,2), (1,2,0,3)]
            sage: tableA = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
            sage: tableB = [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]
            sage: tilingA = WangTiling(tableA, tiles)
            sage: tilingB = WangTiling(tableB, tiles)
            sage: t = tilingA.concatenate_top(tilingB)
            sage: t
            A wang tiling of a 3 x 8 rectangle
            sage: t.table()
            [[0, 0, 0, 0, 1, 1, 1, 1], 
             [0, 0, 0, 0, 1, 1, 1, 1], 
             [0, 0, 0, 0, 1, 1, 1, 1]]

        """
        table = [colA+colB for colA,colB in zip(self.table(), other.table())]
        return WangTiling(table, self._tiles, color=self._color)

    def diff(self, other):
        r"""
        Return a Wang tiling where positions where the tile in self differs
        from the tile in other are replaced by ``None``.

        INPUT:

        - ``other`` -- a Wang tiling

        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,2), (1,2,0,3)]
            sage: tableA = [[0, 1, 1, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tilingA = WangTiling(tableA, tiles)
            sage: tableB = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 1, 0, 1]]
            sage: tilingB = WangTiling(tableB, tiles)
            sage: d = tilingA.diff(tilingB)
            sage: d
            A wang tiling of a 3 x 4 rectangle
            sage: d.table()
            [[0, 1, None, None], [1, 0, None, 0], [0, 1, 0, 1]]

        """
        table = []
        for colA,colB in zip(self.table(), other.table()):
            col = [(a if a==b else None) for a,b in zip(colA, colB)]
            table.append(col)
        return WangTiling(table, self._tiles, color=self._color)


    def apply_matrix_transformation(self, M):
        r"""
        INPUT:

        - ``M`` -- matrix in GL(2,Z)

        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: M = matrix(2, (1,1,0,1))
            sage: tiling_M = tiling.apply_matrix_transformation(M)
            sage: tiling_M.table()
            [[0, None, None, None],
             [1, 1, None, None],
             [0, 0, 0, None],
             [None, 1, 1, 1],
             [None, None, 0, 0],
             [None, None, None, 1]]

        """
        from sage.modules.free_module_element import vector
        D = defaultdict(lambda:None)
        for i,column in enumerate(self._table):
            for j,value in enumerate(column):
                new_pos = tuple(M*vector((i, j)))
                D[new_pos] = value
        D_keys = D.keys()
        X,Y = zip(*D_keys)
        min_X = min(X)
        max_X = max(X)
        min_Y = min(Y)
        max_Y = max(Y)
        table = []
        for x in range(min_X, max_X+1):
            column = [D[(x,y)] for y in range(min_Y, max_Y+1)]
            table.append(column)
        return WangTiling(table, self._tiles, self._color)

    def slide(self, shift, x0=None, y0=1):
        r"""
        INPUT:

        - ``shift`` -- integer
        - ``x0`` -- integer or None, every tile at (x,y) such that x>=x0
          will be shifted by (0,shift)
        - ``y0`` -- integer or None, every tile at (x,y) such that y>=y0
          will be shifted by (shift,0)

        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.slide(0).table()
            [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling.slide(3).table()
            [[0, None, None, None],
             [1, None, None, None],
             [0, None, None, None],
             [None, 1, 0, 1],
             [None, 0, 1, 0],
             [None, 1, 0, 1]]
            sage: tiling.slide(2, x0=2).table()
            [[0, 1, 0, 1, None, None], [1, 0, 1, 0, None, None], [None, None, 0, 1, 0, 1]]
            sage: tiling.slide(-2, x0=2).table()
            [[None, None, 0, 1, 0, 1], [None, None, 1, 0, 1, 0], [0, 1, 0, 1, None, None]]

        """
        if x0 is None:
            bottom = [col[:y0] for col in self._table]
            top = [col[y0:] for col in self._table]
            top_height = self.height()-y0
            if shift >= 0:
                top = [[None]*top_height]*shift + top
                bottom += [[None]*y0]*shift
            else:
                top += [[None]*top_height]*(-shift)
                bottom = [[None]*y0]*(-shift) + bottom
            table = [b+t for (b,t) in zip(bottom, top)]
            return WangTiling(table, self._tiles, self._color)
        else:
            return self.transpose().slide(shift, x0=None, y0=x0).transpose()


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

    def plot_points_on_torus(self, M, pointsize=5, color_dict=None,
            start=None):
        r"""
        Plot points modulo some values in x and y.

        INPUT

        - ``M`` -- M is the matrix projection to `\mathbb{R}^2/\mathbb{Z}^2`
        - ``pointsize`` -- positive real number (default:``5``)
        - ``color_dict`` -- dict, tile index -> color or None
          (default:``None``)
        - ``start`` -- None or vector

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
        if start is None:
            start = vector((0,0))
        else:
            start = vector(start)
        PTS = {a:[] for a in alphabet}
        for i,column in enumerate(self._table):
            for j,a in enumerate(column):
                x,y = M*vector((i,j)) - start
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
            G += text(str(a), (px,py), color='black')

        # color of points in legend
        for d,a in enumerate(alphabet):
            color = color_dict[a]
            G += points([(-.1,0)], color=color, size=pointsize, legend_label=str(d))

        return G

    def pattern_occurrences(self, pattern, avoid_border=0):
        r"""
        Return the set of occurrences of the given pattern in the tiling.

        INPUT

        - ``pattern`` -- dict
        - ``avoid_border`` -- integer (default: 0), the size of the border
          to avoid during the computation

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: tiling.pattern_occurrences({(0,0):0})
            {(0, 0), (0, 2), (1, 1), (1, 3), (2, 0), (2, 2)}
            sage: tiling.pattern_occurrences({(0,0):1})
            {(0, 1), (0, 3), (1, 0), (1, 2), (2, 1), (2, 3)}
            sage: tiling.pattern_occurrences({(0,0):1, (1,0):1})
            set()
            sage: tiling.pattern_occurrences({(0,0):1, (1,0):1, (0,1):1})
            set()
            sage: tiling.pattern_occurrences({(0,0):1, (1,0):0, (0,1):0})
            {(0, 1), (1, 0), (1, 2)}

        The positions depends on the relative position of the pattern::

            sage: tiling.pattern_occurrences({(0,-1):1})
            {(0, 2), (0, 4), (1, 1), (1, 3), (2, 2), (2, 4)}
            sage: tiling.pattern_occurrences({(-1,-1):1})
            {(1, 2), (1, 4), (2, 1), (2, 3), (3, 2), (3, 4)}
            sage: tiling.pattern_occurrences({(-100,-100):1})
            {(100, 101), (100, 103), (101, 100), (101, 102), (102, 101), (102, 103)}

        The x coordinates of the pattern corresponds to the x coordinates
        when you plot it::

            sage: tiles = [(0,3,0,4), (1,4,1,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: tiling.pattern_occurrences({(0,0):1})
            {(0, 1), (0, 3), (1, 1), (1, 3), (2, 1), (2, 3)}
            sage: tiling.pattern_occurrences({(0,0):1, (1,0):1})
            {(0, 1), (0, 3), (1, 1), (1, 3)}
            sage: tiling.pattern_occurrences({(0,0):1, (0,1):1})
            set()
            sage: tiling.tikz().pdf(view=False)   # not tested

        When avoiding the border::

            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: tiling.pattern_occurrences({(0,0):0}, avoid_border=1)
            {(1, 1)}
        """
        xmin = min(x for (x,y) in pattern)
        xmax = max(x for (x,y) in pattern)
        ymin = min(y for (x,y) in pattern)
        ymax = max(y for (x,y) in pattern)
        S = set()
        for i in range(0-xmin+avoid_border, self.width()-xmax-avoid_border):
            for j in range(0-ymin+avoid_border, self.height()-ymax-avoid_border):
                if all(self._table[i+x][j+y] == pattern[(x,y)] for (x,y) in pattern):
                    S.add((i,j))
        return S

    def number_of_occurrences(self, pattern, avoid_border=0):
        r"""
        Return the number of occurrences of the given pattern in the tiling.

        INPUT

        - ``pattern`` -- dict
        - ``avoid_border`` -- integer (default: 0), the size of the border
          to avoid during the computation

        EXAMPLES::

            sage: from slabbe import WangTileSolver
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: tiling.number_of_occurrences({(0,0):0})
            6
            sage: tiling.number_of_occurrences({(0,0):1})
            6
            sage: tiling.number_of_occurrences({(0,0):1, (1,0):1})
            0
            sage: tiling.number_of_occurrences({(0,0):1, (1,0):1, (0,1):1})
            0
            sage: tiling.number_of_occurrences({(0,0):1, (1,0):0, (0,1):0})
            3

        The pattern is translation invariant::

            sage: tiling.number_of_occurrences({(0,-1):1})
            6
            sage: tiling.number_of_occurrences({(-1,-1):1})
            6
            sage: tiling.number_of_occurrences({(-100,-100):1})
            6

        The x coordinates of the pattern corresponds to the x coordinates
        when you plot it::

            sage: tiles = [(0,3,0,4), (1,4,1,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: tiling.number_of_occurrences({(0,0):1})
            6
            sage: tiling.number_of_occurrences({(0,0):1, (1,0):1})
            4
            sage: tiling.number_of_occurrences({(0,0):1, (0,1):1})
            0
            sage: tiling.tikz().pdf(view=False)   # not tested

        When avoiding the border::

            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: tiling.number_of_occurrences({(0,0):0}, avoid_border=1)
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

    def number_of_occurrences_with_shape(self, shape, avoid_border=0):
        r"""
        Return the number of occurrences of every pattern having a given
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
            sage: tiling.number_of_occurrences_with_shape([(0,0)])
            Counter({(0,): 12})

        ::

            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.number_of_occurrences_with_shape([(0,0)])
            Counter({(0,): 6, (1,): 6})
            sage: c = tiling.number_of_occurrences_with_shape([(0,0), (1,0), (0,1)])
            sage: sorted(c.items())
            [((0, 1, 1), 3), ((1, 0, 0), 3)] 

        When avoiding the border::

            sage: tiling.number_of_occurrences_with_shape([(0,0)], avoid_border=1)
            Counter({(0,): 1, (1,): 1})
            sage: tiling.number_of_occurrences_with_shape([(0,0)], avoid_border=2)
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
        C = self.number_of_occurrences_with_shape([(0,0)], avoid_border=avoid_border)
        s = sum(C.values())
        return {k:QQ((v,s)) for k,v in C.items()}

    def tile_positions(self, M):
        r"""
        Return the list of positions where tile of M appear.

        INPUT:

        - ``M`` -- subset of tile indices

        EXAMPLES::

            sage: from slabbe import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.tile_positions([0])
            [(0, 0), (0, 2), (1, 1), (1, 3), (2, 0), (2, 2)]

        TESTS::

            sage: tiling.tile_positions([])
            []

        """
        set_M = set(M)
        L = []
        for i in range(self.width()):
            for j in range(self.height()):
                if self._table[i][j] in set_M:
                    L.append((i,j))
        return L

    @rename_keyword(fontsize='font')
    def tikz(self, color=None, font=r'\normalsize', rotate=None,
            id=True, id_color='', id_format='{}', label=True,
            label_shift=.2, label_color='black', scale=1, size=1,
            edges=True, draw_H=None, draw_V=None, extra_before='',
            extra_after=''):
        r"""
        Return a tikzpicture showing one solution.

        INPUT:

        - ``color`` -- None or dict from tile values -> tikz colors
        - ``font`` -- string (default: ``r'\normalsize'``
        - ``rotate`` -- list or ``None`` (default:``None``) list of four angles
          in degrees like ``(0,0,0,0)``, the rotation angle to apply to each
          label of Wang tiles. If ``None``, it performs a 90 degres rotation
          for left and right labels taking more than one character.
        - ``id`` -- boolean (default: ``True``), presence of the tile id
        - ``id_color`` -- string (default: ``''``) 
        - ``id_format`` -- string (default: ``r'{}'``) to be called with
          ``id_format.format(key)``
        - ``edges`` -- bool (default: ``True``) 
        - ``label`` -- boolean (default: ``True``) or a tuple of four
          boolean, presence of the color labels
        - ``label_shift`` -- number (default: ``.2``) translation distance
          of the label from the edge
        - ``label_color`` -- string (default: ``'black'``)
        - ``scale`` -- number or 2-tuple (default: ``1``), tikzpicture
          scale. If given a 2-tuple, it is interpreted as (xscale, yscale)
        - ``size`` -- number (default: ``1``) size of tiles
        - ``draw_H`` -- dict (default: ``None``) from tile values -> tikz
          draw commands. If ``None`` the values of the dict get replaced by
          straight lines, more precisely by ``r'\draw {} -- ++ (1,0);'``.
          Dict values must be strings ``s`` such that ``s.format((x,y))``
          works.
        - ``draw_V`` -- dict (default: ``None``) from tile values -> tikz
          draw commands. If ``None`` the values of the dict get replaced by
          straight lines, more precisely by ``r'\draw {} -- ++ (0,1);'``.
          Dict values must be strings ``s`` such that ``s.format((x,y))``
          works.
        - ``extra_before`` -- string (default: ``''``) extra lines of tikz
          code to add at the start
        - ``extra_after`` -- string (default: ``''``) extra lines of tikz
          code to add at the end

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
            \begin{tikzpicture}
            [scale=1]
            \tikzstyle{every node}=[font=\normalsize]
            % tile at position (x,y)=(0, 0)
            \node[] at (0.5, 0.5) {0};
            ...
            ... 97 lines not printed (3571 characters in total) ...
            ...
            \node[rotate=0,black] at (2.8, 3.5) {0};
            \node[rotate=0,black] at (2.5, 3.8) {0};
            \node[rotate=0,black] at (2.2, 3.5) {0};
            \node[rotate=0,black] at (2.5, 3.2) {0};
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

        ::

            sage: m = matrix(2,[1,1,0,1])
            sage: t = WangTiling(table, tiles, color).apply_matrix_transformation(m).tikz()

        Using puzzle boundary instead of colors::

            sage: tiles = [(0,3,1,2), (1,2,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: t = WangTiling(table, tiles)
            sage: draw_H = {0:r'\draw {} -- ++ (1/2,.2) -- ++ (1/2,-.2);',
            ....:           1:r'\draw {} -- ++ (1/2,.2) -- ++ (1/2,-.2);',
            ....:           2:r'\draw {} -- ++ (1/2,.2) -- ++ (1/2,-.2);',
            ....:           3:r'\draw {} -- ++ (1/2,.2) -- ++ (1/2,-.2);'}
            sage: v = r'\draw {} -- ++ (0,.4) -- ++ (.2,0) -- ++ (0,.2) -- ++ (-.2,0) -- ++ (0,.4);'
            sage: draw_V = {0:v, 1:v, 2:v, 3:v}
            sage: tikz = t.tikz(label=False, draw_H=draw_H, draw_V=draw_V)
        """
        if color is None:
            color = self._color
        lines = []
        lines.append(r'\begin{tikzpicture}')
        if isinstance(scale, tuple):
            xscale, yscale = scale
            lines.append('[xscale={},yscale={}]'.format(xscale,yscale))
        else:
            lines.append('[scale={}]'.format(scale))
        lines.append(r'\tikzstyle{{every node}}=[font={}]'.format(font))
        if extra_before:
            lines.append(extra_before)
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
                this_id = i if id else None
                position = (j,k)
                tile = self._tiles[i]
                right_edges = edges and (j == W - 1 or self._table[j+1][k] is None)
                top_edges = edges and (k == H - 1 or self._table[j][k+1] is None)
                more_lines = tile_to_tikz(tile, position, color=color,
                        id=this_id, id_color=id_color, id_format=id_format,
                        sizex=size, sizey=size, rotate=rotate, label=label,
                        label_shift=label_shift, label_color=label_color,
                        right_edges=right_edges, top_edges=top_edges,
                        left_edges=edges, bottom_edges=edges,
                        draw_H=draw_H, draw_V=draw_V)
                lines.extend(more_lines)
        if extra_after:
            lines.append(extra_after)
        lines.append(r'\end{tikzpicture}')
        from slabbe import TikzPicture
        return TikzPicture('\n'.join(lines))

