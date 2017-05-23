# -*- coding: utf-8 -*-
r"""
Wang tile solver

This uses Coin solver which can be installed with::

    sage -i cbc sagelib

EXAMPLES::

    sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: tiling = W.solve()
    sage: _ = tiling.tikz().pdf()

::

    sage: tiles = [(1/2,1/2,1/2,1/2), (1,1,1,1), (2,2,2,2)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: tiling = W.solve()
    sage: _ = tiling.tikz().pdf()

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
    sage: _ = tiling.tikz().pdf()

Rao-Jeandel::

    sage: tiles = [(2,4,2,1), (2,2,2,0), (1,1,3,1), (1,2,3,2), (3,1,3,3),
    ....: (0,1,3,1), (0,0,0,1), (3,1,0,2), (0,2,1,2), (1,2,1,4), (3,3,1,2)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: tiling = W.solve()
    sage: _ = tiling.tikz().pdf()

.. TODO::

 - Reduce the problem to SAT
                                         
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
from sage.misc.cachefunc import cached_method
from sage.numerical.mip import MixedIntegerLinearProgram

class WangTileSolver(object):
    r"""
    INPUT:

    - ``tiles`` -- list of tiles, a tile is a 4-tuple (right color, top
        color, left color, bottom color)
    - ``width`` -- integer
    - ``height`` -- integer
    - ``preassigned`` -- None or list of 4 dict or the form ``[{}, {}, {}, {}]``
        right, top, left, bottom colors preassigned to some positions (on
        the border or inside)
    - ``solver`` -- None or string

    EXAMPLES::

        sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
        sage: W = WangTileSolver(tiles, 3, 3)
        sage: tiling = W.solve()
        sage: tiling._table
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    With color 2 preassigned to the right part of tile at position (1,1)::

        sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
        sage: right = {(1,1):2}
        sage: W = WangTileSolver(tiles,3,3,preassigned=[right,{},{},{}])
        sage: tiling = W.solve()
        sage: tiling._table
        [[2, 2, 2], [2, 2, 2], [2, 2, 2]]

    When constraints are inconsistent::

        sage: right = {(1,1):1, (2,2):0}
        sage: W = WangTileSolver(tiles,3,3,preassigned=[right,{},{},{}])
        sage: W.solve()
        Traceback (most recent call last):
        ...
        MIPSolverException: CBC : The problem or its dual has been proven infeasible!

    TESTS:

    Colors must be convertable to float::

        sage: tiles = [('a','a','a','a'), ('b','b','b','b')]
        sage: W = WangTileSolver(tiles,3,4)
        sage: _ = W.tikz().pdf()
        Traceback (most recent call last):
        ...
        ValueError: could not convert string to float: a
    """
    def __init__(self, tiles, width, height, preassigned=None, color=None, solver='Coin'):
        r"""
        See class for documentation.

        EXAMPLES::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles, 3, 4)
        """
        self._tiles = tiles
        self._width = width
        self._height = height
        if preassigned is None:
            preassigned = [{}, {}, {}, {}]
        self._preassigned = preassigned
        self._color = color
        self._solver = solver

    def milp(self):
        r"""
        EXAMPLES::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: p,x = W.milp()
        """
        tiles = self._tiles
        indices = range(len(tiles))

        # x[i,j,k] == 1 iff tile i is at position (j,k)
        p = MixedIntegerLinearProgram(solver=self._solver)
        x = p.new_variable(binary=True)

        # exactly one tile at each position (j,k)
        for j in range(self._width):
            for k in range(self._height):
                S = p.sum(x[i,j,k] for i in indices)
                name = "one tile at {}".format((j,k))
                p.add_constraint(S==1, name=name)

        # matching vertical colors
        for j in range(self._width-1):
            for k in range(self._height):
                A = p.sum(tiles[i][0]*x[i,j,k] for i in indices)
                B = p.sum(tiles[i][2]*x[i,j+1,k] for i in indices)
                name = "matching right of {}".format((j,k))
                p.add_constraint(A==B, name=name)

        # matching horizontal colors
        for j in range(self._width):
            for k in range(self._height-1):
                A = p.sum(tiles[i][1]*x[i,j,k] for i in indices)
                B = p.sum(tiles[i][3]*x[i,j,k+1] for i in indices)
                name = "matching top of {}".format((j,k))
                p.add_constraint(A==B, name=name)

        # matching preassigned color constraints
        legend = {0:'right',1:'top',2:'left',3:'bottom'}
        for angle, D in enumerate(self._preassigned):
            for j,k in D:
                A = p.sum(tiles[i][angle]*x[i,j,k] for i in indices)
                name = "preassigned {} of {}".format(legend[angle], (j,k))
                p.add_constraint(A==D[(j,k)], name=name)

        p.set_objective(x[0,0,0])
        return p, x

    @cached_method
    def solve(self):
        r"""
        Return a dictionary associating to each tile a list of positions
        where to find this tile.

        EXAMPLES::

            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: table = tiling._table
            sage: table
            [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]

        The tile at position (1,3) is::

            sage: table[1][3]
            0
        """
        p,x = self.milp()
        p.solve()
        soln = p.get_values(x)
        support = [key for key in soln if soln[key]]
        assert len(support) == self._width * self._height, "yoo"
        table = [[0]*self._height for _ in range(self._width)]
        for i,j,k in support:
            table[j][k] = i
        return WangTiling(table, self._tiles, color=self._color)

class WangTiling(object):
    r"""
    INPUT:

    - ``table`` -- list of lists
    - ``tiles`` -- list of tiles, a tile is a 4-tuple (right color, top
        color, left color, bottom color)
    - ``color`` -- dict (default: None)

    EXAMPLES::

        sage: tiles = [(0,3,1,4), (1,4,0,3)]
        sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
        sage: tiling = WangTiling(table, tiles)
        sage: tiling
        A wang tiling of a 3 x 4 rectangle
    """
    def __init__(self, table, tiles, color=None):
        r"""
        See class for documentation.

        EXAMPLES::

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

            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.height()
            4
        """
        return len(table[0])

    def width(self):
        r"""
        EXAMPLES::

            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.width()
            3
        """
        return len(table)

    def horizontal_words_list(self, side=3):
        r"""
        Return a list of horizontal words of colors appearing on a given
        side.

        INPUT

        - ``side`` -- integer in [0,1,2,3], 3 is for bottom

        EXAMPLES::

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

    def tikz(self, color=None):
        r"""
        Return a tikzpicture showing one solution.

        INPUT:

        - ``color`` -- None or dict

        EXAMPLES::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: t = tiling.tikz()
            sage: t
            \documentclass[tikz]{standalone}
            \usepackage{amsmath}
            \begin{document}
            \begin{tikzpicture}
            % tile at position (0,0)
            \draw (0, 0) -- (1, 0);
            \draw (0, 0) -- (0, 1);
            \node[left]  at (1, 0.5) {0};
            ...
            ... 76 lines not printed (2365 characters in total) ...
            ...
            \node[left]  at (3, 3.5) {0};
            \node[below] at (2.5, 4) {0};
            \node[right] at (2, 3.5) {0};
            \node[above] at (2.5, 3) {0};
            \end{tikzpicture}
            \end{document}
            sage: _ = t.pdf()

        With colors::

            sage: tiles = [(0,2,1,3), (1,3,0,2)]
            sage: color = {0:'white',1:'red',2:'blue',3:'green'}
            sage: W = WangTileSolver(tiles,3,4,color=color)
            sage: tiling = W.solve()
            sage: t = tiling.tikz()

        With colors, alternatively::

            sage: tiles = [(0,2,1,3), (1,3,0,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: tiling = W.solve()
            sage: color = {0:'white',1:'red',2:'blue',3:'green'}
            sage: t = tiling.tikz(color=color)
        """
        if color is None:
            color = self._color
        lines = []
        lines.append(r'\begin{tikzpicture}')
        for j in range(self.width()):
            for k in range(self.height()):
                i = self._table[j][k]
                tile = self._tiles[i]
                lines.append('% tile at position {}'.format((j,k)))
                if color:
                    tri = r'\fill[{}] {} -- {} -- {};'
                    c = (j+.5,k+.5)
                    lines.append(tri.format(color[tile[0]],(j+1,k),c,(j+1,k+1)))
                    lines.append(tri.format(color[tile[1]],(j,k+1),c,(j+1,k+1)))
                    lines.append(tri.format(color[tile[2]],(j,k),c,(j,k+1)))
                    lines.append(tri.format(color[tile[3]],(j,k),c,(j+1,k)))
                lines.append(r'\draw {} -- {};'.format((j,k), (j+1,k)))
                lines.append(r'\draw {} -- {};'.format((j,k), (j,k+1)))
                lines.append(r'\node[left]  at {} {{{}}};'.format((j+1,k+.5), tile[0]))
                lines.append(r'\node[below] at {} {{{}}};'.format((j+.5,k+1), tile[1]))
                lines.append(r'\node[right] at {} {{{}}};'.format((j,k+.5), tile[2]))
                lines.append(r'\node[above] at {} {{{}}};'.format((j+.5,k), tile[3]))
        lines.append(r'\end{tikzpicture}')
        from slabbe import TikzPicture
        return TikzPicture('\n'.join(lines))

