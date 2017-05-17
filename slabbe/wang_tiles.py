# -*- coding: utf-8 -*-
r"""
Wang tile solver

This uses Coin solver which can be installed with::

    sage -i cbc sagelib

EXAMPLES::

    sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: _ = W.tikz().pdf()

::

    sage: tiles = [(1/2,1/2,1/2,1/2), (1,1,1,1), (2,2,2,2)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: _ = W.tikz().pdf()

::

    sage: tiles = [(0,3,1,4), (1,4,0,3)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: table = W.solve()
    sage: table
    [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]

Kari-Culik (Here 0' is replaced by 10)::

    sage: divide_by_2 = [(10,0,10,10), (10,1,10,2), (1/2,0,10,1), (1/2,10,10,1),
    ....:     (1/2,0,1/2,10), (1/2,1,1/2,2), (10,1,1/2,1)]
    sage: times_3 = [(1,2,0,1), (2,1,0,1), (2,2,1,1), (0,1,1,0), (0,2,2,0),
    ....:     (1,1,2,0)]
    sage: tiles = divide_by_2 + times_3
    sage: W = WangTileSolver(tiles,3,4)
    sage: _ = W.tikz().pdf()

Rao-Jeandel::

    sage: tiles = [(2,4,2,1), (2,2,2,0), (1,1,3,1), (1,2,3,2), (3,1,3,3),
    ....: (0,1,3,1), (0,0,0,1), (3,1,0,2), (0,2,1,2), (1,2,1,4), (3,3,1,2)]
    sage: W = WangTileSolver(tiles,3,4)
    sage: _ = W.tikz().pdf()

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
from sage.numerical.mip import MixedIntegerLinearProgram

class WangTileSolver(object):
    r"""
    INPUT:

    - ``tiles`` -- list of tiles, a tile is a 4-tuple (right color, top
        color, left color, bottom color)
    - ``width`` -- integer
    - ``height`` -- integer
    - ``precolor`` -- None or list of 4 dict or the form ``[{}, {}, {}, {}]``
        right, top, left, bottom colors preassigned to some positions (on
        the border or inside)

    EXAMPLES::

        sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
        sage: W = WangTileSolver(tiles, 3, 3)
	sage: W.solve()
	[[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    With color 2 preassigned to the right part of tile at position (1,1)::

	sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
	sage: right = {(1,1):2}
	sage: W = WangTileSolver(tiles,3,3,precolor=[right,{},{},{}])
	sage: W.solve()
	[[2, 2, 2], [2, 2, 2], [2, 2, 2]]

    When constraints are inconsistent::

	sage: right = {(1,1):1, (2,2):0}
	sage: W = WangTileSolver(tiles,3,3,precolor=[right,{},{},{}])
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
    def __init__(self, tiles, width, height, precolor=None):
        r"""
        See class for documentation.

        EXAMPLES::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles, 3, 4)
        """
        self._tiles = tiles
        self._width = width
        self._height = height
        if precolor is None:
            precolor = [{}, {}, {}, {}]
        self._precolor = precolor

    def milp(self, solver='Coin'):
        r"""
        INPUT:

        - ``solver`` -- None or string

        EXAMPLES::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: p,x = W.milp()
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
        for angle, D in enumerate(self._precolor):
            for j,k in D:
                A = p.sum(tiles[i][angle]*x[i,j,k] for i in indices)
                name = "preassigned {} of {}".format(legend[angle], (j,k))
                p.add_constraint(A==D[(j,k)], name=name)

        p.set_objective(x[0,0,0])
        return p, x

    def solve(self, solver='Coin'):
        r"""
        Return a dictionary associating to each tile a list of positions
        where to find this tile.

        INPUT:

        - ``solver`` --

        EXAMPLES::

            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: table = W.solve()
            sage: table
            [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]

        The tile at position (1,3) is::

            sage: table[1][3]
            0
        """
        p,x = self.milp(solver)
        p.solve()
        soln = p.get_values(x)
        support = [key for key in soln if soln[key]]
        assert len(support) == self._width * self._height, "yoo"
        table = [[0]*self._height for _ in range(self._width)]
        for i,j,k in support:
            table[j][k] = i
        return table

    def tikz(self, solver='Coin', color=None):
        r"""
        Return a tikzpicture showing one solution.

        INPUT:

        - ``solver`` --
        - ``color`` -- None or dict

        EXAMPLES::

            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: t = W.tikz()
            \documentclass[tikz]{standalone}
            \usepackage{amsmath}
            \begin{document}
            \begin{tikzpicture}
            % tile at position (0,0)
            \draw (0, 0) -- (1, 0);
            \draw (0, 0) -- (0, 1);
            \node[left]  at (1, 0.5) {0};
            ...
            ... 76 lines not printed (2353 characters in total) ...
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
            sage: W = WangTileSolver(tiles,3,4)
            sage: color = {0:'white',1:'red',2:'blue',3:'green'}
            sage: t = W.tikz(color=color)
        """
        table = self.solve(solver)
        lines = []
        lines.append(r'\begin{tikzpicture}')
        for j in range(self._width):
            for k in range(self._height):
                i = table[j][k]
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






