# -*- coding: utf-8 -*-
r"""
Wang tile solver

This uses MILP solvers like Coin or Gurobi. Coin can be installed with::

    sage -i cbc sagelib

EXAMPLES::

    sage: from slabbe.wang_tiles import WangTileSolver
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
from collections import Counter

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

    EXAMPLES::

        sage: from slabbe.wang_tiles import WangTileSolver
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
    def __init__(self, tiles, width, height, preassigned=None, color=None):
        r"""
        See class for documentation.

        EXAMPLES::

            sage: from slabbe.wang_tiles import WangTileSolver
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

            sage: from slabbe.wang_tiles import WangTileSolver
            sage: tiles = [(0,0,0,0), (1,1,1,1), (2,2,2,2)]
            sage: W = WangTileSolver(tiles,3,4)
            sage: p,x = W.milp()
            sage: p
            Mixed Integer Program  ( maximization, 36 variables, 29 constraints )
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
        for angle, D in enumerate(self._preassigned):
            for j,k in D:
                A = p.sum(tiles[i][angle]*x[i,j,k] for i in indices)
                name = "preassigned {} of {}".format(legend[angle], (j,k))
                p.add_constraint(A==D[(j,k)], name=name)

        p.set_objective(x[0,0,0])
        return p, x

    def solve(self, solver=None, solver_parameters=None):
        r"""
        Return a dictionary associating to each tile a list of positions
        where to find this tile.

        INPUT:

        - ``solver`` -- string or None (default: ``None``), other possible
          values are ``'Coin'`` or ``'Gurobi'``
        - ``solver_parameters`` -- dict (default: ``{}``), parameters given
          to the solver using method ``solver_parameter``. For a list of
          available parameters for example for the Gurobi backend, see
          dictionary ``parameters_type`` in the file
          ``sage/numerical/backends/gurobi_backend.pyx``

        OUTPUT:

            a wang tiling object

        EXAMPLES::

            sage: from slabbe.wang_tiles import WangTileSolver
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

        REFERENCES:

            How do I set solver_parameter to make Gurobi use more than one
            processor?, https://ask.sagemath.org/question/37726/
        """
        p,x = self.milp(solver=solver)
        if solver_parameters is None:
            solver_parameters = {}
        for key, value in solver_parameters.items():
            p.solver_parameter(key, value)
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

    .. NOTE::

        ``table[x][y]`` refers to the tile at position `(x,y)` using the
        cartesian coordinates. Thus, it is **not** using the matrix-like
        coordinates.

    EXAMPLES::

        sage: from slabbe.wang_tiles import WangTiling
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

            sage: from slabbe.wang_tiles import WangTileSolver
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

            sage: from slabbe.wang_tiles import WangTiling
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

            sage: from slabbe.wang_tiles import WangTiling
            sage: tiles = [(0,3,1,4), (1,4,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: tiling = WangTiling(table, tiles)
            sage: tiling.width()
            3
        """
        return len(self._table)

    def horizontal_words_list(self, side=3):
        r"""
        Return a list of horizontal words of colors appearing on a given
        side.

        INPUT

        - ``side`` -- integer in [0,1,2,3], 3 is for bottom

        EXAMPLES::

            sage: from slabbe.wang_tiles import WangTileSolver
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

    def number_of_occurences(self, pattern, avoid_border=0):
        r"""
        Return the number of occurences of the given pattern in the tiling.

        INPUT

        - ``pattern`` -- dict
        - ``avoid_border`` -- integer (default: 0), the size of the border
          to avoid during the computation

        EXAMPLES::

            sage: from slabbe.wang_tiles import WangTileSolver
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
            sage: tiling.tikz().pdf()   # not tested

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

            sage: from slabbe.wang_tiles import WangTiling
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

    def tikz(self, color=None, fontsize=r'\normalsize', rotate=(0,0,0,0),
            space=.2, scale=1, transformation_matrix=None):
        r"""
        Return a tikzpicture showing one solution.

        INPUT:

        - ``color`` -- None or dict from tile values -> tikz colors
        - ``fontsize`` -- string (default: ``r'\normalsize'``
        - ``rotate`` -- list (default:``(0,0,0,0)``) of four angles in
          degrees, the rotation angle to apply to each label of Wang tiles
        - ``space`` -- number (default: ``.2``) translation distance of the
          label from the edge
        - ``scale`` -- number (default: ``1``), tikzpicture scale
        - ``transformation_matrix`` -- matrix (default: ``None``), a matrix
          to apply to the coordinate before drawing, it can be in
          ``SL(2,ZZ)`` or not.

        .. TODO::

            - Fix the top and right lines when using the transformation
              matrix.

        EXAMPLES::

            sage: from slabbe.wang_tiles import WangTileSolver
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
            sage: tiling = W.solve('Gurobi')
            sage: color = {0:'white',1:'red',2:'blue',3:'green'}
            sage: t = tiling.tikz(color=color)

        Using some blank tiles::

            sage: tiles = [(0,3,1,2), (1,2,0,3)]
            sage: table = [[0, 1, None, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: color = {0:'white',1:'red',2:'blue',3:'green'}
            sage: tiling = WangTiling(table, tiles, color)
            sage: t = tiling.tikz()

        Testing the options::

            sage: tiles = [(0,3,1,2), (1,2,0,3)]
            sage: table = [[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1]]
            sage: color = {0:'white',1:'red',2:'blue',3:'green'}
            sage: t = WangTiling(table, tiles, color).tikz(fontsize=r'\Huge')
            sage: t = WangTiling(table, tiles, color).tikz(rotate=(0,90,0,0))
            sage: t = WangTiling(table, tiles, color).tikz(space=.05)
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
        s = space # because it is shorter to write below
        lines = []
        lines.append(r'\begin{{tikzpicture}}[scale={}]'.format(scale))
        W = self.width()
        H = self.height()
        # missing lines on the top
        for j in range(W):
            lines.append(r'\draw {} -- {};'.format((j,H), (j+1,H)))
        # missing lines on the right
        for k in range(H):
            lines.append(r'\draw {} -- {};'.format((W,k), (W,k+1)))
        # the tiles with borders left and below
        for j in range(W):
            for k in range(H):
                i = self._table[j][k]
                if i is None:
                    # this is a blank tile
                    continue
                x,y = transformation_matrix*vector((j,k))
                tile = self._tiles[i]
                lines.append('% tile at position (j,k)={} or (x,y)={}'.format((j,k), (x,y)))
                if color:
                    triangle = r'\fill[{}] {} -- {} -- {};'
                    c = (x+.5,y+.5)
                    lines.append(triangle.format(color[tile[0]],(x+1,y),c,(x+1,y+1)))
                    lines.append(triangle.format(color[tile[1]],(x,y+1),c,(x+1,y+1)))
                    lines.append(triangle.format(color[tile[2]],(x,y),c,(x,y+1)))
                    lines.append(triangle.format(color[tile[3]],(x,y),c,(x+1,y)))
                lines.append(r'\draw {} -- {};'.format((x,y), (x+1,y)))
                lines.append(r'\draw {} -- {};'.format((x,y), (x,y+1)))
                node_str = r'\node[rotate={},font={}] at {} {{{}}};'
                lines.append(node_str.format(rotate[0],fontsize,(x+1-s,y+.5),  tile[0]))
                lines.append(node_str.format(rotate[1],fontsize,(x+.5, y+1-s), tile[1]))
                lines.append(node_str.format(rotate[2],fontsize,(x+s,  y+.5),  tile[2]))
                lines.append(node_str.format(rotate[3],fontsize,(x+.5, y+s),   tile[3]))
        lines.append(r'\end{tikzpicture}')
        from slabbe import TikzPicture
        return TikzPicture('\n'.join(lines))

