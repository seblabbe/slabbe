# -*- coding: utf-8 -*-
r"""
2d substitutions

EXAMPLES::

    sage: from slabbe import Substitution2d
    sage: A = [[0,1],[2,3]]
    sage: B = [[4,5]]
    sage: d = {0:A, 1:B}
    sage: Substitution2d(d)
    Substitution 2d: {0: [[0, 1], [2, 3]], 1: [[4, 5]]}
"""
#*****************************************************************************
#       Copyright (C) 2017 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
import itertools

from sage.misc.decorators import rename_keyword

class Substitution2d(object):
    r"""
    INPUT:

    - ``d`` -- dict, key -> value, where each value is a table such that
      table[x][y] refers to the tile at position (x,y) in cartesian
      coordinates (*not* in the matrix-like coordinates)

    EXAMPLES::

        sage: from slabbe import Substitution2d
        sage: A = [[0,1],[2,3]]
        sage: B = [[4,5]]
        sage: d = {0:A, 1:B}
        sage: s = Substitution2d(d)
        sage: s
        Substitution 2d: {0: [[0, 1], [2, 3]], 1: [[4, 5]]}
    """
    def __init__(self, d):
        r"""
        See class documentation.

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[2,3]]
            sage: B = [[4,5]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
        """
        self._d = d

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[2,3]]
            sage: B = [[4,5]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: s
            Substitution 2d: {0: [[0, 1], [2, 3]], 1: [[4, 5]]}
        """
        return "Substitution 2d: {}".format(self._d)

    def _latex_(self, ncolumns=8, align='l', variableA=None,
            variableB=None):
        r"""
        INPUT:

        - ``ncolumns`` -- integer
        - ``align`` -- character (default:``'l'``), latex alignment symbol
          ``'l'``, ``'r'`` or ``'c'``.
        - ``variableA`` -- string or ``None``
        - ``variableB`` -- string or ``None``

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[2,3]]
            sage: B = [[4,5]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: latex(s)
            \begin{array}{llllllll}
            0\mapsto \left(\begin{array}{rr}
            1 & 3 \\
            0 & 2
            \end{array}\right)
            ,&
            1\mapsto \left(\begin{array}{r}
            5 \\
            4
            \end{array}\right)
            .
            \end{array}

        ::

            sage: s._latex_(2, 'c', variableA='x', variableB='y')
            \begin{array}{cc}
            x_{0}\mapsto \left(\begin{array}{rr}
            y_{1} & y_{3} \\
            y_{0} & y_{2}
            \end{array}\right)
            ,&
            x_{1}\mapsto \left(\begin{array}{r}
            y_{5} \\
            y_{4}
            \end{array}\right)
            .
            \end{array}

        """
        from sage.matrix.constructor import matrix
        from sage.misc.latex import latex
        from sage.misc.latex import LatexExpr
        from sage.calculus.var import var
        from slabbe.matrices import map_coefficients_to_variable_index
        lines = []
        lines.append(r'\begin{{array}}{{{}}}'.format(align*ncolumns))
        for i,(key,table) in enumerate(self._d.items()):
            M = matrix.column([col[::-1] for col in table])
            if variableA:
                key = var('{}_{}'.format(variableA, key))
            if variableB:
                M = map_coefficients_to_variable_index(M, variableB)
            lines.append(r'{}\mapsto {}'.format(latex(key), latex(M)))
            if (i+1) == len(self._d):
                lines.append(r'.')
            elif (i+1) % ncolumns == 0:
                lines.append(r',\\')
            else:
                lines.append(r',&')
        lines.append(r'\end{array}')
        return LatexExpr('\n'.join(lines))

    @classmethod
    def from_1d_row_substitution(self, s):
        r"""
        INPUT:

        - ``s`` -- dict

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: fibo = {0:[0,1], 1:[0]}
            sage: s = Substitution2d.from_1d_row_substitution(fibo)
            sage: s
            Substitution 2d: {0: [[0], [1]], 1: [[0]]}
        """
        d = {key:[[a] for a in row] for key,row in s.items()}
        return Substitution2d(d)

    @classmethod
    def from_1d_column_substitution(self, s):
        r"""
        INPUT:

        - ``s`` -- dict

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: fibo = {0:[0,1], 1:[0]}
            sage: s = Substitution2d.from_1d_column_substitution(fibo)
            sage: s
            Substitution 2d: {0: [[0, 1]], 1: [[0]]}
        """
        d = {key:[col] for key,col in s.items()}
        return Substitution2d(d)

    @classmethod
    def from_permutation(self, d):
        r"""
        INPUT:

        - ``d`` -- dict

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: s = Substitution2d.from_permutation({4:0, 5:1})
            sage: s
            Substitution 2d: {4: [[0]], 5: [[1]]}

        ::

            sage: A = [[5,6],[7,8]]
            sage: B = [[6,5],[9,8]]
            sage: t = Substitution2d({0:A, 1:B})
            sage: t
            Substitution 2d: {0: [[5, 6], [7, 8]], 1: [[6, 5], [9, 8]]}
            sage: t*s
            Substitution 2d: {4: [[5, 6], [7, 8]], 5: [[6, 5], [9, 8]]}

        ::

            sage: u = Substitution2d.from_permutation({5:0, 6:1, 7:2, 8:3, 9:4})
            sage: u
            Substitution 2d: {8: [[3]], 9: [[4]], 5: [[0]], 6: [[1]], 7: [[2]]}
            sage: u * t
            Substitution 2d: {0: [[0, 1], [2, 3]], 1: [[1, 0], [4, 3]]}
        """
        d = {key:[[val]] for key,val in d.items()}
        return Substitution2d(d)

    def call_on_row(self, row):
        r"""
        INPUT:

        - ``row`` -- list

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[2,3]]
            sage: B = [[4,5]]
            sage: C = [[6,7,8]]
            sage: d = {0:A, 1:B, 2:C}
            sage: s = Substitution2d(d)
            sage: row = [0,1,1,0]
            sage: s.call_on_row(row)
            [[0, 1], [2, 3], [4, 5], [4, 5], [0, 1], [2, 3]]
            sage: s.call_on_row([2])
            [[6, 7, 8]]

        TESTS::

            sage: s.call_on_row([])
            []
            sage: s.call_on_row([1,2])
            Traceback (most recent call last):
            ...
            ValueError: the image of the row contains columns of different height (=set([2, 3]))
        """
        if not row:
            return []
        columns = []
        for r in row:
            for column in self._d[r]:
                columns.append(column)
        s = set(map(len,columns))
        if len(s) != 1:
            raise ValueError("the image of the row contains columns of "
                    "different height (={})".format(s))
        return columns

    def call_on_column(self, column):
        r"""
        INPUT:

        - ``column`` -- list

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[2,3]]
            sage: B = [[4],[5]]
            sage: C = [[6,7,8]]
            sage: d = {0:A, 1:B, 2:C}
            sage: s = Substitution2d(d)
            sage: s.call_on_column([0])
            [[0, 1], [2, 3]]
            sage: s.call_on_column([0,1])
            [[0, 1, 4], [2, 3, 5]]
            sage: s.call_on_column([0,1,1,0,0])
            [[0, 1, 4, 4, 0, 1, 0, 1], [2, 3, 5, 5, 2, 3, 2, 3]]

        TESTS::

            sage: s.call_on_column([])
            []
            sage: s.call_on_column([0,2])
            Traceback (most recent call last):
            ...
            ValueError: the image of 2 in the column (=[0, 2]) has width 1
            but the image of another has width 2

        """
        if not column:
            return []
        width = len(self._d[column[0]])
        rep = [[] for _ in range(width)]
        for a in column:
            image_a = self._d[a]
            if not len(image_a) == width:
                raise ValueError("the image of {} in the column (={}) has"
                        " width {} but the image of another"
                        " has width {}".format(a, column, len(image_a), width))
            for i,col in enumerate(image_a):
                rep[i].extend(col)
        return rep

    def __call__(self, table, order=1):
        r"""
        INPUT:

        - ``table`` -- list of list, such that table[x][y] refers to the
          tile at position (x,y) in cartesian coordinates (*not* in the
          matrix-like coordinates)
        -  ``order`` - integer or plus ``Infinity`` (default: 1)

        TODO: implement another call on table with matrix like coordinates
        
        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[2,3]]
            sage: B = [[4,5],[6,7]]
            sage: C = [[8,9]]
            sage: d = {0:A, 1:B, 2:C}
            sage: s = Substitution2d(d)
            sage: table = [[0,1],[1,1]]
            sage: s(table)
            [[0, 1, 4, 5], [2, 3, 6, 7], [4, 5, 4, 5], [6, 7, 6, 7]]

        ::

            sage: A = [[0,1],[2,0]]
            sage: B = [[2,1],[2,0]]
            sage: C = [[1,2],[1,1]]
            sage: d = {0:A, 1:B, 2:C}
            sage: s = Substitution2d(d)
            sage: s([[0]])
            [[0, 1], [2, 0]]
            sage: s([[0]],2)
            [[0, 1, 2, 1], [2, 0, 2, 0], [1, 2, 0, 1], [1, 1, 2, 0]]
            sage: s([[0]],3)
            [[0, 1, 2, 1, 1, 2, 2, 1],
             [2, 0, 2, 0, 1, 1, 2, 0],
             [1, 2, 0, 1, 1, 2, 0, 1],
             [1, 1, 2, 0, 1, 1, 2, 0],
             [2, 1, 1, 2, 0, 1, 2, 1],
             [2, 0, 1, 1, 2, 0, 2, 0],
             [2, 1, 2, 1, 1, 2, 0, 1],
             [2, 0, 2, 0, 1, 1, 2, 0]]

        TESTS::

            sage: A = [[0,1],[2,3]]
            sage: B = [[4,5],[6,7]]
            sage: C = [[8,9]]
            sage: d = {0:A, 1:B, 2:C}
            sage: s = Substitution2d(d)
            sage: s([])
            []
            sage: s([[0,1], [1,2]])
            Traceback (most recent call last):
            ...
            ValueError: the image of 2 in the column (=[1, 2]) has width 1
            but the image of another has width 2
        """
        if order == 1:
            columns = []
            for col in table:
                col_image = self.call_on_column(col)
                columns.extend(col_image)
            return columns
        elif order > 1:
            return self(self(table, order-1))
        elif order == 0:
            return table
        else:
            raise TypeError("order (%s) must be a positive integer" % order)

    def __mul__(self, other):
        r"""
        INPUT:

        - ``other`` -- substitution 2d

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[0,1]]
            sage: B = [[1,0],[1,1]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: s
            Substitution 2d: {0: [[0, 1], [0, 1]], 1: [[1, 0], [1, 1]]}
            sage: s * s
            Substitution 2d: {0: [[0, 1, 1, 0], [0, 1, 1, 1], [0, 1, 1, 0], [0, 1, 1, 1]], 1: [[1, 0, 0, 1], [1, 1, 0, 1], [1, 0, 1, 0], [1, 1, 1, 1]]}
            sage: s * s * s
            Substitution 2d: {0: [[0, 1, 1, 0, 1, 0, 0, 1], [0, 1, 1, 1, 1, 1, 0, 1], [0, 1, 1, 0, 1, 0, 1, 0], [0, 1, 1, 1, 1, 1, 1, 1], [0, 1, 1, 0, 1, 0, 0, 1], [0, 1, 1, 1, 1, 1, 0, 1], [0, 1, 1, 0, 1, 0, 1, 0], [0, 1, 1, 1, 1, 1, 1, 1]], 1: [[1, 0, 0, 1, 0, 1, 1, 0], [1, 1, 0, 1, 0, 1, 1, 1], [1, 0, 1, 0, 0, 1, 1, 0], [1, 1, 1, 1, 0, 1, 1, 1], [1, 0, 0, 1, 1, 0, 0, 1], [1, 1, 0, 1, 1, 1, 0, 1], [1, 0, 1, 0, 1, 0, 1, 0], [1, 1, 1, 1, 1, 1, 1, 1]]}

        ::

            sage: t = Substitution2d.from_permutation({0:4, 1:5})
            sage: t * s
            Substitution 2d: {0: [[4, 5], [4, 5]], 1: [[5, 4], [5, 5]]}
            sage: s * t
            Traceback (most recent call last):
            ...
            ValueError: codomain alphabet of other (=set([4, 5])) must be
            included in domain alphabet of self (=set([0, 1]))
        """
        if not self.domain_alphabet() >= other.codomain_alphabet():
            raise ValueError("codomain alphabet of other (={}) must be included"
                " in domain alphabet of self (={})".format(other.codomain_alphabet(),
                                                            self.domain_alphabet()))
        d = {}
        for a,image_a in other._d.items():
            d[a] = self(image_a)
        return Substitution2d(d)


    def domain_alphabet(self):
        r"""
        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[5,6],[7,8]]
            sage: B = [[6,5],[9,8]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: s.domain_alphabet()
            {0, 1}
        """
        return set(self._d.keys())

    def codomain_alphabet(self):
        r"""
        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[5,6],[7,8]]
            sage: B = [[6,5],[9,8]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: s.codomain_alphabet()
            {5, 6, 7, 8, 9}
        """
        s = set()
        for a,image_a in self._d.items():
            for column in image_a:
                s.update(column)
        return s

    def incidence_matrix(self):
        r"""
        Return the incidence matrix of self.

        Some default ordering (sorted) is used for the domain and codomain
        alphabet.

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[2,3]]
            sage: B = [[4,5]]
            sage: C = [[6,7,8]]
            sage: d = {0:A, 1:B, 2:C}
            sage: s = Substitution2d(d)
            sage: s.incidence_matrix()
            [1 0 0]
            [1 0 0]
            [1 0 0]
            [1 0 0]
            [0 1 0]
            [0 1 0]
            [0 0 1]
            [0 0 1]
            [0 0 1]
        """
        from collections import Counter
        from sage.matrix.constructor import matrix
        domain_alphabet = sorted(self.domain_alphabet())
        codomain_alphabet = sorted(self.codomain_alphabet())
        L = []
        for a in domain_alphabet:
            c = Counter()
            image_a = self._d[a]
            for column in image_a:
                c.update(column)
            v = [c[b] for b in codomain_alphabet]
            L.append(v)
        return matrix.column(L)

    def desubstitute(self, tiles, function=None):
        r"""
        Return the Wang tile set obtained from the desubstitution of the
        given Wang tile set.

        INPUT:

        - ``tiles`` -- list of Wang tiles, each tile being a 4-tuple of
          (east, north, west, south) colors
        - ``fn`` -- a function (default: ``None``) to apply to the
          new colors which are tuple of previous colors

        OUTPUT:

            dict, key -> tile

        """
        raise NotImplementedError('method desubstitute moved to wang tile')

    @rename_keyword(fontsize='font')
    def wang_tikz(self, domain_tiles, codomain_tiles, color=None, size=1,
            scale=1, font=r'\normalsize', rotate=None, label_shift=.2,
            transformation_matrix=None, ncolumns=4, tabular='tabular',
            align='l', direction='right'):
        r"""
        Return the tikz code showing what the substitution A->B* does on
        Wang tiles.

        INPUT:

        - ``domain_tiles`` -- tiles of the domain
        - ``codomain_tiles`` -- tiles of the codomain
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
        - ``direction`` -- string (default: ``'right'``) or ``'down'``

        OUTPUT:

            dict, key -> tile

        EXAMPLES::

            sage: from slabbe import WangTileSet, Substitution2d
            sage: A = [[0,1,2],[1,0,0]]
            sage: B = [[0,1,2]]
            sage: d = {4:A, 5:B}
            sage: s = Substitution2d(d)
            sage: codomain_tiles = [(0,3,1,4), (1,4,0,3), (5,6,7,8)]
            sage: W = WangTileSet(codomain_tiles)
            sage: fn = lambda colors:''.join(map(str, colors))
            sage: domain_tiles = W.desubstitute(s, fn)
            sage: output = s.wang_tikz(domain_tiles, codomain_tiles, rotate=(90,0,90,0))
            sage: view(output)    # not tested

        Applying a transformation matrix::

            sage: M = matrix(2, [1,1,0,1])
            sage: output = s.wang_tikz(domain_tiles, codomain_tiles, 
            ....:                    transformation_matrix=M)
            sage: view(output)    # not tested

        Down direction::

            sage: output = s.wang_tikz(domain_tiles, codomain_tiles,
            ....:                      direction='down')
            sage: view(output)    # not tested
        """
        from slabbe.wang_tiles import tile_to_tikz, WangTileSet, WangTiling
        if isinstance(codomain_tiles, WangTileSet):
            codomain_tiles = codomain_tiles._tiles
        if isinstance(domain_tiles, WangTileSet):
            domain_tiles = domain_tiles._tiles

        for a in self.codomain_alphabet():
            try:
                codomain_tiles[a]
            except IndexError:
                raise ValueError("codomain_alphabet={}, but tiles are"
                        " {}".format(self.codomain_alphabet(),
                        codomain_tiles))
        for a in self.domain_alphabet():
            try:
                domain_tiles[a]
            except IndexError:
                raise ValueError("domain_alphabet={}, but tiles are"
                        " {}".format(self.domain_alphabet(),
                        domain_tiles))

        if not direction in ['right', 'down']:
            raise ValueError("direction(={}) must be 'right' or"
                    " 'down'".format(direction))

        lines = []
        lines.append(r'\begin{{{}}}{{{}}}'.format(tabular, align*ncolumns))
        for i,a in enumerate(self._d):

            lines.append(r'\begin{tikzpicture}')
            lines.append(r'[scale={}]'.format(scale))
            lines.append(r'\tikzstyle{{every node}}=[font={}]'.format(font))

            desubstituted_tile = domain_tiles[a] 
            new_lines = tile_to_tikz(desubstituted_tile, (0,0), color=color,
                    size=size, rotate=rotate, label_shift=label_shift,
                    top_right_edges=True)
            lines.extend(new_lines)

            if direction == 'right':
                lines.append(r'\node at (1.5,.5) {$\mapsto$};')
            elif direction == 'down':
                lines.append(r'\node[rotate=-90] at (.5,-.5) {$\mapsto$};')

            image_a = self._d[a]
            tiling = WangTiling(image_a, codomain_tiles, color)
            tikz = tiling.tikz(color=color, font=font, rotate=rotate,
                    label_shift=label_shift, scale=scale,
                    transformation_matrix=transformation_matrix)

            if direction == 'right':
                xshift = 2.0 + .5 * len(image_a)
                yshift = .5
            elif direction == 'down':
                xshift = .5
                yshift = -1.0 - .5 * len(image_a[0])
            lines.append(r'\node at ({},{}) {{{}}};'.format(xshift, yshift,
                                             tikz.tikz_picture_code()))

            lines.append(r'\end{tikzpicture}')
            if (i+1) == len(self._d):
                lines.append(r'.')
            elif (i+1) % ncolumns == 0:
                lines.append(r',\\')
            else:
                lines.append(r',&')
        lines.append(r'\end{{{}}}'.format(tabular))
        from sage.misc.latex import LatexExpr
        return LatexExpr('\n'.join(lines))

    def list_2x2_factors(self):
        r"""
        Return the list of 2x2 factors in the associated substitutive shift.

        INPUT:

        - ``self`` -- expansive and primitive 2d substitution

        OUTPUT:

            list of tables

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[0,1]]
            sage: B = [[1,0],[1,1]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: s.list_2x2_factors()
            [[[0, 1], [0, 1]],
             [[1, 0], [1, 1]],
             [[1, 1], [1, 0]],
             [[1, 1], [1, 1]],
             [[1, 1], [0, 1]],
             [[1, 1], [0, 0]],
             [[0, 1], [1, 1]],
             [[1, 0], [0, 1]],
             [[0, 0], [1, 0]],
             [[0, 1], [1, 0]],
             [[1, 0], [1, 0]],
             [[1, 0], [0, 0]]]

        """
        from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
        alphabet = self.domain_alphabet()
        a = next(iter(alphabet))
        table = [[a]]
        while len(table) < 2 or len(table[0]) < 2:
            table = self(table)
        shape = [(0,0), (0,1), (1,0), (1,1)]
        seeds = set_of_factors(table, shape)
        def children(factor):
            table = [[factor[0], factor[1]],
                     [factor[2], factor[3]]]
            image = self(table)
            S = set_of_factors(image, shape)
            return S
        R = RecursivelyEnumeratedSet(seeds, children)
        return [ [[a, b], [c, d]] for (a,b,c,d) in R]

    _matrix_ = incidence_matrix

def set_of_factors(table, shape, avoid_border=0):
    r"""
    Return the set of factors of given shape in the table.

    INPUT

    - ``table`` -- list of lists
    - ``shape`` -- list, list of coordinates
    - ``avoid_border`` -- integer (default: 0), the size of the border
        to avoid during the computation

    OUTPUT:

        set of tuple of integers

    EXAMPLES::

        sage: from slabbe.substitution_2d import set_of_factors
        sage: table = [[0,1,2], [3,4,5], [6,7,8]]
        sage: set_of_factors(table, shape=[(0,0), (1,0), (0,1), (1,1)])
        {(0, 3, 1, 4), (1, 4, 2, 5), (3, 6, 4, 7), (4, 7, 5, 8)}
    """
    xmin = min(x for (x,y) in shape)
    xmax = max(x for (x,y) in shape)
    ymin = min(y for (x,y) in shape)
    ymax = max(y for (x,y) in shape)
    width = len(table)
    height = len(table[0])
    S = set()
    for i in range(0-xmin+avoid_border, width-xmax-avoid_border):
        for j in range(0-ymin+avoid_border, height-ymax-avoid_border):
            pattern = tuple(table[i+x][j+y] for (x,y) in shape)
            S.add(pattern)
    return S

