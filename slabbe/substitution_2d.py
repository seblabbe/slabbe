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
from sage.misc.superseded import deprecated_function_alias

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

    def __eq__(self, other):
        r"""
        INPUT:

        - ``other`` -- substitution 2d

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[2,3]]
            sage: B = [[4,5]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: s == Substitution2d(d)
            True
        """
        return isinstance(other, Substitution2d) and self._d == other._d

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
        from sage.modules.free_module_element import vector
        from sage.misc.latex import latex
        from sage.misc.latex import LatexExpr
        from sage.calculus.var import var
        from slabbe.matrices import map_coefficients_to_variable_index
        lines = []
        lines.append(r'\begin{{array}}{{{}}}'.format(align*ncolumns))
        for i,(key,table) in enumerate(self._d.items()):
            M = matrix.column([col[::-1] for col in table])
            if M.nrows() == 1:
                # do not latex array for horizontal vectors
                M = vector(M)
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
        d = {key:[[a for a in col]] for key,col in s.items()}
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
            Substitution 2d: {5: [[0]], 6: [[1]], 7: [[2]], 8: [[3]], 9: [[4]]}
            sage: u * t
            Substitution 2d: {0: [[0, 1], [2, 3]], 1: [[1, 0], [4, 3]]}
        """
        d = {key:[[d[key]]] for key in sorted(d)}
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
            ValueError: the image of the row contains columns of different height (={2, 3})
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

    def call_on_column(self, column, heights=None):
        r"""
        INPUT:

        - ``column`` -- list
        - ``heights`` -- None or list (default: ``None``)

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

        It can compute the image of columns with ``None`` as entries::

            sage: s.call_on_column([0,None], heights=[2,3])
            [[0, 1, None, None, None], [2, 3, None, None, None]]
            sage: s.call_on_column([0,None], heights=[2,2])
            [[0, 1, None, None], [2, 3, None, None]]
            sage: s.call_on_column([None], heights=[3])
            [[None, None, None]]

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

        # compute the width of the image
        for a in column:
            if not a is None:
                width = len(self._d[a])
                break
        else:
            width = 1

        # compute the image
        rep = [[] for _ in range(width)]
        for i,a in enumerate(column):
            if a is None:
                if heights is None:
                    raise ValueError("the {}-th element of given column is None,"
                            " so you must provide the list ``heights`` as"
                            " input".format(i))
                height = heights[i]
                for j in range(width):
                    rep[j].extend([None]*height)
            else:
                image_a = self._d[a]
                if not len(image_a) == width:
                    raise ValueError("the image of {} in the column (={}) has"
                            " width {} but the image of another"
                            " has width {}".format(a, column, len(image_a), width))
                if heights and not len(image_a[0]) == heights[i]:
                    raise ValueError("the image of {} in the column (={}) has"
                            " height {} but inputs ``heights`` says it should be"
                            " {}".format(a, column, len(image_a[0]), heights[i]))
                for j,col in enumerate(image_a):
                    rep[j].extend(col)
        return rep

    def __call__(self, table, order=1):
        r"""
        INPUT:

        - ``table`` -- list of list, such that table[x][y] refers to the
          tile at position (x,y) in cartesian coordinates (*not* in the
          matrix-like coordinates)
        -  ``order`` - integer or plus ``Infinity`` (default: 1)

        TODO: implement another call on table with matrix like coordinates

        OUTPUT:

            list of columns

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

        Works if some ``None`` entries are involved::

            sage: s([[None]])
            [[None]]
            sage: s([[None],[0]])
            [[None, None], [0, 1], [2, 0]]
            sage: s([[None,0],[0,2]])
            [[None, None, 0, 1], [None, None, 2, 0], [0, 1, 1, 2], [2, 0, 1, 1]]

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
        if not table:
            return []

        if order == 1:
            # compute the heights of the image of each row
            # in case their are some None entries involved
            heights = []
            for i in range(len(table[0])):
                for col in table:
                    if col[i] is not None:
                        height = len(self._d[col[i]][0])
                        break
                else:
                    height = 1
                heights.append(height)
            # compute the image
            columns = []
            for col in table:
                col_image = self.call_on_column(col, heights)
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
            ValueError: codomain alphabet of other (={4, 5}) must be
            included in domain alphabet of self (={0, 1})
        """
        if not self.domain_alphabet() >= other.codomain_alphabet():
            raise ValueError("codomain alphabet of other (={}) must be included"
                " in domain alphabet of self (={})".format(other.codomain_alphabet(),
                                                            self.domain_alphabet()))
        d = {}
        for a,image_a in other._d.items():
            d[a] = self(image_a)
        return Substitution2d(d)

    def reversal(self):
        r"""
        Return the reversal of self.

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[1,2],[3,4]]
            sage: B = [[5,6],[7,8]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: s.reversal()
            Substitution 2d: {0: [[4, 3], [2, 1]], 1: [[8, 7], [6, 5]]}
        """
        d = {}
        for a,image_a in self._d.items():
            d[a] = [col[::-1] for col in reversed(image_a)]
        return Substitution2d(d)

    def inverse(self):
        r"""
        Return the inverse of self (when self is a permutation).

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: d = {0:7, 1:8}
            sage: s = Substitution2d.from_permutation(d)
            sage: s
            Substitution 2d: {0: [[7]], 1: [[8]]}
            sage: s.inverse()
            Substitution 2d: {7: [[0]], 8: [[1]]}

        TESTS::

            sage: s = Substitution2d({8: [[1]], 7: [[0,1]]})
            sage: s.inverse()
            Traceback (most recent call last):
            ...
            ValueError: self must be a permutation but image of 7 is [[0, 1]]

        """
        d = {}
        for a,image_a in self._d.items():
            if len(image_a) > 1 or len(image_a[0]) > 1:
                raise ValueError('self must be a permutation but image '
                                 'of {} is {}'.format(a, image_a))
            d[image_a[0][0]] = a 
        return Substitution2d.from_permutation(d)

    def letter_to_letter_dict(self, pos=(0,0)):
        r"""
        Return the inverse of self (when self is a permutation).

        INPUT:

        - ``pos`` -- tuple (default:``(0,0)``), tuple of two integers

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[2,3]]
            sage: B = [[4,5]]
            sage: s = Substitution2d({0:A, 1:B})
            sage: s
            Substitution 2d: {0: [[0, 1], [2, 3]], 1: [[4, 5]]}
            sage: s.letter_to_letter_dict(pos=(0,0))
            {0: 0, 1: 4}

        """
        x,y = pos
        return {a:image_a[x][y] for a,image_a in self._d.items()}

    def apply_matrix_transformation(self, M):
        r"""
        INPUT:

        - ``M`` -- matrix in SL(2,Z)

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[0,1]]
            sage: B = [[1,0],[1,1]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: M = matrix(2, (1,1,0,1))
            sage: s
            Substitution 2d: {0: [[0, 1], [0, 1]], 1: [[1, 0], [1, 1]]}
            sage: s.apply_matrix_transformation(M)
            Substitution 2d: {0: [[0, None], [0, 1], [None, 1]], 1: [[1, None], [1, 0], [None, 1]]}

        """
        from slabbe.wang_tiles import WangTiling
        d = {}
        for a,image_a in self._d.items():
            d[a] = WangTiling(image_a,tiles=[]).apply_matrix_transformation(M).table()
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

        Blank ``None`` are ignored::

            sage: A = [[5,6],[7,8]]
            sage: B = [[6,5],[9,None]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: s.codomain_alphabet()
            {5, 6, 7, 8, 9}
        """
        s = set()
        for a,image_a in self._d.items():
            for column in image_a:
                s.update(column)
        if None in s:
            s.remove(None)
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
    def wang_tikz(self, domain_tiles, codomain_tiles, domain_color=None,
            codomain_color=None, size=1, scale=1, font=r'\normalsize',
            rotate=None, label_shift=.2, id=True, edges=True, 
            ncolumns=4, direction='right', extra_space=1):
        r"""
        Return the tikz code showing what the substitution A->B* does on
        Wang tiles.

        INPUT:

        - ``domain_tiles`` -- tiles of the domain
        - ``codomain_tiles`` -- tiles of the codomain
        - ``domain_color`` -- dict (default: ``None``) from tile values -> tikz colors
        - ``codomain_color`` -- dict (default: ``None``) from tile values -> tikz colors
        - ``size`` -- number (default: ``1``), size of the tile
        - ``scale`` -- number (default: ``1``), scale of tikzpicture
        - ``font`` -- string (default: ``r'\normalsize'``
        - ``rotate`` -- list or ``None`` (default:``None``) list of four angles
          in degrees like ``(0,0,0,0)``, the rotation angle to apply to each
          label of Wang tiles. If ``None``, it performs a 90 degres rotation
          for left and right labels taking more than one character.
        - ``label_shift`` -- number (default: ``.2``) translation distance
          of the label from the edge
        - ``id`` -- boolean (default: ``True``), presence of the tile id
        - ``ncolumns`` -- integer (default: ``4``)
        - ``edges`` -- bool (default: ``True``) 
        - ``direction`` -- string (default: ``'right'``) or ``'down'``
        - ``extra_space`` -- number (default: ``1``), space between the
          tile and its image

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
            sage: tikz = s.wang_tikz(domain_tiles, codomain_tiles, rotate=(90,0,90,0))
            sage: _ = tikz.pdf(view=False)      # long time

        Applying a transformation matrix::

            sage: M = matrix(2, [1,1,0,1])
            sage: sM = s.apply_matrix_transformation(M)
            sage: tikz = sM.wang_tikz(domain_tiles, codomain_tiles)
            sage: _ = tikz.pdf(view=False)      # long time

        Down direction::

            sage: tikz = s.wang_tikz(domain_tiles, codomain_tiles,
            ....:                      direction='down')
            sage: _ = tikz.pdf(view=False)      # long time
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
        lines.append(r'\begin{tikzpicture}')
        lines.append(r'[scale={}]'.format(scale))
        lines.append(r'\tikzstyle{{every node}}=[font={}]'.format(font))
        lines.append(r'\matrix{')

        for i,a in enumerate(self._d):

            lines.append(r'\node{')
            lines.append(r'\begin{tikzpicture}')
            lines.append(r'[scale={}]'.format(scale))
            lines.append(r'\tikzstyle{{every node}}=[font={}]'.format(font))

            desubstituted_tile = domain_tiles[a] 
            this_id = a if id else None
            new_lines = tile_to_tikz(desubstituted_tile, (0,0),
                    color=domain_color, id=this_id, sizex=size, sizey=size,
                    rotate=rotate, label_shift=label_shift,
                    right_edges=edges, top_edges=edges, left_edges=edges,
                    bottom_edges=edges)
            new_lines.insert(0, r'\begin{tikzpicture}')
            new_lines.append(r'\end{tikzpicture}')
            new_lines = '\n'.join(new_lines)
            lines.append(r'\node (A) at (0,0) {{{}}};'.format(new_lines))

            image_a = self._d[a]
            tiling = WangTiling(image_a, codomain_tiles, codomain_color)
            tiling_tikz = tiling.tikz(color=codomain_color, font=font,
                    rotate=rotate, label_shift=label_shift, scale=scale,
                    edges=edges, id=id, size=size)

            size_image_x = len(image_a)
            size_image_y = len(image_a[0])

            if direction == 'right':
                xshift = extra_space + .5 * size_image_x
                yshift = 0
            elif direction == 'down':
                xshift = 0
                yshift = -extra_space - .5 * size_image_y
            lines.append(r'\node (B) at ({},{}) {{{}}};'.format(xshift, yshift,
                                             tiling_tikz.tikz_picture_code()))

            lines.append(r'\draw[-to,very thick] (A) edge (B);')

            lines.append(r'\end{tikzpicture}')
            lines.append(r'};')

            if (i+1) == len(self._d):
                lines.append(r'\\')
            elif (i+1) % ncolumns == 0:
                lines.append(r'\\')
            else:
                lines.append(r'&')

        lines.append(r'};')  # end of \matrix{
        lines.append(r'\end{tikzpicture}')

        from slabbe import TikzPicture
        return TikzPicture('\n'.join(lines), usetikzlibrary=['positioning'])

    def wang_tiles_codomain_tikz(self, codomain_tiles, color=None,
            size=1, scale=1, font=r'\normalsize', rotate=None,
            id=True, label=True, label_shift=.2, edges=True,
            ncolumns=4, direction='right'):
        r"""
        Return the tikz code of the image of the letters as a table of
        tikz tilings.

        INPUT:

        - ``domain_tiles`` -- tiles of the domain
        - ``codomain_tiles`` -- tiles of the codomain
        - ``domain_color`` -- dict (default: ``None``) from tile values -> tikz colors
        - ``codomain_color`` -- dict (default: ``None``) from tile values -> tikz colors
        - ``size`` -- number (default: ``1``), size of the tile
        - ``scale`` -- number (default: ``1``), scale of tikzpicture
        - ``font`` -- string (default: ``r'\normalsize'``
        - ``rotate`` -- list or ``None`` (default:``None``) list of four angles
          in degrees like ``(0,0,0,0)``, the rotation angle to apply to each
          label of Wang tiles. If ``None``, it performs a 90 degres rotation
          for left and right labels taking more than one character.
        - ``id`` -- boolean (default: ``True``), presence of the tile id
        - ``label`` -- boolean (default: ``True``) 
        - ``label_shift`` -- number (default: ``.2``) translation distance
          of the label from the edge
        - ``edges`` -- bool (default: ``True``) 
        - ``ncolumns`` -- integer (default: ``4``)

        OUTPUT:

            tikzpicture

        EXAMPLES::

            sage: from slabbe import WangTileSet, Substitution2d
            sage: A = [[0,1,2],[1,0,0]]
            sage: B = [[0,1,2]]
            sage: d = {4:A, 5:B}
            sage: s = Substitution2d(d)
            sage: codomain_tiles = [(0,3,1,4), (1,4,0,3), (5,6,7,8)]
            sage: W = WangTileSet(codomain_tiles)
            sage: t = s.wang_tiles_codomain_tikz(W)
            sage: _ = t.pdf(view=False)

        """
        from slabbe.wang_tiles import tile_to_tikz, WangTileSet, WangTiling
        for a in self.codomain_alphabet():
            try:
                codomain_tiles[a]
            except IndexError:
                raise ValueError("codomain_alphabet={}, but tiles are"
                        " {}".format(self.codomain_alphabet(),
                        codomain_tiles))

        lines = []
        lines.append(r'\begin{tikzpicture}')
        lines.append(r'[scale={}]'.format(scale))
        lines.append(r'\tikzstyle{{every node}}=[font={}]'.format(font))
        lines.append(r'\matrix{')
        for i,a in enumerate(self._d):

            image_a = self._d[a]
            tiling = WangTiling(image_a, codomain_tiles, color)
            tikz = tiling.tikz(color=color, font=font, id=id, label=label,
                    rotate=rotate, label_shift=label_shift, scale=scale,
                    edges=edges, size=size)

            lines.append(r'\node{{{}}};'.format(tikz.tikz_picture_code()))

            if (i+1) == len(self._d):
                lines.append(r'\\')
            elif (i+1) % ncolumns == 0:
                lines.append(r'\\')
            else:
                lines.append(r'&')

        lines.append(r'};')  # end of \matrix{
        lines.append(r'\end{tikzpicture}')

        from slabbe import TikzPicture
        return TikzPicture('\n'.join(lines))

    def list_2x2_factors(self, F=None):
        r"""
        Return the list of 2x2 factors in the associated substitutive
        shift. If a list of factors ``F`` is given, it restrict to the
        factors inside the image of ``F``.

        INPUT:

        - ``self`` -- expansive and primitive 2d substitution
        - ``F`` -- list of factors in the domain or ``None``, if given the
          output is restricted to the factors in ``F``

        OUTPUT:

            list of tables

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[0,1]]
            sage: B = [[1,0],[1,1]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: sorted(s.list_2x2_factors())
            [[[0, 0], [1, 0]],
             [[0, 1], [0, 1]],
             [[0, 1], [1, 0]],
             [[0, 1], [1, 1]],
             [[1, 0], [0, 0]],
             [[1, 0], [0, 1]],
             [[1, 0], [1, 0]],
             [[1, 0], [1, 1]],
             [[1, 1], [0, 0]],
             [[1, 1], [0, 1]],
             [[1, 1], [1, 0]],
             [[1, 1], [1, 1]]]

        Restricting to the images of some factors::

            sage: sorted(s.list_2x2_factors([A]))
            [[[0, 1], [0, 1]], [[1, 0], [1, 1]], [[1, 1], [1, 0]], [[1, 1], [1, 1]]]
            sage: sorted(s.list_2x2_factors([B]))
            [[[0, 0], [1, 0]],
             [[0, 1], [0, 1]],
             [[0, 1], [1, 0]],
             [[0, 1], [1, 1]],
             [[1, 0], [0, 1]],
             [[1, 0], [1, 1]],
             [[1, 1], [1, 0]]]
            sage: sorted(s.list_2x2_factors([A,B]))
            [[[0, 0], [1, 0]],
             [[0, 1], [0, 1]],
             [[0, 1], [1, 0]],
             [[0, 1], [1, 1]],
             [[1, 0], [0, 1]],
             [[1, 0], [1, 1]],
             [[1, 1], [1, 0]],
             [[1, 1], [1, 1]]]
            sage: s.list_2x2_factors([])
            []

        """
        if not self.codomain_alphabet() <= self.domain_alphabet():
            raise ValueError("codomain alphabet (={}) is not a subset of the"
                    " domain alphabet (={})".format(self.codomain_alphabet(),
                                              self.domain_alphabet()))
        if F is None:
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
        else:
            shape = [(0,0), (0,1), (1,0), (1,1)]
            S = set()
            for table in F:
                S.update(set_of_factors(self(table), shape))
            return [ [[a, b], [c, d]] for (a,b,c,d) in S]

    def lines_alphabet(self, direction='horizontal'):
        r"""
        Return the possible alphabets on lines, i.e., the possible alphabet
        of letters that we see on a given line.

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[0,1]]
            sage: B = [[1,0],[1,1]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: sorted(s.lines_alphabet())
            [(0,), (0, 1), (1,)]
            sage: sorted(s.lines_alphabet(direction='vertical'))
            [(0, 1), (1,)]

        """
        if direction == 'horizontal':
            def children(node):
                T = self.call_on_row(node)
                return [tuple(sorted(set(row))) for row in zip(*T)]
        elif direction == 'vertical':
            def children(node):
                T = self.call_on_column(node)
                return [tuple(sorted(set(column))) for column in T]
        else:
            raise ValueError("direction (={}) must be 'vertical' or"
                    " 'horizontal'".format(direction))

        from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
        seeds = [(a,) for a in self.domain_alphabet()]
        R = RecursivelyEnumeratedSet(seeds, children)
        G = R.to_digraph(multiedges=False)
        result = []
        for s in G.strongly_connected_components_subgraphs():
            if s.num_edges() > 0:
                for v in s.vertices():
                    result.append(v)
        return result

    def prolongable_seeds_graph(self):
        r"""
        Return the directed graph of 2x2 factors where (u,v) is an edge if
        v is the seed at the origin of the image of u under self.

        OUTPUT:

            list of tuple (2x2 matrix, integer)

        EXAMPLES::

            sage: d = {0: [[17]],
            ....:  1: [[16]],
            ....:  2: [[15], [11]],
            ....:  3: [[13], [9]],
            ....:  4: [[17], [8]],
            ....:  5: [[16], [8]],
            ....:  6: [[15], [8]],
            ....:  7: [[14], [8]],
            ....:  8: [[14, 6]],
            ....:  9: [[17, 3]],
            ....:  10: [[16, 3]],
            ....:  11: [[14, 2]],
            ....:  12: [[15, 7], [11, 1]],
            ....:  13: [[14, 6], [11, 1]],
            ....:  14: [[13, 7], [9, 1]],
            ....:  15: [[12, 6], [9, 1]],
            ....:  16: [[18, 5], [10, 1]],
            ....:  17: [[13, 4], [9, 1]],
            ....:  18: [[14, 2], [8, 0]]}
            sage: from slabbe import Substitution2d
            sage: omega = Substitution2d(d)
            sage: G = omega.prolongable_seeds_graph()
            sage: G
            Looped digraph on 50 vertices

        """
        from sage.matrix.constructor import matrix
        from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
        seeds = self.list_2x2_factors(F=None)
        seeds = [matrix.column(c[::-1] for c in columns) for columns in seeds]
        for m in seeds: m.set_immutable()
        def children(m):
            b,d,a,c = m.list()
            A = self([[a]])[-1][-1] # bottom left
            B = self([[b]])[-1][0] # top left
            C = self([[c]])[0][-1] # bottom right
            D = self([[d]])[0][0] # top right
            M = matrix.column(((B,A),(D,C)))
            M.set_immutable()
            return [M]
        R = RecursivelyEnumeratedSet(seeds, children)
        return R.to_digraph(multiedges=False)

    prolongable_origins = deprecated_function_alias(123456, prolongable_seeds_graph)
    def prolongable_seeds_list(self):
        r"""
        Return the list of seed which are prolongable for some power of
        self.

        OUTPUT:

            list of cycles

        EXAMPLES::

            sage: d = {0: [[17]],
            ....:  1: [[16]],
            ....:  2: [[15], [11]],
            ....:  3: [[13], [9]],
            ....:  4: [[17], [8]],
            ....:  5: [[16], [8]],
            ....:  6: [[15], [8]],
            ....:  7: [[14], [8]],
            ....:  8: [[14, 6]],
            ....:  9: [[17, 3]],
            ....:  10: [[16, 3]],
            ....:  11: [[14, 2]],
            ....:  12: [[15, 7], [11, 1]],
            ....:  13: [[14, 6], [11, 1]],
            ....:  14: [[13, 7], [9, 1]],
            ....:  15: [[12, 6], [9, 1]],
            ....:  16: [[18, 5], [10, 1]],
            ....:  17: [[13, 4], [9, 1]],
            ....:  18: [[14, 2], [8, 0]]}
            sage: from slabbe import Substitution2d
            sage: omega = Substitution2d(d)
            sage: omega.prolongable_seeds_list()
            [[
             [ 9 14]  [17 13]
             [ 1  6], [16 15]
             ],
             [
             [ 9 14]  [17 13]
             [ 8 16], [ 6  5]
             ],
             [
             [10 12]  [16 15]
             [ 9 14], [ 3  7]
             ],
             [
             [10 14]  [16 13]
             [11 17], [ 2  4]
             ]]

        """
        G = self.prolongable_seeds_graph()
        return [cycle[:-1] for cycle in G.all_simple_cycles()]

    _matrix_ = incidence_matrix
    def relabel_domain(self, other):
        r"""
        Return a permutation p such that self*p == other, if it exists.

        INPUT:

        - ``other`` -- substitution 2d

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[0,1]]
            sage: B = [[1,0],[1,1]]
            sage: s = Substitution2d({0:A, 1:B})
            sage: t = Substitution2d({7:A, 8:B})
            sage: s.relabel_domain(t)
            Substitution 2d: {7: [[0]], 8: [[1]]}

        TESTS::

            sage: s = Substitution2d({0:A, 1:B})
            sage: s.relabel_domain(s)
            Substitution 2d: {0: [[0]], 1: [[1]]}

        ::

            sage: s = Substitution2d({0:A, 1:B})
            sage: t = Substitution2d({7:A, 8:B, 9:[[4]]})
            sage: t.relabel_domain(s)
            Traceback (most recent call last):
            ...
            ValueError: image of letter 9 is [[4]] and is not in other
            sage: s.relabel_domain(t)
            Traceback (most recent call last):
            ...
            AssertionError: problem: self * p == other not satisfied
        """
        table_to_tuple = lambda table:tuple(tuple(col) for col in table)
        other_inv = {table_to_tuple(val):key for key,val in other._d.items()}
        p = {}
        for key,val in self._d.items():
            val_ = table_to_tuple(val)
            if val_ in other_inv:
                p[other_inv[val_]] = key 
            else:
                raise ValueError('image of letter {} is {} and is not in other'.format(key, val))
        p = Substitution2d.from_permutation(p)
        assert self * p == other, "problem: self * p == other not satisfied"
        return p

    def fixed_point_tikz(self, seed, niterations=3):
        r"""
        Return a tikz representation of a fixed point defined by the give seed.
        In the image, rectangular boxes indicate the i-th image of each seed.

        INPUT:

        - ``seed`` -- 2x2 matrix
        - ``niterations`` -- (default:``3``), number of iterations

        OUTPUT

            tikz picture

        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[1,1],[1,1]]
            sage: B = [[0,0]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: tikz = s.fixed_point_tikz([[0,0],[0,0]])
            sage: _ = tikz.pdf()                              # not tested

        The substitution ``s`` is not prolongable, so the boxes in the
        image obtained from the square of ``s`` might better::

            sage: s2 = s*s
            sage: tikz = s2.fixed_point_tikz([[0,0],[0,0]])
            sage: _ = tikz.pdf()                              # not tested

        """
        if isinstance(seed, list):
            [[C,B],[D,A]] = seed
        else:
            # we assume it is a 2x2 matrix
            B,A,C,D = seed.list()
        quadrants = {}
        quadrants[1] = [[[A]]]
        quadrants[2] = [[[B]]]
        quadrants[3] = [[[C]]]
        quadrants[4] = [[[D]]]
        for i in [1,2,3,4]:
            for order in range(niterations):
                previous = quadrants[i][-1]
                quadrants[i].append(self(previous))
        width_heigth = lambda quad : (len(quad),len(quad[0]))
        dimq1x,dimq1y = width_heigth(quadrants[1][-1])
        dimq2x,dimq2y = width_heigth(quadrants[2][-1])
        dimq3x,dimq3y = width_heigth(quadrants[3][-1])
        dimq4x,dimq4y = width_heigth(quadrants[4][-1])
        lines = []
        lines.append(r"\begin{tikzpicture}")
        for i,col in enumerate(quadrants[1][-1]):
            for j,a in enumerate(col):
                coords = (i+.5,j+.5)
                lines.append(r"\node at {} {{{}}};".format(coords,a))
        for i,col in enumerate(quadrants[2][-1]):
            for j,a in enumerate(col):
                coords = (i-dimq2x+.5,j+.5)
                lines.append(r"\node at {} {{{}}};".format(coords,a))
        for i,col in enumerate(quadrants[3][-1]):
            for j,a in enumerate(col):
                coords = (i-dimq3x+.5,j-dimq3y+.5)
                lines.append(r"\node at {} {{{}}};".format(coords,a))
        for i,col in enumerate(quadrants[4][-1]):
            for j,a in enumerate(col):
                coords = (i+.5,j-dimq4y+.5)
                lines.append(r"\node at {} {{{}}};".format(coords,a))
        # axes
        lines.append(r"\draw {} -- {};".format((0,-dimq4y), (0,dimq1y)))
        lines.append(r"\draw {} -- {};".format((-dimq2x,0), (dimq1x,0)))
        # boxes
        for order in range(niterations):
            w,h = width_heigth(quadrants[3][order])
            lower_left = (-w, -h)
            upper_right = width_heigth(quadrants[1][order])
            lines.append(r"\draw {} rectangle {};".format(lower_left, upper_right))
        lines.append(r"\end{tikzpicture}")

        from slabbe import TikzPicture
        return TikzPicture('\n'.join(lines))

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

