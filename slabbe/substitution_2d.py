# -*- coding: utf-8 -*-
r"""
2d substitutions

EXAMPLES::

"""
#*****************************************************************************
#       Copyright (C) 2017 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import itertools

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

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from slabbe import Substitution2d
            sage: A = [[0,1],[2,3]]
            sage: B = [[4,5]]
            sage: d = {0:A, 1:B}
            sage: s = Substitution2d(d)
            sage: latex(s)
            \left\{0 : \left(\begin{array}{rr}
            1 & 3 \\
            0 & 2
            \end{array}\right), 1 : \left(\begin{array}{r}
            5 \\
            4
            \end{array}\right)\right\}
        """
        from sage.matrix.constructor import matrix
        from sage.misc.latex import latex
        d = {key:matrix.column([col[::-1] for col in table]) 
                            for key,table in self._d.items()}
        return latex(d)

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

    def __call__(self, table):
        r"""
        INPUT:

        - ``table`` -- list of list, such that table[x][y] refers to the
          tile at position (x,y) in cartesian coordinates (*not* in the
          matrix-like coordinates)

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


        TESTS::

            sage: s([])
            []
            sage: s([[0,1], [1,2]])
            Traceback (most recent call last):
            ...
            ValueError: the image of 2 in the column (=[1, 2]) has width 1
            but the image of another has width 2
        """
        columns = []
        for col in table:
            col_image = self.call_on_column(col)
            columns.extend(col_image)
        return columns

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

