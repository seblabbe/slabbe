# -*- coding: utf-8 -*-
r"""
d-dimensional Sturmian Configurations

The `d`-dimensional Sturmian configuration is a function
`\ZZ^d\to\{0,1,\dots,d\}` as follows.
Given `\boldsymbol{\alpha}=(\alpha_1,\dots,\alpha_d)\in\RR^d`, we define

.. MATH::

    \begin{array}{rccl}
    s_{\boldsymbol{\alpha}\rho}:&\ZZ^d & \to & \{0,1,\dots,d\}\\
    &\boldsymbol{n} & \mapsto &
    \sum_{i=1}^d \left(\lfloor\alpha_i+\boldsymbol{n}\cdot\boldsymbol{\alpha}+\rho\rfloor
    -\lfloor\boldsymbol{n}\cdot\boldsymbol{\alpha}+\rho\rfloor\right),\\
    \end{array}

and

.. MATH::

    \begin{array}{rccl}
    s'_{\boldsymbol{\alpha}\rho}:&\ZZ^d & \to & \{0,1,\dots,d\}\\
    &\boldsymbol{n} & \mapsto &
    \sum_{i=1}^d \left(\lceil\alpha_i+\boldsymbol{n}\cdot\boldsymbol{\alpha}+\rho\rceil
    -\lceil\boldsymbol{n}\cdot\boldsymbol{\alpha}+\rho\rceil\right),\\
    \end{array}

When `d=1`, this corresponds to Sturmian sequences, or more precisely, the
lower mechanical word and the upper mechanical word. When `d=2`, this
definition is equivalent to discrete planes as defined in
\cite{MR1782038,MR1906478}.

EXAMPLES::

    sage: z = polygen(QQ, 'z')
    sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
    sage: phi = K.gen()

    sage: from slabbe import dSturmianConfiguration
    sage: c_ceil = dSturmianConfiguration((phi^-1, sqrt(2)-1), 0, ceil_or_floor='ceil')
    sage: c_floor = dSturmianConfiguration((phi^-1, sqrt(2)-1), 0, ceil_or_floor='floor')
    sage: c_ceil.rectangular_subword_matrix(((-1,10), (-1,10)))
    [0 2 0 2 1 0 2 1 0 2 0]
    [2 0 2 1 0 2 1 0 2 0 2]
    [0 2 1 0 2 0 2 2 0 2 1]
    [2 1 0 2 0 2 1 0 2 1 0]
    [1 0 2 0 2 1 0 2 1 0 2]
    [0 2 0 2 1 0 2 0 2 2 0]
    [2 0 2 1 0 2 0 2 1 0 2]
    [0 2 1 0 2 0 2 1 0 2 1]
    [2 1 0 2 0 2 1 0 2 0 2]
    [0 2 2 0 2 1 0 2 0 2 1]
    [2 1 0 2 1 0 2 0 2 1 0]
    sage: c_floor.rectangular_subword_matrix(((-1,10), (-1,10)))
    [0 2 0 2 1 0 2 1 0 2 0]
    [2 0 2 1 0 2 1 0 2 0 2]
    [0 2 1 0 2 0 2 2 0 2 1]
    [2 1 0 2 0 2 1 0 2 1 0]
    [1 0 2 0 2 1 0 2 1 0 2]
    [0 2 0 2 1 0 2 0 2 2 0]
    [2 0 2 1 0 2 0 2 1 0 2]
    [0 2 1 0 2 0 2 1 0 2 1]
    [2 1 0 2 0 2 1 0 2 0 2]
    [1 0 2 0 2 1 0 2 0 2 1]
    [2 2 0 2 1 0 2 0 2 1 0]
    sage: window = ((0,30),(0,30))
    sage: sorted(c_floor.rectangular_subwords_matrix((2,2), window))
    [
    [0 2]  [0 2]  [1 0]  [1 0]  [2 0]  [2 1]  [2 1]  [2 2]
    [2 0], [2 1], [0 2], [2 2], [0 2], [0 2], [1 0], [1 0]
    ]
    sage: len(_)
    8

.. TODO::

    Add a method to plot the configuration as a discrete plane with
    rhombus.

"""
#*****************************************************************************
#       Copyright (C) 2020 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
from sage.modules.free_module_element import vector
from sage.functions.other import floor, ceil
from sage.matrix.constructor import matrix

class dSturmianConfiguration(object):
    r"""
    INPUT:

    - ``normal_vector`` -- tuple of coordinates
    - ``rho`` -- real number
    - ``ceil_or_floor`` -- string (default: ``'floor'``) ``'ceil'`` or
      ``'floor'``

    EXAMPLES::

        sage: from slabbe import dSturmianConfiguration
        sage: dSturmianConfiguration((.34, .72), 0)
        2-dim Sturmian Configuration: alpha=(0.340000000000000, 0.720000000000000), rho=0 using floor

    """
    def __init__(self, normal_vector, rho, ceil_or_floor='floor'):
        r"""
        Constructor.

        EXAMPLES::

            sage: from slabbe import dSturmianConfiguration
            sage: dSturmianConfiguration((.34, .72), 0)
            2-dim Sturmian Configuration: alpha=(0.340000000000000, 0.720000000000000), rho=0 using floor

        """
        self._normal_vector = vector(normal_vector)
        self._rho = rho
        self._d = len(self._normal_vector)
        self._ceil_of_floor = ceil_or_floor

    def __repr__(self):
        r"""
        Return the string representation.

        EXAMPLES::

            sage: from slabbe import dSturmianConfiguration
            sage: dSturmianConfiguration((.34, .72), 0)
            2-dim Sturmian Configuration: alpha=(0.340000000000000, 0.720000000000000), rho=0 using floor

        """
        return ("{}-dim Sturmian Configuration: "
                "alpha={}, rho={} using {}".format(self._d,
                    self._normal_vector, self._rho, self._ceil_of_floor))

    def __call__(self, n):
        r"""
        Return the letter at position ``n``.

        INPUT:

        - ``n`` -- tuple of coordinates

        EXAMPLES::

            sage: from slabbe import dSturmianConfiguration
            sage: c = dSturmianConfiguration((.34, .72), 0)
            sage: c((3,4))
            2

        """
        n = vector(n)
        scalar_product = n*self._normal_vector
        if self._ceil_of_floor == 'floor':
            return sum((floor(self._normal_vector[i] + scalar_product + self._rho) -
                        floor(scalar_product + self._rho)) for i in range(self._d))
        elif self._ceil_of_floor == 'ceil':
            return sum((ceil(self._normal_vector[i] + scalar_product + self._rho) -
                        ceil(scalar_product + self._rho)) for i in range(self._d))
        else:
            raise ValueError('ceil_or_floor (={})'.format(self._ceil_of_floor))

    def rectangular_subword(self, window):
        r"""
        Return the rectangular subword appearing in the
        given retangular window.

        INPUT:

        - ``window`` -- tuple of 2-tuples

        OUTPUT:

        list of list with euclidean coordinates

        EXAMPLES::

            sage: from slabbe import dSturmianConfiguration
            sage: c = dSturmianConfiguration((.34, .72), 0)
            sage: window = ((0,3),(0,4))
            sage: sorted(c.rectangular_subword(window))
            [[0, 2, 1, 0], [1, 0, 2, 1], [2, 1, 0, 2]]

        """
        (xmin,xmax), (ymin,ymax) = window
        return [[self((a,b)) for b in range(ymin, ymax)]
                for a in range(xmin, xmax)]

    def rectangular_subword_matrix(self, window):
        r"""
        Return the rectangular subword appearing in the
        given retangular window (as a matrix).

        INPUT:

        - ``window`` -- tuple of 2-tuples ``(start, stop)``

        OUTPUT:

        matrix

        EXAMPLES::

            sage: from slabbe import dSturmianConfiguration
            sage: c = dSturmianConfiguration((.34, .72), 0)
            sage: c.rectangular_subword_matrix(((0,3),(0,4)))
            [0 1 2]
            [1 2 0]
            [2 0 1]
            [0 1 2]

        """
        table = self.rectangular_subword(window)
        columns = [col[::-1] for col in table]
        return matrix.column(columns)

    def pattern_number_occurrences(self, shape, window, avoid_border=0):
        r"""
        Return the number of occurrences of every pattern having a given
        shape inside of a rectangular window box.

        INPUT:

        - ``shape`` -- list, list of coordinates
        - ``window`` -- tuple of 2-tuples
        - ``avoid_border`` -- integer (default: 0), the size of the border
            to avoid during the computation

        OUTPUT

        a dict where each key is a tuple giving the tiles at each
        coordinate of the shape (in the same order) and values are integers

        EXAMPLES::

            sage: z = polygen(QQ, 'z')
            sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: phi = K.gen()

        ::

            sage: from slabbe import dSturmianConfiguration
            sage: c = dSturmianConfiguration((phi^-1, phi^-2), 0)
            sage: shape = [(0,0), (1,0), (0,1), (1,1)]
            sage: window = ((0,10),(0,10))
            sage: c.pattern_number_occurrences(shape, window)
            Counter({(0, 2, 1, 0): 20, (1, 0, 2, 1): 18, (2, 1, 0, 2): 18, 
                     (2, 0, 0, 2): 13, (0, 2, 2, 0): 12})

        Totally irrational normal vector::

            sage: c = dSturmianConfiguration((phi^-1, sqrt(2)), 0)
            sage: shape = [(0,0), (1,0), (0,1), (1,1)]
            sage: window = ((0,10),(0,10))
            sage: len(c.pattern_number_occurrences(shape, window))
            8

        """
        table = self.rectangular_subword(window)
        width = len(table)
        height = len(table[0])

        shape_xmin = min(x for (x,y) in shape)
        shape_xmax = max(x for (x,y) in shape)
        shape_ymin = min(y for (x,y) in shape)
        shape_ymax = max(y for (x,y) in shape)

        from collections import Counter
        C = Counter()
        for i in range(0-shape_xmin+avoid_border, width-shape_xmax-avoid_border):
            for j in range(0-shape_ymin+avoid_border, height-shape_ymax-avoid_border):
                pattern = tuple(table[i+x][j+y] for (x,y) in shape)
                C[pattern] += 1
        return C

    def rectangular_subwords(self, sizes, window):
        r"""
        Return the list of rectangular subword appearing in the
        configuration.

        INPUT:

        - ``sizes`` -- tuple of integers
        - ``window`` -- tuple of 2-tuples

        OUTPUT:

        list of euclidean coordinate tables

        EXAMPLES::

            sage: z = polygen(QQ, 'z')
            sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: phi = K.gen()

        ::

            sage: from slabbe import dSturmianConfiguration
            sage: c = dSturmianConfiguration((phi^-1, sqrt(2)), 0)
            sage: c.rectangular_subwords((2,3), ((0,30),(0,30)))
            [[[1, 2, 3], [3, 1, 2]],
             [[2, 3, 1], [1, 2, 3]],
             [[3, 1, 3], [2, 3, 1]],
             [[1, 3, 1], [3, 1, 3]],
             [[3, 1, 2], [1, 3, 1]],
             [[1, 2, 3], [3, 1, 3]],
             [[2, 3, 2], [1, 3, 1]],
             [[3, 2, 3], [3, 1, 2]],
             [[3, 1, 3], [2, 3, 2]],
             [[1, 3, 1], [3, 2, 3]],
             [[3, 1, 3], [1, 3, 1]]]

        """
        import itertools
        shape = list(itertools.product(*[range(size) for size in sizes]))
        C = self.pattern_number_occurrences(shape, window)
        patterns = C.keys()
        L = []
        for pattern in patterns:
            table = [[None for _ in range(sizes[1])] for _ in range(sizes[0])]
            for i,(x,y) in enumerate(shape):
                table[x][y] = pattern[i] 
            L.append(table)
        return L

    def rectangular_subwords_matrix(self, sizes, window):
        r"""
        Return the list of rectangular subword appearing under the form of
        matrices.

        INPUT:

        - ``sizes`` -- tuple of integers
        - ``window`` -- tuple of 2-tuples

        OUTPUT:

        list of matrices

        EXAMPLES::

            sage: z = polygen(QQ, 'z')
            sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: phi = K.gen()

            sage: from slabbe import dSturmianConfiguration
            sage: c = dSturmianConfiguration((phi^-1, sqrt(2)), 0)
            sage: c.rectangular_subwords_matrix((2,3), ((0,30),(0,30)))
            [
            [3 2]  [1 3]  [3 1]  [1 3]  [2 1]  [3 3]  [2 1]  [3 2]  [3 2]  [1 3]
            [2 1]  [3 2]  [1 3]  [3 1]  [1 3]  [2 1]  [3 3]  [2 1]  [1 3]  [3 2]
            [1 3], [2 1], [3 2], [1 3], [3 1], [1 3], [2 1], [3 3], [3 2], [1 3],
            <BLANKLINE>
            [3 1]
            [1 3]
            [3 1]
            ]
            sage: len(_)
            11

        """
        L = self.rectangular_subwords(sizes, window)
        return [matrix.column([col[::-1] for col in table]) for table in L]

