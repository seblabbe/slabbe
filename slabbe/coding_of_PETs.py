# -*- coding: utf-8 -*-
r"""
Coding of Polyhedron exchange transformations (PETs)

Coding of Z^2-actions given by a tuple of Polyhedron exchange
transformations (PETs) and one polyhedron partition

EXAMPLES:

A polyhedron partition::

    sage: from slabbe import PolyhedronPartition
    sage: from slabbe import PolyhedronExchangeTransformation as PET

AUTHORS:

- Sébastien Labbé, January 2020, initial version
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
import itertools
from copy import copy
from sage.misc.cachefunc import cached_method
from sage.geometry.polyhedron.constructor import Polyhedron
from slabbe import PolyhedronPartition
from slabbe import PolyhedronExchangeTransformation as PET

class PETsCoding(object):
    r"""
    Coding of a tuple of commuting PETs by a partition

    INPUT:

    - ``PETs`` -- tuple of PolyhedronExchangeTransformation
    - ``partition`` -- polyhedron partition

    EXAMPLES::

        sage: from slabbe import PolyhedronPartition
        sage: h = 1/3
        sage: p = Polyhedron([(0,h),(0,1),(h,1)])
        sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
        sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
        sage: s = Polyhedron([(h,0), (1,0), (1,h)])
        sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
        sage: from slabbe import PolyhedronExchangeTransformation as PET
        sage: base = identity_matrix(2)
        sage: Re1 = PET.toral_translation(base, vector((2/3, 0)))
        sage: Re2 = PET.toral_translation(base, vector((0, 1/4)))
        sage: from slabbe import PETsCoding
        sage: PETsCoding((Re1,Re2), P)
        Coding of PETs (Polyhedron Exchange Transformation of
        Polyhedron partition of 2 atoms with 2 letters
        with translations {0: (2/3, 0), 1: (-1/3, 0)}, Polyhedron Exchange
        Transformation of
        Polyhedron partition of 2 atoms with 2 letters
        with translations {0: (0, 1/4), 1: (0, -3/4)}) by partition
        Polyhedron partition of 4 atoms with 4 letters

    """
    def __init__(self, PETs, partition):
        self._PETs = PETs
        self._partition = partition

    def __repr__(self):
        return "Coding of PETs {} by partition {}".format(self._PETs,
                self._partition)

    def ambient_space(self):
        r"""
        TODO: Maybe we want to make the union with the ambient space of the PETs?
        """
        return self._partition.ambient_space()

    def configuration(self, x0):
        raise NotImplementedError

    def pattern(self, x0, sizes):
        r"""
        Return the pattern obtained as the coding of the orbit of some
        starting point by the application of the PETs a certain number of
        times given by the tuple of sizes.

        TODO: add a input direction when the point lies in more than one
        atoms

        INPUT:

        - ``x0`` -- point in the domain of the partition
        - ``sizes`` -- tuple of integers

        OUTPUT:

            list of lists (using cartesian coordinates)

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: from slabbe import PolyhedronExchangeTransformation as PET
            sage: base = identity_matrix(2)
            sage: Re1 = PET.toral_translation(base, vector((2/3, 0)))
            sage: Re2 = PET.toral_translation(base, vector((0, 1/4)))
            sage: from slabbe import PETsCoding
            sage: X_P_R = PETsCoding((Re1,Re2), P)
            sage: X_P_R.pattern((1/7,1/7), (3,5))
            [[1, 1, 0, 0, 1], [3, 2, 2, 2, 3], [2, 2, 2, 2, 2]]

        When the point lies on the boundary, it currently raises an error::

            sage: X_P_R.pattern((0,0), (3,5))
            Traceback (most recent call last):
            ...
            ValueError: polyhedron p whose vertices are (A vertex at (0, 3/4),)
            lies in more than one atoms (=[0, 1])

        """
        if len(sizes) != 2:
            raise NotImplementedError("we assume len(sizes) is 2 for now")

        Re1 = self._PETs[0]
        Re2 = self._PETs[1]

        x0 = self.ambient_space()(x0)

        horizontal_orbit_x0 = [x0]
        for i in range(sizes[0]-1):
            p = Re1(horizontal_orbit_x0[-1])
            horizontal_orbit_x0.append(p)

        table = []
        for p in horizontal_orbit_x0:
            a = self._partition.code(Polyhedron([p]))
            column = [a]
            for j in range(sizes[1]-1):
                p = Re2(p)
                a = self._partition.code(Polyhedron([p]))
                column.append(a)
            table.append(column)
        return table

    def cylinder(self, pattern):
        r"""
        Return the coding region of the pattern.

        INPUT:

        - ``pattern`` -- list of lists or dict of positions to code

        OUTPUT:

            polyhedron partition (containing probably only one atom, or
            more to handle the case of union of polyhedrons)

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: from slabbe import PolyhedronExchangeTransformation as PET
            sage: base = identity_matrix(2)
            sage: Re1 = PET.toral_translation(base, vector((2/3, 0)))
            sage: Re2 = PET.toral_translation(base, vector((0, 1/4)))
            sage: from slabbe import PETsCoding
            sage: X_P_R = PETsCoding((Re1,Re2), P)
            sage: pattern = [[1, 1, 0, 0, 1], [3, 2, 2, 2, 3], [2, 2, 2, 2, 2]]
            sage: C = X_P_R.cylinder(pattern)
            sage: C
            Polyhedron partition of 1 atoms with 1 letters
            sage: atom = C.atoms()[0]
            sage: atom
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 6 vertices
            sage: atom.vertices()
            (A vertex at (1/9, 1/18),
             A vertex at (5/24, 1/4),
             A vertex at (0, 0),
             A vertex at (1/6, 1/4),
             A vertex at (1/18, 7/36),
             A vertex at (0, 1/12))
            sage: v = vector((1/7, 1/7))
            sage: v.set_immutable()
            sage: v in atom
            True

        """
        if isinstance(pattern, dict):
            raise NotImplementedError
        elif not isinstance(pattern, (list, tuple)):
            raise TypeError("pattern(={}) must be a list "
                    "of lists or a dict".format(pattern))

        Re1 = self._PETs[0]
        Re2 = self._PETs[1]

        Re1_inv = Re1.inverse()
        Re2_inv = Re2.inverse()

        region = self._partition
        for i,column in enumerate(pattern):
            for j,a in enumerate(column):
                translated_back = Re1_inv(self._partition[a], niterations=i)
                translated_back = Re2_inv(translated_back, niterations=j)
                region = region.refinement(translated_back)

        return region


    def partition_for_patterns(self, sizes):
        r"""
        Return the coding region of the pattern.

        INPUT:

        - ``pattern`` -- list of lists or dict of the form
          ``{positions:code}``

        OUTPUT:

            polyhedron partition (containing probably only one atom, or
            more to handle the case of union of polyhedrons)

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: from slabbe import PolyhedronExchangeTransformation as PET
            sage: base = identity_matrix(2)
            sage: Re1 = PET.toral_translation(base, vector((2/3, 0)))
            sage: Re2 = PET.toral_translation(base, vector((0, 1/4)))
            sage: from slabbe import PETsCoding
            sage: X_P_R = PETsCoding((Re1,Re2), P)
            sage: X_P_R.partition_for_patterns((2,2))
            (Polyhedron partition of 24 atoms with 24 letters,
             {0: [[0, 0], [2, 2]],
              1: [[0, 1], [2, 2]],
              2: [[0, 1], [2, 3]],
              3: [[1, 0], [2, 2]],
              4: [[1, 0], [3, 2]],
              5: [[1, 1], [2, 2]],
              6: [[1, 1], [3, 2]],
              7: [[1, 1], [3, 3]],
              8: [[1, 1], [2, 3]],
              9: [[2, 2], [0, 0]],
              10: [[2, 2], [1, 0]],
              11: [[2, 2], [1, 1]],
              12: [[2, 2], [2, 2]],
              13: [[2, 2], [0, 1]],
              14: [[2, 2], [1, 1]],
              15: [[2, 2], [2, 2]],
              16: [[2, 3], [0, 1]],
              17: [[2, 3], [1, 1]],
              18: [[2, 3], [2, 2]],
              19: [[2, 3], [2, 3]],
              20: [[3, 2], [1, 1]],
              21: [[3, 2], [2, 2]],
              22: [[3, 2], [3, 2]],
              23: [[3, 3], [3, 2]]})
            sage: X_P_R.partition_for_patterns((1,3))
            (Polyhedron partition of 18 atoms with 18 letters,
             {0: [[0, 0, 0]],
              1: [[0, 0, 1]],
              2: [[0, 1, 0]],
              3: [[0, 1, 1]],
              4: [[1, 0, 0]],
              5: [[1, 0, 1]],
              6: [[1, 1, 0]],
              7: [[1, 1, 1]],
              8: [[1, 1, 1]],
              9: [[1, 1, 1]],
              10: [[2, 2, 2]],
              11: [[2, 2, 2]],
              12: [[2, 2, 2]],
              13: [[2, 2, 3]],
              14: [[2, 3, 2]],
              15: [[2, 3, 3]],
              16: [[3, 2, 2]],
              17: [[3, 3, 2]]})

        """
        if len(sizes) != 2:
            raise NotImplementedError("we assume len(sizes) is 2 for now")

        Re1 = self._PETs[0]
        Re2 = self._PETs[1]

        Re1_inv = Re1.inverse()
        Re2_inv = Re2.inverse()

        P = self._partition
        key_to_column = {k:[k] for (k,atoms) in P}
        #print("k",key_to_column)

        for j in range(sizes[1]-1):
            P,d = self._partition.refinement(Re2_inv(P), certificate=True)
            #print("d",d)
            key_to_column = {k:[d[k][0]]+key_to_column[d[k][1]] for k in d}
            #print("k",key_to_column)

        Q = P
        key_to_word = {k:[col] for (k,col) in key_to_column.items()}
        for i in range(sizes[0]-1):
            Q,d = P.refinement(Re1_inv(Q), certificate=True)
            key_to_word = {k:[key_to_column[d[k][0]]]+key_to_word[d[k][1]] for k in d}

        return Q, key_to_word



