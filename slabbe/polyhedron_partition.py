# -*- coding: utf-8 -*-
r"""
Polyhedron partitions

EXAMPLES:

A polyhedron partition::

    sage: from slabbe import PolyhedronPartition
    sage: h = 1/3
    sage: p = Polyhedron([(0,h),(0,1),(h,1)])
    sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
    sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
    sage: s = Polyhedron([(h,0), (1,0), (1,h)])
    sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
    sage: P.is_pairwise_disjoint()
    True
    sage: list(P)
    [(0, A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices),
     (1, A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices),
     (2, A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices),
     (3, A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices)]
    sage: G = P.plot()

AUTHORS:

- Sébastien Labbé, November 2017, initial version of polyhedron partitions
"""
#*****************************************************************************
#       Copyright (C) 2017-2019 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
import itertools
from sage.misc.cachefunc import cached_method
from sage.geometry.polyhedron.constructor import Polyhedron

def find_unused_key(d, sequence):
    r"""
    Return the first key in sequence which is not in d.
    
    EXAMPLES::

        sage: from slabbe.polyhedron_partition import find_unused_key
        sage: d = {3:32, 0:21, 1:4, 5:5}
        sage: find_unused_key(d, NN)
        2
        sage: d[2] = 1234
        sage: find_unused_key(d, NN)
        4
        sage: d[4] = 1234
        sage: find_unused_key(d, NN)
        6
    """
    for a in sequence:
        if a not in d:
            return a

def is_union_convex(t):
    r"""
    Return whether the union of the polyhedrons is convex.

    INPUT:

    - ``t`` -- list of polyhedron

    EXAMPLES::

        sage: from slabbe.polyhedron_partition import is_union_convex
        sage: h = 1/2
        sage: p = Polyhedron([(0,h),(0,1),(h,1)])
        sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
        sage: r = Polyhedron([(h,0), (1,0), (1,h)])
        sage: is_union_convex((p,q,r))
        True
        sage: is_union_convex((p,q))
        True
        sage: is_union_convex((p,r))
        False

    Here we need to consider the three at the same time to get a convex
    union::

        sage: h = 1/5
        sage: p = Polyhedron([(0,0),(h,1-h),(0,1)])
        sage: q = Polyhedron([(0,1), (h,1-h), (1,1)])
        sage: r = Polyhedron([(0,0), (h,1-h), (1,1), (1,0)])
        sage: is_union_convex((p,q))
        False
        sage: is_union_convex((p,r))
        False
        sage: is_union_convex((q,r))
        False
        sage: is_union_convex((p,q,r))
        True
    """
    if not t:
        return True
    base_ring = t[0].base_ring()
    vertices = sum((p.vertices() for p in t), tuple())
    r = Polyhedron(vertices, base_ring=base_ring)
    return r.volume() == sum(p.volume() for p in t)

class PolyhedronPartition(object):
    r"""
    Return a partition into polyhedron.

    Note: Many atoms may share the same key.

    INPUT:

    - ``atoms`` -- list of polyhedron or dict of key -> polyhedron or list
      of (key, polyhedron)
    - ``base_ring`` -- base ring (default: ``None``) of the vertices

    EXAMPLES::

        sage: from slabbe import PolyhedronPartition
        sage: h = 1/2
        sage: p = Polyhedron([(0,h),(0,1),(h,1)])
        sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
        sage: r = Polyhedron([(h,0), (1,0), (1,h)])
        sage: P = PolyhedronPartition([p,q,r])
        sage: P
        Polyhedron partition of 3 atoms with 3 letters

    ::

        sage: P.is_pairwise_disjoint()
        True
        sage: P.volume()
        1
        sage: G = P.plot()

    From a dict::

        sage: PolyhedronPartition(dict(a=p,b=q,c=r))
        Polyhedron partition of 3 atoms with 3 letters

    From a list of (key, polyhedron)::

        sage: PolyhedronPartition([(9,p),(8,q),(9,r)])
        Polyhedron partition of 3 atoms with 2 letters
    """
    def __init__(self, atoms, base_ring=None):
        r"""
        See class for documentation.

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
        """
        if isinstance(atoms, list):
            if all(hasattr(p, 'vertices') for p in atoms):
                self._items = list(enumerate(atoms))
            elif all(isinstance(t, tuple) for t in atoms):
                self._items = atoms
            else:
                raise TypeError('atoms (={}) must be a list of polyhedron or a'
                        ' list of tuples'.format(atoms))
        elif isinstance(atoms, dict):
            self._items = atoms.items()
        else:
            raise TypeError('atoms (={}) must be a list or a'
                    ' dict'.format(atoms))
        if base_ring is None:
            if len(self) == 0:
                from sage.rings.integer_ring import ZZ
                base_ring == ZZ
            else:
                base_ring = next(iter(self))[1].base_ring()
        self._base_ring = base_ring

    @classmethod
    def jeandel_rao_tilings_partition(cls, backend=None):
        r"""
        This construct the polygon partition associated to Jeandel-Rao
        tilings introduced in [Lab2019]_.

        INPUT:

        - ``backend`` -- string, polyhedron backend

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: P0 = PolyhedronPartition.jeandel_rao_tilings_partition()
            sage: P0.is_pairwise_disjoint()
            True
            sage: P0.volume()
            4*phi + 1

        The volume is consistent with::

            sage: z = polygen(QQ, 'z')
            sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: phi = K.gen()
            sage: phi * (phi + 3)
            4*phi + 1

        REFERENCES:

        .. [Lab2019] S. Labbé. A Markov partition for Jeandel-Rao aperiodic
           Wang tilings. March 2019. https://arxiv.org/abs/1903.06137
        """
        from sage.rings.polynomial.polynomial_ring import polygen
        from sage.rings.rational_field import QQ
        from sage.rings.number_field.number_field import NumberField

        # the golden mean
        z = polygen(QQ, 'z')
        K = NumberField(z**2-z-1, 'phi', embedding=QQ(1.6))
        phi = K.gen()

        # x and y coordinates
        xcoords = [0, phi**-2, phi**-1, 1, phi]
        ycoords = [0, 1, 2, phi+1, phi+2, phi+3]

        # the vertices
        A = [(x,0) for x in xcoords]
        B = [(x,1) for x in xcoords]
        C = [(x,2) for x in xcoords]
        D = [(x,phi+1) for x in xcoords]
        E = [(x,phi+2) for x in xcoords]
        F = [(x,phi+3) for x in xcoords]

        # atoms corresponding to Jeandel-Rao tile numbers
        L = [
            (0, (A[0], A[2], B[2])), 
            (0, (A[2], A[3], B[3])),
            (0, (A[3], A[4], B[4])),
            (1, (A[0], B[0], B[2])),
            (1, (A[2], B[2], B[3])),
            (1, (A[3], B[3], B[4])),
            (2, (C[0], D[0], E[2])),
            (3, (B[3], C[3], E[4], C[4])),
            (4, (C[0], E[2], E[3])),
            (4, (E[2], E[3], F[3])),
            (5, (D[0], E[0], E[2])),
            (5, (E[0], E[1], F[1])),
            (5, (E[1], E[2], F[3])),
            (6, (E[0], F[0], F[1])),
            (6, (E[1], F[1], F[3])),
            (6, (E[3], F[3], F[4])),
            (7, (D[3], E[3], E[4])),
            (7, (E[3], E[4], F[4])),
            (7, (B[0], C[0], E[3])),
            (8, (B[2], C[2], E[4])),
            (9, (B[0], B[2], C[2])),
            (9, (B[2], B[3], C[3])),
            (9, (B[3], B[4], C[4])),
            (10, (B[0], D[3], E[3])),
            ]
        L = [(key, Polyhedron(vertices, base_ring=K, backend=backend)) 
             for (key,vertices) in L]
        return PolyhedronPartition(L)

    @classmethod
    def self_similar_19_tiles_partition(cls):
        r"""
        This construct the polygon partition introduced in [Lab2019]_
        associated to the self-similar 19 Wang tiles [Lab2018]_.

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: PU = PolyhedronPartition.self_similar_19_tiles_partition()
            sage: PU.is_pairwise_disjoint()
            True
            sage: PU.volume()
            1

        REFERENCES:

        .. [Lab2018] S. Labbé. A self-similar aperiodic set of 19 Wang
           tiles. Geom. Dedicata, 2018.
           https://doi.org/10.1007/s10711-018-0384-8.
        """
        from sage.rings.polynomial.polynomial_ring import polygen
        from sage.rings.rational_field import QQ
        from sage.rings.number_field.number_field import NumberField

        # the golden mean
        z = polygen(QQ, 'z')
        K = NumberField(z**2-z-1, 'phi', embedding=QQ(1.6))
        phi = K.gen()

        # the partition vertices
        L = [
         (0, [(phi - 1, -2*phi + 4), (phi - 1, phi - 1), (-2*phi + 4, phi - 1)]),
         (1,
          [(phi - 1, -2*phi + 4),
           (phi - 1, 1),
           (1, 1),
           (-2*phi + 4, phi - 1),
           (1, phi - 1)]),
         (2, [(0, 1), (-phi + 2, phi - 1), (2*phi - 3, phi - 1)]),
         (3, [(0, 1), (0, phi - 1), (2*phi - 3, phi - 1)]),
         (4, [(phi - 1, phi - 1), (-phi + 2, 1), (phi - 1, -2*phi + 4)]),
         (5, [(phi - 1, 1), (-phi + 2, 1), (phi - 1, -2*phi + 4)]),
         (6, [(0, 1), (-phi + 2, 1), (-phi + 2, phi - 1)]),
         (7, [(-phi + 2, 1), (-phi + 2, phi - 1), (phi - 1, phi - 1)]),
         (8, [(1, 0), (phi - 1, 0), (phi - 1, -phi + 2)]),
         (9, [(-2*phi + 4, phi - 1), (phi - 1, phi - 1), (1, -phi + 2), (1, 0)]),
         (10, [(1, -phi + 2), (-2*phi + 4, phi - 1), (1, phi - 1)]),
         (11, [(1, 0), (phi - 1, -phi + 2), (phi - 1, phi - 1)]),
         (12, [(-phi + 2, -phi + 2), (-phi + 2, phi - 1), (2*phi - 3, phi - 1)]),
         (13,
          [(-phi + 2, phi - 1),
           (phi - 1, -phi + 2),
           (-phi + 2, -phi + 2),
           (phi - 1, 0)]),
         (14,
          [(2*phi - 3, phi - 1), (0, phi - 1), (-phi + 2, -phi + 2), (-phi + 2, 0)]),
         (15, [(phi - 1, 0), (-phi + 2, 0), (-phi + 2, -phi + 2)]),
         (16, [(0, 0), (-phi + 2, 0), (0, -phi + 2)]),
         (17, [(0, -phi + 2), (-phi + 2, 0), (0, phi - 1)]),
         (18, [(phi - 1, -phi + 2), (-phi + 2, phi - 1), (phi - 1, phi - 1)])
         ]

        L = [(key, Polyhedron(vertices, base_ring=K)) for (key,vertices) in L]
        return PolyhedronPartition(L)

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: next(iter(P))
            (0,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull
             of 3 vertices)
        """
        return iter(self._items)

    def atoms(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: P.atoms()
            [A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 6 vertices,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices]
        """
        return [atom for (key,atom) in self._items]

    def base_ring(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: P.base_ring()
            Rational Field
        """
        return self._base_ring

    def ambient_space(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: P.ambient_space()
            Vector space of dimension 2 over Rational Field
        """
        return next(iter(self))[1].ambient_space()

    @cached_method
    def cached_atoms_set(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: P.cached_atoms_set()
            {A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices,
             A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 6 vertices}
        """
        return set(self.atoms())

    def __eq__(self, other):
        r"""
        Return whether two partitions are the same.

        The coding is not considered.

        INPUT:

        - ``other`` -- a partition

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: Q = PolyhedronPartition([p,q])
            sage: R = PolyhedronPartition({'asd':q, 'yo':r, 'foo':p})
            sage: P == P
            True
            sage: P == Q
            False
            sage: P == R
            True
        """
        return (isinstance(other, PolyhedronPartition) and
                self.cached_atoms_set() == other.cached_atoms_set())

    def __contains__(self, p):
        r"""
        Return whether a polyhedron is an atom of the partition.

        INPUT:

        - ``p`` -- a polyhedron

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: p in P
            True

        ::

            sage: Q = PolyhedronPartition([p,q])
            sage: r in Q
            False
        """
        return p in self.cached_atoms_set()

    def __len__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: len(P)
            3
        """
        return len(self._items)

    def __getitem__(self, key):
        r"""
        Return the list of atoms associated with the given key.

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([(0,p),(1,q),(0,r)])
            sage: P[0]
            Polyhedron partition of 2 atoms with 1 letters
            sage: P[1]
            Polyhedron partition of 1 atoms with 1 letters

        """
        return PolyhedronPartition([(k,atom) for (k,atom) in self if k == key])

    def alphabet(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([(3,p), (5,q), (9,r)])
            sage: P.alphabet()
            {3, 5, 9}
            sage: P = PolyhedronPartition([(3,p), (5,q), (3,r)])
            sage: P.alphabet()
            {3, 5}
        """
        return set(key for (key,atom) in self)

    def alphabet_size(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([(3,p), (5,q), (9,r)])
            sage: P.alphabet_size()
            3
            sage: P = PolyhedronPartition([(3,p), (5,q), (3,r)])
            sage: P.alphabet_size()
            2
        """
        return len(self.alphabet())

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: PolyhedronPartition([p,q,r])
            Polyhedron partition of 3 atoms with 3 letters
        """
        return ("Polyhedron partition of {} atoms "
                "with {} letters".format(len(self),self.alphabet_size()))

    def __rmul__(self, factor):
        r"""
        Returns the partition of the induced transformation on the domain.

        INPUT:

        - ``factor`` -- number

        OUTPUT:

            a polyhedron partition

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: 4 * P
            Polyhedron partition of 3 atoms with 3 letters
            sage: -4 * P
            Polyhedron partition of 3 atoms with 3 letters

        TESTS::

            sage: (4.5 * P).base_ring()
            Real Double Field
        """
        return PolyhedronPartition({key:factor*p for (key,p) in self})

    def __neg__(self):
        r"""
        Returns the miror image of the partition.

        OUTPUT:

            a polyhedron partition

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: -P
            Polyhedron partition of 3 atoms with 3 letters
        """
        return PolyhedronPartition({key:-p for (key,p) in self})

    def rename_keys(self, d):
        r"""
        Return a polyhedron partition whose keys are the images under a map.
        
        INPUT:

        - ``d`` -- dict, function old key -> new key

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: Q = P.rename_keys({0:'b', 1:'a', 2:'z'})
            sage: Q
            Polyhedron partition of 3 atoms with 3 letters
            sage: sorted(key for key,p in Q)
            ['a', 'b', 'z']

        It does not have to be injective::

            sage: Q = P.rename_keys({0:'b', 1:'a', 2:'b'})
            sage: sorted(key for key,p in Q)
            ['a', 'b', 'b']
        """
        return PolyhedronPartition([(d[key],p) for (key,p) in self])

    def keys_permutation(self, other):
        r"""
        Return a relabelling permutation of the keys for self to look like
        other.

        .. NOTE::

            currently, the code works only if the coding of self and other
            is injective, i.e., no two polyhedron are coded by the same
            letter.

        INPUT:

        - ``other`` -- a polyhedron partition (with injective coding)

        OUTPUT:

            dict, key -> key

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({4:p, 1:q, 2:r})
            sage: Q = PolyhedronPartition({0:p, 5:q})
            sage: d = P.keys_permutation(Q)
            sage: d
            {1: 5, 2: 1, 4: 0}
            sage: P.rename_keys(d)
            Polyhedron partition of 3 atoms with 3 letters
        """
        if not isinstance(other, PolyhedronPartition):
            raise TypeError("other (of type={}) must a polyhedron"
                    " partition".format(type(other)))

        if self.alphabet_size() != len(self):
            raise NotImplementedError('keys_permutation method is'
                ' implemented only if the coding of self (={}) is'
                ' injective'.format(self))
        if other.alphabet_size() != len(other):
            raise NotImplementedError('keys_permutation method is'
                ' implemented only if the coding of self (={}) is'
                ' injective'.format(self))

        d = {}
        atoms_not_in_other = []
        for self_key,p in self:
            if p in other:
                d[self_key] = other.code(p)
            else:
                atoms_not_in_other.append((self_key,p))
        forbidden_keys = set(key for (key,q) in other._items)
        for self_key,p in atoms_not_in_other:
            new_key = find_unused_key(forbidden_keys, itertools.count())
            d[self_key] = new_key
            forbidden_keys.add(new_key)
        return d

    def keys_permutation_lexicographic(self):
        r"""
        Return a permutation relabelling of the keys for self in increasing
        order for the lexicographic order of the centers of the polyhedrons.

        OUTPUT:

            dict, key -> key

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({4:p, 1:q, 2:r})
            sage: d = P.keys_permutation_lexicographic()
            sage: d
            {1: 1, 2: 2, 4: 0}
            sage: P.rename_keys(d)
            Polyhedron partition of 3 atoms with 3 letters

        ::

            sage: Q = PolyhedronPartition({0:p, 5:q})
            sage: Q.keys_permutation_lexicographic()
            {0: 0, 5: 1}

        It works when the partition has two atoms coded by the same key::

            sage: P = PolyhedronPartition([(0,p), (0,q), (3,r)])
            sage: d = P.keys_permutation_lexicographic()
            sage: d
            {0: 0, 3: 1}
            sage: P.rename_keys(d).alphabet()
            {0, 1}

        """
        L = [(key, p.center()) for (key,p) in self]
        L.sort(key=lambda t:t[1])
        d = {}
        it = itertools.count()
        for key,center in L:
            if key not in d:
                d[key] = next(it)
        return d

    def apply_linear_map(self, M):
        r"""
        INPUT:

        - ``M`` -- a matrix

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})

        Vertical symmetry::

            sage: M = diagonal_matrix((-1,1))
            sage: P = P.apply_linear_map(M)
            sage: P = P.translation((1,0))
            sage: P
            Polyhedron partition of 4 atoms with 4 letters

        """
        L = [(key, M*p) for key,p in self]
        return PolyhedronPartition(L)

    def translation(self, displacement):
        """
        Return the translated partition of polyhedron.

        INPUT:

        - ``displacement`` -- a displacement vector or a list/tuple of
          coordinates that determines a displacement vector.

        OUTPUT:

        The translated partition.

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: P.translation((1,1))
            Polyhedron partition of 3 atoms with 3 letters
        """
        return PolyhedronPartition([(key,p.translation(displacement)) for (key,p) in self])

    def volume(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: P.volume()
            1

        TESTS::

            sage: PolyhedronPartition([], base_ring=ZZ).volume()
            0
        """
        return sum(p.volume() for p in self.atoms())

    def volume_dict(self, normalize=False):
        r"""
        INPUT

        - ``normalize`` -- boolean (default:``False``), whether to
          normalize the sum of the whole volume to 1

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: P.volume_dict()
            {0: 1/8, 1: 3/4, 2: 1/8}
            sage: (2*P).volume_dict()
            {0: 1/2, 1: 3, 2: 1/2}
        """
        from collections import Counter
        d = Counter()
        for (key,atom) in self._items:
            d[key] += atom.volume()
        if normalize:
            volume = self.volume()
            return {key:v/volume for (key,v) in d.items()}
        else:
            return dict(d)

    def plot(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: P.plot()
            Graphics object consisting of 21 graphics primitives
        """
        from sage.plot.graphics import Graphics
        from sage.plot.text import text
        G = Graphics()
        for key,P in self:
            G += P.plot(fill='white')
            G += text(key, P.center())
        return G

    def edges(self):
        r"""
        Return the edges of partition (one copy of each edge).

        .. NOTE::

            If there are vertices of atoms on the interior of the edge of
            another atom, then, the overlapping edges will be repeated.

        OUTPUT:

        - set of sorted pair of immutable vectors

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: sorted(P.edges())
            [((0, 0), (0, 1/2)),
             ((0, 0), (1/2, 0)),
             ((0, 1/2), (0, 1)),
             ((0, 1/2), (1/2, 1)),
             ((0, 1), (1/2, 1)),
             ((1/2, 0), (1, 0)),
             ((1/2, 0), (1, 1/2)),
             ((1/2, 1), (1, 1)),
             ((1, 0), (1, 1/2)),
             ((1, 1/2), (1, 1))]

        Irrational partition::

            sage: z = polygen(QQ, 'z') #z = QQ['z'].0 # same as
            sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: phi = K.gen()
            sage: h = 1/phi^2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s}, base_ring=K)
            sage: sorted(P.edges())
            [((0, 0), (0, -phi + 2)),
             ((0, 0), (-phi + 2, 0)),
             ((0, -phi + 2), (0, 1)),
             ((0, -phi + 2), (-phi + 2, 1)),
             ((0, 1), (-phi + 2, 1)),
             ((-phi + 2, 0), (-phi + 2, 1)),
             ((-phi + 2, 0), (1, 0)),
             ((-phi + 2, 0), (1, -phi + 2)),
             ((-phi + 2, 1), (1, 1)),
             ((1, 0), (1, -phi + 2)),
             ((1, -phi + 2), (1, 1))]
        """
        edges = set()
        for key,P in self:
            proj = P.projection()
            for (a,b) in proj.lines:
                a = proj.coords[a]
                b = proj.coords[b]
                a.set_immutable()
                b.set_immutable()
                sorted_edge = tuple(sorted((a,b)))
                edges.add(sorted_edge)
        return edges

    def tikz(self, fontsize=r'\normalsize', scale=1, 
            label_format = r'{}',
            extra_code=''):
        r"""
        INPUT:

        - ``fontsize`` -- string (default: ``r'\normalsize'``
        - ``scale`` -- number (default: ``1``)
        - ``label_format`` -- string (default: ``r'{}'``) to be called with
          ``label_format.format(key)``
        - ``extra_code`` -- string (default: ``''``)

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: _ = P.tikz().pdf(view=False)

        Irrational partition::

            sage: z = polygen(QQ, 'z') #z = QQ['z'].0 # same as
            sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: phi = K.gen()
            sage: h = 1/phi^2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s}, base_ring=K)
            sage: _ = P.tikz().pdf(view=False)

        Testing the options::

            sage: _ = P.tikz(fontsize=r'\scriptsize').pdf(view=False)
            sage: _ = P.tikz(scale=2).pdf(view=False)
            sage: _ = P.tikz(label_format=r'$a_{{{}}}$').pdf(view=False)
        """
        from slabbe import TikzPicture
        lines = []
        lines.append(r'\begin{tikzpicture}')
        lines.append('[scale={}]'.format(scale))
        # edges
        for (a,b) in self.edges():
            a = a.n(digits=5)
            b = b.n(digits=5)
            line = r'\draw {} -- {};'.format(a,b)
            lines.append(line)
        # node key
        node_format = r'\node[font={}] at {} {{{}}};'
        for key,P in self:
            lines.append(r'% atom with key {}'.format(key))
            label = label_format.format(key)
            lines.append(node_format.format(fontsize, 
                                            P.center().n(digits=5),
                                            label))
        lines.append(extra_code)
        lines.append(r'\end{tikzpicture}')
        return TikzPicture('\n'.join(lines))

    def is_pairwise_disjoint(self):
        r"""
        Return whether atoms of the partition are pairwise disjoint.

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])
            sage: P.is_pairwise_disjoint()
            True
        """
        for ((keyi,p),(keyj,q)) in itertools.permutations(self, 2):
            volume = p.intersection(q).volume()
            if volume > 0:
                raise ValueError('Intersection of atom with key {} (={})'
                    ' and atom with key {} (={}) is not of zero volume'
                    ' (={})'.format(keyi, p.vertices(), keyj, q.vertices(),
                        volume))
        return True

    def merge_atoms(self, d):
        r"""
        Return the polyhedron partition obtained by merging atoms having
        the same image under the dictionnary.

        INPUT:

        - ``d`` -- dict

        OUTPUT:

            a polyhedron partition

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r})
            sage: P.merge_atoms({0:4, 1:4, 2:5})
            Polyhedron partition of 2 atoms with 2 letters
            sage: P.merge_atoms({0:4, 1:5, 2:4})
            Polyhedron partition of 3 atoms with 2 letters

        When pair of atoms are not convex, it needs to merge 3 or more
        atoms::

            sage: h = 1/5
            sage: p = Polyhedron([(0,0),(h,1-h),(0,1)])
            sage: q = Polyhedron([(0,1), (h,1-h), (1,1)])
            sage: r = Polyhedron([(0,0), (h,1-h), (1,1), (1,0)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r})
            sage: P.merge_atoms({0:4, 1:4, 2:4})
            Polyhedron partition of 1 atoms with 1 letters
        """
        from collections import defaultdict
        from sage.misc.misc import exists

        to_merge = defaultdict(set)
        for key,p in self:
            to_merge[d[key]].add(p)

        final_atoms = []
        for new_key,atoms in to_merge.items():
            subsets = (s for k in range(2,len(atoms)+1)
                         for s in itertools.permutations(atoms, k))
            answer,t = exists(subsets, is_union_convex)
            while answer:
                vertices = sum((p.vertices() for p in t), tuple())
                base_ring = t[0].base_ring()
                r = Polyhedron(vertices, base_ring=base_ring)
                for p in t:
                    atoms.remove(p)
                atoms.add(r)
                subsets = (s for k in range(2,len(atoms)+1)
                             for s in itertools.permutations(atoms, k))
                answer,t = exists(subsets, is_union_convex)
            for atom in atoms:
                final_atoms.append((new_key, atom))

        return PolyhedronPartition(final_atoms)

    def domain(self):
        r"""
        Return the domain of the partition.

        OUTPUT:

            a polyhedron

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/3
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: P.domain()
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
            sage: P.domain().vertices()
            (A vertex at (0, 0),
             A vertex at (0, 1),
             A vertex at (1, 0),
             A vertex at (1, 1))
        """
        union = self.merge_atoms({a:0 for a in self.alphabet()})
        if not len(union) == 1:
            raise NotImplementedError("non convex domain (={})".format(union))
        [(key,polyhedron)] = union
        return polyhedron

    def code(self, p):
        r"""
        Returns in which atom the polyhedron lives in.

        INPUT:

        - ``p`` -- a polyhedron

        OUTPUT:

            integer (for the i-th atom)

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: P.code(p)
            0
            sage: P.code(q)
            1
            sage: t = Polyhedron([(0, 8/9), (0, 1), (1/9, 1)])
            sage: P.code(t)
            0

        TESTS::

            sage: t = Polyhedron([(0, 1/9), (0, 1), (1/9, 1)])
            sage: P.code(t)
            Traceback (most recent call last):
            ...
            ValueError: polyhedron p whose vertices are (A vertex at (0,
            1), A vertex at (0, 1/9), A vertex at (1/9, 1)) lies in no atom
        """
        if not hasattr(p, 'vertices'):
            raise TypeError('p (={}) must be a polyhedron'.format(p))
        L = [i for i,atom in self if p <= atom]
        if len(L) == 1:
            return L[0]
        elif len(L) > 1:
            raise ValueError("polyhedron p whose vertices are {} lies "
                    "in more than one atoms (={})".format(p.vertices(), L))
        else:
            raise ValueError("polyhedron p whose vertices are {} lies "
                    "in no atom".format(p.vertices()))


    def refine_by_hyperplane(self, ieq):
        r"""
        Refine the partition with the two half spaces of each side of an
        hyperplane.

        INPUT:

        - ``ieq`` -- list, an inequality. An entry equal to "[-1,7,3,4]"
          represents the inequality 7x_1+3x_2+4x_3>= 1.

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq = [-4, 5, 1]
            sage: P.refine_by_hyperplane(ieq)
            Polyhedron partition of 6 atoms with 6 letters
        """
        half = Polyhedron(ieqs=[ieq])
        half_partition = PolyhedronPartition([half])
        other_half = Polyhedron(ieqs=[[-a for a in ieq]])
        other_half_partition = PolyhedronPartition([other_half])
        A = self.refinement(half_partition)
        B = self.refinement(other_half_partition)
        return PolyhedronPartition(A.atoms()+B.atoms())

    def refinement(self, other, key_fn=None):
        r"""
        Return the polyhedron partition obtained by the intersection of the
        atoms of self with the atoms of other.

        Only atoms of positive volume are kept.

        INPUT:

        - ``other`` -- a polyhedron partition
        - ``key_fn`` -- function to apply on pairs of labels, or None

        OUTPUT:

            a polyhedron partition

        .. TODO::

            Avoid a quadratic enumeration of all pairs.

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: g = 1/5
            sage: t1 = Polyhedron([(g,g), (g,1-g), (1-g,g) ])
            sage: t2 = Polyhedron([(g,1-g), (1-g,g), (1-g,1-g)])
            sage: Q = PolyhedronPartition([t1,t2])
            sage: P.refinement(Q)
            Polyhedron partition of 8 atoms with 8 letters
        """
        if not isinstance(other, PolyhedronPartition):
            raise TypeError("other (of type={}) must a polyhedron"
                    " partition".format(type(other)))
        L = []
        for ((a,p),(b,q)) in itertools.product(self, other):
            p_q = p.intersection(q)
            if p_q.is_full_dimensional():
                if key_fn is None:
                    L.append(p_q)
                else:
                    new_key = key_fn(a,b)
                    L.append((new_key, p_q))
        return PolyhedronPartition(L)
        
