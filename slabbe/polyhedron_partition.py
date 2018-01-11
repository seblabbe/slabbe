# -*- coding: utf-8 -*-
r"""
Polyhedron partition and induction

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

Applying a rationnal rotation::

    sage: from slabbe import rotation_mod
    sage: u = rotation_mod(0, 2/3, 1, QQ)
    sage: Q = P.apply_transformation(u)
    sage: Q
    Polyhedron partition of 4 atoms with 4 letters

Inducing an irrationnal rotation on a subdomain::

    sage: z = polygen(QQ, 'z') #z = QQ['z'].0 # same as
    sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
    sage: phi = K.gen()
    sage: h = 1/phi^2
    sage: p = Polyhedron([(0,h),(0,1),(h,1)])
    sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
    sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
    sage: s = Polyhedron([(h,0), (1,0), (1,h)])
    sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s}, base_ring=K)
    sage: u = rotation_mod(0, 1/phi, 1, K)
    sage: u_inv = rotation_mod(0, 1/phi^2, 1, K)
    sage: ieq = [h, -1, 0]   # x0 <= h
    sage: P1,sub01 = P.induced_partition(u, u_inv, ieq)
    sage: P1
    Polyhedron partition of 7 atoms with 7 letters
    sage: sub01
    {0: [0, 2],
     1: [0, 2, 2],
     2: [1, 2],
     3: [1, 2, 2],
     4: [1, 3],
     5: [1, 3, 2],
     6: [1, 3, 3]}

AUTHORS:

- Sébastien Labbé, November 2017, initial version
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
from __future__ import absolute_import, print_function
import itertools
from copy import copy
from sage.misc.cachefunc import cached_method
from sage.geometry.polyhedron.constructor import Polyhedron

def rotation_mod(i, angle, mod, base_ring):
    r"""
    Return a rotation function acting on polyhedron.

    INPUT:

    - ``i`` -- integer, coordinate of the rotation
    - ``angle`` -- number, angle of rotation
    - ``mod`` -- number, modulo 
    - ``base_ring`` -- ring, base ring for the vertices of the polyhedron

    OUTPUT:

        a function defined on polyhedron

    EXAMPLES::

        sage: from slabbe import rotation_mod
        sage: z = polygen(QQ, 'z') #z = QQ['z'].0 # same as
        sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
        sage: phi = K.gen()
        sage: p = ((-phi + 2, phi - 1), (-phi + 2, 1), (phi - 1, 1))
        sage: p = Polyhedron(p, base_ring=K)

    A rotation modulo phi on the x coordinate::

        sage: t0 = rotation_mod(0, 1, phi, K)
        sage: t0(p).vertices()
        (A vertex at (-phi + 3, phi - 1),
         A vertex at (-phi + 3, 1),
         A vertex at (phi, 1))

    The inverse map::

        sage: t0_inv = rotation_mod(0, 1/phi, phi, K)
        sage: t0(p) == p
        False
        sage: t0_inv(t0(p)) == p
        True

    A rotation modulo 1 on the y coordinate::

        sage: t1 = rotation_mod(1, 1/phi^2, 1, K)
        sage: t1(p).vertices()
        (A vertex at (-phi + 2, 0),
         A vertex at (-phi + 2, -phi + 2),
         A vertex at (phi - 1, -phi + 2))
    """
    def trans(p):
        if all(v[i] <= mod-angle for v in p.vertices()):
            L = [tuple(vj+angle if j==i else vj 
                    for (j,vj) in enumerate(v))
                    for v in p.vertices()]
        else:
            L = [tuple(vj+angle-mod if j==i else vj 
                    for (j,vj) in enumerate(v))
                    for v in p.vertices()]
        return Polyhedron(L, base_ring=base_ring)
    return trans

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

        .. NOTE::

            currently, the code works only if the coding of self and other
            is injective, i.e., no two polyhedron are coded by the same
            letter.

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
            sage: Q = PolyhedronPartition({0:p, 5:q})
            sage: Q.keys_permutation_lexicographic()
            {0: 0, 5: 1}
        """
        if self.alphabet_size() != len(self):
            raise NotImplementedError('keys_permutation method is'
                ' implemented only if the coding of self (={}) is'
                ' injective'.format(self))

        L = [(key, p.center()) for (key,p) in self]
        L.sort(key=lambda t:t[1])
        d = {key:i for i,(key,center) in enumerate(L)}
        return d

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

    def tikz(self, fontsize=r'\normalsize', scale=1):
        r"""
        INPUT:

        - ``fontsize`` -- string (default: ``r'\normalsize'``
        - ``scale`` -- number (default: ``1``)

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
        for key,P in self:
            lines.append(r'% atom with key {}'.format(key))
            node_str = r'\node[font={}] at {} {{{}}};'
            lines.append(node_str.format(fontsize, 
                                         P.center().n(digits=5),
                                         key))
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

    def apply_transformation(self, trans):
        r"""
        INPUT:

        - ``trans`` -- a function: polyhedron -> polyhedron

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, rotation_mod
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: u = rotation_mod(0, 2/3, 1, QQ)
            sage: Q = P.apply_transformation(u)
            sage: Q
            Polyhedron partition of 4 atoms with 4 letters

        ::

            sage: u = rotation_mod(0, 2/3, 1, QQ)
            sage: u_inv = rotation_mod(0, 1/3, 1, QQ)
            sage: R = P.apply_transformation(u).apply_transformation(u_inv)
            sage: P == R
            True
        """
        L = []
        for key,p in self:
            trans_p = trans(p)
            assert p.volume() == trans_p.volume()
            L.append((key, trans_p))
        return PolyhedronPartition(L)

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

    def refinement(self, other):
        r"""
        Return the polyhedron partition obtained by the intersection of the
        atoms of self with the atoms of other.

        Only atoms of positive volume are kept.

        INPUT:

        - ``other`` -- a polyhedron partition

        OUTPUT:

            a polyhedron partition

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, rotation_mod
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
        for (p,q) in itertools.product(self.atoms(), other.atoms()):
            p_q = p.intersection(q)
            if p_q.volume() > 0:
                L.append(p_q)
        return PolyhedronPartition(L)
        
    def induced_out_partition(self, trans, ieq):
        r"""
        Returns the output partition obtained as the induction of the given
        transformation on the domain given by an inequality.

        Note: the output partition corresponds to the arrival partition in
        the domain, not the initial one.

        INPUT:

        - ``trans`` -- a function: polyhedron -> polyhedron
        - ``ieq`` -- list, an inequality. An entry equal to "[-1,7,3,4]"
          represents the inequality 7x_1+3x_2+4x_3>= 1.

        OUTPUT:

            dict of polyhedron partitions with keys giving the return time

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, rotation_mod
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: u = rotation_mod(0, 1/3, 1, QQ)
            sage: ieq = [h, -1, 0]   # x0 <= h
            sage: P.induced_out_partition(u, ieq)
            {3: Polyhedron partition of 4 atoms with 4 letters}

        ::

            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq2 = [1/2, -1, 0]   # x0 <= 1/2
            sage: d = P.induced_out_partition(u, ieq2)
            sage: d
            {1: Polyhedron partition of 2 atoms with 2 letters,
             2: Polyhedron partition of 3 atoms with 3 letters,
             3: Polyhedron partition of 4 atoms with 4 letters}
            sage: Q = PolyhedronPartition(d[1].atoms()+d[2].atoms()+d[3].atoms())
            sage: Q.is_pairwise_disjoint()
            True

        ::

            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq3 = [-1/2, 1, 0]   # x0 >= 1/2
            sage: P.induced_out_partition(u, ieq3)
            {2: Polyhedron partition of 3 atoms with 3 letters,
             3: Polyhedron partition of 4 atoms with 4 letters}

        It is an error if the induced region is empty::

            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq4 = [-1/2, -1, 0]   # x0 <= -1/2
            sage: P.induced_out_partition(u, ieq4)
            Traceback (most recent call last):
            ...
            ValueError: Inequality An inequality (-2, 0) x - 1 >= 0 does
            not intersect P (=Polyhedron partition of 4 atoms with 4
            letters)

        The whole domain::

            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq5 = [1/2, 1, 0]   # x0 >= -1/2
            sage: P.induced_out_partition(u, ieq5)
            {1: Polyhedron partition of 4 atoms with 4 letters}

        An irrational rotation::

            sage: z = polygen(QQ, 'z') #z = QQ['z'].0 # same as
            sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: phi = K.gen()
            sage: h = 1/phi^2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s}, base_ring=K)
            sage: u = rotation_mod(0, 1/phi, 1, K)
            sage: ieq = [phi^-4, -1, 0]   # x0 <= phi^-4
            sage: d = P.induced_out_partition(u, ieq)
            sage: d
            {5: Polyhedron partition of 6 atoms with 6 letters,
             8: Polyhedron partition of 9 atoms with 9 letters}
        """
        # good side of the hyperplane
        half = Polyhedron(ieqs=[ieq])
        half_part = PolyhedronPartition([half])
        # the other side of the hyperplane
        other_half = Polyhedron(ieqs=[[-a for a in ieq]])
        other_half_part = PolyhedronPartition([other_half])
        # initial refinement
        P = self.refinement(half_part)
        if len(P) == 0:
            raise ValueError("Inequality {} does not intersect P "
                    "(={})".format(half.inequalities()[0], self))
        level = 1
        ans = {}
        P = P.apply_transformation(trans)
        while len(P):
            P_returned = P.refinement(half_part)
            if P_returned:
                ans[level] = P_returned
            # for what is remaining we do:
            P = P.refinement(other_half_part)
            P = P.refinement(self)
            P = P.apply_transformation(trans)
            level += 1
        return ans

    def induced_in_partition(self, trans, trans_inv, ieq):
        r"""
        Returns the partition of the induced transformation on the domain.
        given by an inequality.

        INPUT:

        - ``trans`` -- a function: polyhedron -> polyhedron
        - ``trans_inv`` -- a function: polyhedron -> polyhedron
        - ``ieq`` -- list, an inequality. An entry equal to "[-1,7,3,4]"
          represents the inequality 7x_1+3x_2+4x_3>= 1.

        OUTPUT:

            dict of polyhedron partitions with keys giving the return time

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, rotation_mod
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: u = rotation_mod(0, 1/3, 1, QQ)
            sage: u_inv = rotation_mod(0, 2/3, 1, QQ)
            sage: ieq = [h, -1, 0]   # x0 <= h
            sage: P.induced_in_partition(u, u_inv, ieq)
            {3: Polyhedron partition of 4 atoms with 4 letters}

        ::

            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq2 = [1/2, -1, 0]   # x0 <= 1/2
            sage: d = P.induced_in_partition(u, u_inv, ieq2)
            sage: d
            {1: Polyhedron partition of 2 atoms with 2 letters,
             2: Polyhedron partition of 3 atoms with 3 letters,
             3: Polyhedron partition of 4 atoms with 4 letters}
        """
        out_partition = self.induced_out_partition(trans, ieq)
        in_partition = {}
        for i,P in out_partition.items():
            for _ in range(i):
                P = P.apply_transformation(trans_inv)
            in_partition[i] = P
        return in_partition

    def induced_partition(self, trans, trans_inv, ieq):
        r"""
        Returns the partition of the induced transformation on the domain.

        INPUT:

        - ``trans`` -- a function: polyhedron -> polyhedron
        - ``trans_inv`` -- a function: polyhedron -> polyhedron
        - ``ieq`` -- list, an inequality. An entry equal to "[-1,7,3,4]"
          represents the inequality 7x_1+3x_2+4x_3>= 1.

        OUTPUT:

            - a polyhedron partition
            - dict, a substitution

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, rotation_mod
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: u = rotation_mod(0, 1/3, 1, QQ)
            sage: u_inv = rotation_mod(0, 2/3, 1, QQ)
            sage: ieq = [h, -1, 0]   # x0 <= h
            sage: Q,sub = P.induced_partition(u, u_inv, ieq)
            sage: Q
            Polyhedron partition of 4 atoms with 4 letters
            sage: sub
            {0: [0, 2, 2], 1: [1, 2, 2], 2: [1, 2, 3], 3: [1, 3, 3]}

        ::

            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq2 = [1/2, -1, 0]   # x0 <= 1/2
            sage: Q,sub = P.induced_partition(u, u_inv, ieq2)
            sage: Q
            Polyhedron partition of 9 atoms with 9 letters
            sage: sub
            {0: [0],
             1: [0, 2, 2],
             2: [1],
             3: [1, 2, 2],
             4: [1, 2, 3],
             5: [1, 3, 3],
             6: [2, 2],
             7: [2, 3],
             8: [3, 3]}

        Irrationnal rotations::

            sage: z = polygen(QQ, 'z') #z = QQ['z'].0 # same as
            sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: phi = K.gen()
            sage: h = 1/phi^2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s}, base_ring=K)
            sage: u = rotation_mod(0, 1/phi, 1, K)
            sage: u_inv = rotation_mod(0, 1/phi^2, 1, K)
            sage: ieq = [h, -1, 0]   # x0 <= h
            sage: P1,sub01 = P.induced_partition(u, u_inv, ieq)
            sage: P1
            Polyhedron partition of 7 atoms with 7 letters
            sage: sub01
            {0: [0, 2],
             1: [0, 2, 2],
             2: [1, 2],
             3: [1, 2, 2],
             4: [1, 3],
             5: [1, 3, 2],
             6: [1, 3, 3]}

        We do the induction on a smaller domain::

            sage: ieq2 = [1/phi^3, -1, 0]   # x0 <= h
            sage: P2,sub02 = P.induced_partition(u, u_inv, ieq2)
            sage: P2
            Polyhedron partition of 10 atoms with 10 letters
            sage: sub02
            {0: [0, 2, 0, 2, 2],
             1: [0, 2, 1, 2, 2],
             2: [0, 2, 2],
             3: [1, 2, 1, 2, 2],
             4: [1, 2, 1, 3, 2],
             5: [1, 2, 2],
             6: [1, 3, 1, 3, 2],
             7: [1, 3, 1, 3, 3],
             8: [1, 3, 2],
             9: [1, 3, 3]}

        We check that inductions commute::

            sage: u1 = rotation_mod(0, phi^-3, phi^-2, K)
            sage: u1_inv = rotation_mod(0, phi^-4, phi^-2, K)
            sage: P2_alt,sub12 = P1.induced_partition(u1, u1_inv, ieq2)
            sage: P2_alt
            Polyhedron partition of 10 atoms with 10 letters
            sage: P2_alt == P2
            True

        Up to a permutation of the alphabet, ``sub02`` and ``sub01*sub12``
        are equal::

            sage: s01 = WordMorphism(sub01)
            sage: s12 = WordMorphism(sub12)
            sage: s02 = WordMorphism(sub02)
            sage: s02
            WordMorphism: 0->02022, 1->02122, 2->022, 3->12122, 4->12132, 5->122, 6->13132, 7->13133, 8->132, 9->133
            sage: s01*s12 == s02
            True

        By chance, the above is true, but in general, we have::

            sage: perm = WordMorphism(P2.keys_permutation(P2_alt))
            sage: perm
            WordMorphism: 0->0, 1->1, 2->2, 3->3, 4->4, 5->5, 6->6, 7->7, 8->8, 9->9
            sage: s01*s12*perm == s02
            True
        """
        in_partition = self.induced_in_partition(trans, trans_inv, ieq)

        # Goal: we want two atoms to have the same key if they have
        # the same behavior under the induction

        # Solution: we construct a dict image of letter -> list of atoms
        from collections import defaultdict
        d = defaultdict(list)
        for return_time,P in in_partition.items():
            for garbage_key,p in P:
                p_copy = copy(p)
                w = []
                for _ in range(return_time):
                    w.append(self.code(p))
                    p = trans(p)
                d[tuple(w)].append(p_copy)

        # We construct the list of (key, atom) and the substitution
        L = []
        substitution = {}
        for key,(w,atoms) in enumerate(sorted(d.items())):
            for atom in atoms:
                L.append((key,atom))
                substitution[key] = list(w)

        return PolyhedronPartition(L), substitution


