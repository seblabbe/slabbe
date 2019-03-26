# -*- coding: utf-8 -*-
r"""
Polyhedron partition, polyhedron exchange transformations and induced transformations

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

    sage: from slabbe import PolyhedronExchangeTransformation as PET
    sage: base = identity_matrix(2)
    sage: translation = vector((2/3, 0))
    sage: u = PET.toral_translation(base, translation)
    sage: Q = u(P)
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
    sage: base = identity_matrix(2)
    sage: translation = vector((1/phi, 0))
    sage: u = PET.toral_translation(base, translation)
    sage: ieq = [h, -1, 0]   # x0 <= h
    sage: P1,sub01 = u.induced_partition(ieq, P)
    sage: P1
    Polyhedron partition of 7 atoms with 7 letters
    sage: sub01
    {0: [0, 2],
     1: [1, 2],
     2: [1, 3],
     3: [0, 2, 2],
     4: [1, 2, 2],
     5: [1, 3, 2],
     6: [1, 3, 3]}

AUTHORS:

- Sébastien Labbé, November 2017, initial version of polyhedron partitions
- Sébastien Labbé, January 2019, added a class for polyhedron exchange transformations
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
from copy import copy
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
            if p_q.volume() > 0:
                if key_fn is None:
                    L.append(p_q)
                else:
                    new_key = key_fn(a,b)
                    L.append((new_key, p_q))
        return PolyhedronPartition(L)
        
class PolyhedronExchangeTransformation(object):
    r"""
    Polyhedron Exchange Transformation (PET).

    INPUT:

    - ``partition`` -- a polyhedron partition
    - ``translations`` -- list or dict

    EXAMPLES::

        sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
        sage: h = 1/3
        sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
        sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
        sage: P = PolyhedronPartition({0:p, 1:q})
        sage: T = {0:(1-h,0), 1:(-h,0)}
        sage: PolyhedronExchangeTransformation(P, T)
        Polyhedron Exchange Transformation of 
        Polyhedron partition of 2 atoms with 2 letters
        with translations {0: (2/3, 0), 1: (-1/3, 0)}

    .. TODO::

        - Code the __pow__ methods.

        X Code the __mul__ methods.

        - Do we want to merge atoms mapped by the same translation?

        X Induction should return the induced transformation somehow.

        - Check mathematically that induced_in_partition and
          induced_out_partition are doing ok

        - Add a ploting function with seperated domain/codomain

    REFERENCES:

    - Schwartz, Richard Evan. The Octagonal PETs. First Edition edition.
      Providence, Rhode Island: American Mathematical Society, 2014.
    """
    def __init__(self, partition, translations):
        r"""

        """
        if isinstance(partition, PolyhedronPartition):
            self._partition = partition
        else:
            raise TypeError('partition(={}) must be a '
                            'PolyhedronPartition'.format(partition))
        self._ambient_space = self._partition.ambient_space()

        if isinstance(translations, list):
            self._translations = {a:self._ambient_space(t) for (a,t) in enumerate(translations)}
        elif isinstance(translations, dict):
            self._translations = {a:self._ambient_space(t) for (a,t) in translations.items()}
        else:
            raise TypeError('translations(={}) must be a '
                            'list or dict'.format(translations))

    @classmethod
    def toral_translation(cls, base, translation, fundamental_domain=None):
        r"""
        Return a polyhedron exchange transformation defined by a translation on
        a d-dimensional torus.

        INPUT:

        - ``base`` -- matrix, the columns are the base of a lattice
        - ``translation`` -- vector, translation vector
        - ``fundamental_domain`` -- polyhedron or ``None`` (default:
          ``None``), if ``None`` the parallelotope defined by ``base`` is
          used.

        OUTPUT:

            a polyhedron exchange transformation on the fundamental domain of
            the lattice

        EXAMPLES::

            sage: from slabbe import PolyhedronExchangeTransformation as PET
            sage: base = diagonal_matrix((1,1))
            sage: translation = vector((1/5, 1/3))
            sage: T = PET.toral_translation(base, translation)
            sage: T
            Polyhedron Exchange Transformation of
            Polyhedron partition of 4 atoms with 4 letters
            with translations {0: (1/5, 1/3), 1: (1/5, -2/3), 2: (-4/5, 1/3), 3: (-4/5, -2/3)}
            sage: T.partition()
            Polyhedron partition of 4 atoms with 4 letters

        Some preliminary definitions::

            sage: z = polygen(QQ, 'z') #z = QQ['z'].0 # same as
            sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: phi = K.gen()
            sage: vertices = ((-phi + 2, phi - 1), (-phi + 2, 1), (phi - 1, 1))
            sage: p = Polyhedron(vertices, base_ring=K)

        A translation +1 modulo phi on the x coordinate::

            sage: base = diagonal_matrix((phi,phi))
            sage: translation = vector((1, 0))
            sage: t0 = PET.toral_translation(base, translation)
            sage: t0
            Polyhedron Exchange Transformation of
            Polyhedron partition of 2 atoms with 2 letters
            with translations {0: (1, 0), 1: (-phi + 1, 0)}
            sage: t0(p).vertices()
            (A vertex at (-phi + 3, phi - 1),
             A vertex at (-phi + 3, 1),
             A vertex at (phi, 1))

        The inverse map::

            sage: t0.inverse()
            Polyhedron Exchange Transformation of
            Polyhedron partition of 2 atoms with 2 letters
            with translations {0: (-1, 0), 1: (phi - 1, 0)}
            sage: t0(p) == p
            False
            sage: t0.inverse()(t0(p)) == p
            True

        A rotation modulo 1 on the y coordinate::

            sage: base = diagonal_matrix((phi,phi))
            sage: translation = vector((0, 1))
            sage: t1 = PET.toral_translation(base, translation)
            sage: t1(p).vertices()
            (A vertex at (-phi + 2, 0),
             A vertex at (-phi + 2, -phi + 2),
             A vertex at (phi - 1, -phi + 2))

        It works if the translation is larger than the fundamental domain::

            sage: base = diagonal_matrix((1,1))
            sage: translation = vector((phi, 0))
            sage: t2 = PET.toral_translation(base, translation)
            sage: t2(p).vertices()
            (A vertex at (0, phi - 1), 
             A vertex at (0, 1), 
             A vertex at (2*phi - 3, 1))

        The domain is the fundamental domain of the given lattice::

            sage: base = diagonal_matrix((phi^-2,1))
            sage: translation = vector((phi^-3, 0))
            sage: t3 = PET.toral_translation(base, translation)
            sage: t3.domain().vertices()
            (A vertex at (-phi + 2, 0),
             A vertex at (-phi + 2, 1),
             A vertex at (0, 0),
             A vertex at (0, 1))

        The fundamental domain can be given as input. For example, it can
        be a translated copy of the base parallelotope::

            sage: base = diagonal_matrix((1,1))
            sage: translation = vector((1/5, 1/3))
            sage: F = polytopes.parallelotope(base)
            sage: T = PET.toral_translation(base, translation, F-vector((1/10,1/10)))

        But it does not always work well yet, for example for other shape
        of fundamental domains::

            sage: m = matrix(2, (1,1,0,1))
            sage: T = PET.toral_translation(base, translation, m*F)
            Traceback (most recent call last):
            ...
            NotImplementedError: Volume of the partition is 41/45 but the
            fundamental domain as volume 1. The code does not handle this
            case properly yet.

        .. TODO::

            Fix the above when the fundamental domain is far from the
            lattice base.
        """
        from sage.geometry.polyhedron.library import polytopes
        from sage.modules.free_module_element import vector
        from sage.functions.other import floor

        # Compute the representent of the translation inside the base
        v = base.inverse() * translation
        v_floor = vector(map(floor, v))
        translation -= base * v_floor

        # The fundamental domain
        base_parallelotope = polytopes.parallelotope(base.columns())
        if fundamental_domain is None:
            fundamental_domain = base_parallelotope
        FD = fundamental_domain # shorcut

        # Computing the partitions and translations
        atoms = {}
        trans = {}
        for vertex in base_parallelotope.vertices():
            t = vertex.vector() - translation
            I = FD.intersection(FD.translation(t))
            if I.volume():
                k = len(atoms)
                atoms[k] = I
                trans[k] = -t
        partition = PolyhedronPartition(atoms)

        if partition.volume() != FD.volume():
            raise NotImplementedError('Volume of the partition is {} but'
                    ' the fundamental domain as volume {}. The code does'
                    ' not handle this case properly yet.'.format(partition.volume(), 
                                                             FD.volume()))

        return PolyhedronExchangeTransformation(partition, trans)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 1/3
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: T = {0:(1-h,0), 1:(-h,0)}
            sage: F = PolyhedronExchangeTransformation(P, T)
            sage: F
            Polyhedron Exchange Transformation of 
            Polyhedron partition of 2 atoms with 2 letters
            with translations {0: (2/3, 0), 1: (-1/3, 0)}
        """
        return ('Polyhedron Exchange Transformation of\n{}\nwith '
               'translations {}').format(self._partition, self._translations)

    def partition(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 1/3
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: d = {0:(1-h,0), 1:(-h,0)}
            sage: T = PolyhedronExchangeTransformation(P, d)
            sage: T.partition()
            Polyhedron partition of 2 atoms with 2 letters
        """
        return self._partition

    def plot(self):
        r"""

        .. TODO::

            return two copy side-to-side of the domain and codomain instead
            of arrows.

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 1/3
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: d = {0:(1-h,0), 1:(-h,0)}
            sage: T = PolyhedronExchangeTransformation(P, d)
            sage: T.plot()
            Graphics object consisting of 16 graphics primitives
        """
        from random import random
        from sage.plot.arrow import arrow
        from sage.modules.free_module_element import vector
        d = self.ambient_space().dimension()

        # computing the range of the domain in each dimension
        V = self.domain().vertices()
        MAX = map(max, *V)
        MIN = map(min, *V)
        H = [b-a for (a,b) in zip(MIN,MAX)]

        G = self.partition().plot()
        for key,p in self.partition():
            t = self._translations[key]
            small_noise_vector = vector([H[i]*.1*random() for i in range(d)])
            center = p.center() + small_noise_vector
            G += arrow(center, center+t, color='green')
        return G

    def domain(self):
        r"""
        Return the domain of the exchange transformation.

        OUTPUT:

            a polyhedron

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 1/3
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: T = {0:(1-h,0), 1:(-h,0)}
            sage: F = PolyhedronExchangeTransformation(P, T)
            sage: F.domain()
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
            sage: F.domain().vertices()
            (A vertex at (0, 0),
             A vertex at (0, 1),
             A vertex at (1, 0),
             A vertex at (1, 1))
        """
        return self.partition().domain()

    def translations(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 1/3
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: T = {0:(1-h,0), 1:(-h,0)}
            sage: F = PolyhedronExchangeTransformation(P, T)
            sage: F.translations()
            {0: (2/3, 0), 1: (-1/3, 0)}
        """
        return self._translations

    def ambient_space(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 1/3
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: T = {0:(1-h,0), 1:(-h,0)}
            sage: F = PolyhedronExchangeTransformation(P, T)
            sage: F.ambient_space()
            Vector space of dimension 2 over Rational Field
        """
        return self._ambient_space

    def merge_atoms_with_same_translation(self):
        r"""
        Return a new partition into convex polyhedrons where atoms mapped
        by the same translation are merged if their union is convex.

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 1/3
            sage: p = Polyhedron([(0,0),(h,0),(h,h),(0,h)])
            sage: q = Polyhedron([(0,h),(h,h),(h,1),(0,1)])
            sage: r = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r})
            sage: d = {0:(1-h,0), 1:(1-h,0), 2:(-h,0)}
            sage: T = PolyhedronExchangeTransformation(P, d)
            sage: T
            Polyhedron Exchange Transformation of
            Polyhedron partition of 3 atoms with 3 letters
            with translations {0: (2/3, 0), 1: (2/3, 0), 2: (-1/3, 0)}
            sage: T.merge_atoms_with_same_translation()
            Polyhedron Exchange Transformation of
            Polyhedron partition of 2 atoms with 2 letters
            with translations {0: (2/3, 0), 2: (-1/3, 0)}

        """
        from collections import defaultdict

        # a dictionary translation vector -> list of keys
        d = defaultdict(list)
        for (key,t) in self.translations().items():
            t.set_immutable()
            d[t].append(key)

        # a dictionary translation vector -> a unique atom key (the minimum)
        d = {t:min(d[t]) for t in d}

        # a dictionary key -> new key
        d = {a:d[t] for a,t in self.translations().items()}

        # merged partition
        P = self.partition().merge_atoms(d)

        # a dictionary of translations with the new keys
        translation_dict = self.translations()
        translation_dict = {a:translation_dict[a] for a in d.values()}

        return PolyhedronExchangeTransformation(P, translation_dict)

    def __call__(self, p, key_fn=None):
        r"""
        Apply the transformation.

        INPUT:

        - ``p`` -- vector or polyhedron or partition
        - ``key_fn`` -- function to apply on pairs of labels, or ``None``

        OUTPUT:

            vector or polyhedron or partition of polyhedron

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 1/3
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: T = {0:(1-h,0), 1:(-h,0)}
            sage: F = PolyhedronExchangeTransformation(P, T)

        Image of a vector::

            sage: F((1/10, 1/10))
            (23/30, 1/10)

        Image of a polyhedron::

            sage: F(p)
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices

        Image of a partition into polyhedron::

            sage: F(P)
            Polyhedron partition of 2 atoms with 2 letters

        TESTS::

            sage: from slabbe import PolyhedronPartition
            sage: from slabbe import PolyhedronExchangeTransformation as PET
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: base = identity_matrix(2)
            sage: translation = vector((2/3, 0))
            sage: u = PET.toral_translation(base, translation)
            sage: uP = u(P); uP
            Polyhedron partition of 4 atoms with 4 letters
            sage: uP.volume()
            1
            sage: uuP = u(uP); uuP
            Polyhedron partition of 6 atoms with 6 letters
            sage: uuP.volume()
            1

        """
        from sage.structure.element import Vector
        from sage.geometry.polyhedron.base import Polyhedron_base

        if isinstance(p, (tuple, Vector)):
            p = self.ambient_space()(p)
            a = self._partition.code(Polyhedron([p]))
            t = self._translations[a]
            return p + t

        elif isinstance(p, Polyhedron_base):
            S = set(i for i,atom in self._partition if p <= atom)
            if len(S) == 1:
                a = next(iter(S))
                t = self._translations[a]
                return p + t
            elif len(S) > 1:
                raise ValueError('image of {} is not well-defined as it'
                ' belongs to many distinct atoms(={}) of the partition'.format(p,S))
            else:
                raise ValueError('image of polyhedron (={}) is not defined as it'
                ' overlaps distinct atoms of the partition'.format(p))

        elif isinstance(p, PolyhedronPartition):
            if key_fn is None:
                # key_fn = lambda a,b:a # the previous default
                def key_fn(a,b):
                    if not isinstance(a, tuple):
                        a = (a,)
                    if not isinstance(b, tuple):
                        b = (b,)
                    return a + b

            p = p.refinement(self._partition, key_fn=key_fn)
            L = []
            for key,atom in p:
                a = self._partition.code(atom)
                t = self._translations[a]
                L.append((key, atom + t))
            return PolyhedronPartition(L)

        else:
            raise TypeError('call undefined on input p(={})'.format(p))

    def image_partition(self):
        r"""
        Return the partition of the image.

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 1/3
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: T = {0:(1-h,0), 1:(-h,0)}
            sage: F = PolyhedronExchangeTransformation(P, T)
            sage: F.image_partition()
            Polyhedron partition of 2 atoms with 2 letters

        """
        return PolyhedronPartition({a:p+self._translations[a] 
                                 for (a,p) in self._partition})

    def inverse(self):
        r"""
        Return the inverse of self.

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 1/3
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: T = {0:(1-h,0), 1:(-h,0)}
            sage: F = PolyhedronExchangeTransformation(P, T)
            sage: F
            Polyhedron Exchange Transformation of 
            Polyhedron partition of 2 atoms with 2 letters
            with translations {0: (2/3, 0), 1: (-1/3, 0)}

        ::

            sage: F.inverse()
            Polyhedron Exchange Transformation of 
            Polyhedron partition of 2 atoms with 2 letters
            with translations {0: (-2/3, 0), 1: (1/3, 0)}

        """
        P = self.image_partition()
        T = {a:-t for (a,t) in self._translations.items()}
        return PolyhedronExchangeTransformation(P, T)

    def __mul__(self, other, key_fn=None):
        r"""
        Return the product of polyhedron exchange transformations.

        INPUT:

        - ``other`` -- polyhedron exchange transformation
        - ``key_fn`` -- function to apply on pairs of labels, or ``None``

        OUTPUT:

        - polyhedron exchange transformation with keys being a tuple of
          previous keys

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 4/5
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: T = {0:(1-h,0), 1:(-h,0)}
            sage: F = PolyhedronExchangeTransformation(P, T)
            sage: F * F
            Polyhedron Exchange Transformation of
            Polyhedron partition of 3 atoms with 3 letters
            with translations {(0, 1): (-3/5, 0), (1, 0): (-3/5, 0), (0, 0): (2/5, 0)}
            sage: F * F * F
            Polyhedron Exchange Transformation of
            Polyhedron partition of 4 atoms with 4 letters
            with translations {(1, 0, 0): (-2/5, 0), (0, 1, 0): (-2/5, 0),
                               (0, 0, 0): (3/5, 0), (0, 0, 1): (-2/5, 0)}

        """
        key_fn_tuple = lambda a,b:(a,b)
        R = self.image_partition().refinement(other.partition(),
                                              key_fn=key_fn_tuple)
        if key_fn is None:
            def key_fn(a,b):
                if not isinstance(a, tuple):
                    a = (a,)
                if not isinstance(b, tuple):
                    b = (b,)
                return a + b

        atoms_dict = {}
        trans_dict = {}
        for (a,b),atom in R:
            key = key_fn(a,b)
            atoms_dict[key] = atom - self._translations[a]
            trans_dict[key] = self._translations[a] + other._translations[b]

        P = PolyhedronPartition(atoms_dict)
        return PolyhedronExchangeTransformation(P, trans_dict)

    def induced_out_partition(self, ieq, partition=None):
        r"""
        Returns the output partition obtained as the induction of the
        transformation on the domain given by an inequality.

        Note: the output partition corresponds to the arrival partition in
        the domain, not the initial one.

        INPUT:

        - ``ieq`` -- list, an inequality. An entry equal to "[-1,7,3,4]"
          represents the inequality 7x_1+3x_2+4x_3>= 1.
        - ``partition`` -- polyhedron partition (default:``None``), if
          None, it uses the domain partition of the transformation

        OUTPUT:

            dict of polyhedron partitions with keys giving the return time

        EXAMPLES::

            sage: from slabbe import PolyhedronExchangeTransformation as PET
            sage: base = identity_matrix(2)
            sage: translation = vector((1/3, 0))
            sage: u = PET.toral_translation(base, translation)
            sage: ieq = [1/2, -1, 0]   # x0 <= 1/2
            sage: d = u.induced_out_partition(ieq)
            sage: [(i, d[i], d[i].alphabet()) for i in d]
            [(1, Polyhedron partition of 1 atoms with 1 letters, {(0,)}),
             (2, Polyhedron partition of 1 atoms with 1 letters, {(0, 1)}),
             (3, Polyhedron partition of 1 atoms with 1 letters, {(0, 0, 1)})]

        ::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq = [h, -1, 0]   # x0 <= h
            sage: d = u.induced_out_partition(ieq, P)
            sage: [(i, d[i], d[i].alphabet()) for i in d]
            [(3,
              Polyhedron partition of 4 atoms with 4 letters,
              {(0, 2, 2), (1, 2, 2), (1, 2, 3), (1, 3, 3)})]

        ::

            sage: ieq2 = [1/2, -1, 0]   # x0 <= 1/2
            sage: d = u.induced_out_partition(ieq2, P)
            sage: [(i, d[i], d[i].alphabet()) for i in d]
            [(1, Polyhedron partition of 2 atoms with 2 letters, {(0,), (1,)}),
             (2, Polyhedron partition of 3 atoms with 3 letters, {(2, 2), (2, 3), (3, 3)}),
             (3,
              Polyhedron partition of 4 atoms with 4 letters,
              {(0, 2, 2), (1, 2, 2), (1, 2, 3), (1, 3, 3)})]
            sage: Q = PolyhedronPartition(d[1].atoms()+d[2].atoms()+d[3].atoms())
            sage: Q.is_pairwise_disjoint()
            True

        ::

            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq3 = [-1/2, 1, 0]   # x0 >= 1/2
            sage: u.induced_out_partition(ieq3, P)
            {1: Polyhedron partition of 2 atoms with 2 letters,
             2: Polyhedron partition of 3 atoms with 3 letters,
             3: Polyhedron partition of 4 atoms with 4 letters}

        It is an error if the induced region is empty::

            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq4 = [-1/2, -1, 0]   # x0 <= -1/2
            sage: u.induced_out_partition(ieq4, P)
            Traceback (most recent call last):
            ...
            ValueError: Inequality An inequality (-2, 0) x - 1 >= 0 does
            not intersect P (=Polyhedron partition of 4 atoms with 4
            letters)

        The whole domain::

            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq5 = [1/2, 1, 0]   # x0 >= -1/2
            sage: d = u.induced_out_partition(ieq5, P)
            sage: [(i, d[i], d[i].alphabet()) for i in d]
            [(1,
              Polyhedron partition of 6 atoms with 4 letters, 
              {(0,), (1,), (2,), (3,)})]

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
            sage: base = identity_matrix(2)
            sage: translation = vector((1/phi, 0))
            sage: u = PET.toral_translation(base, translation)
            sage: ieq = [phi^-4, -1, 0]   # x0 <= phi^-4
            sage: d = u.induced_out_partition(ieq, P)
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

        # Default partition
        if partition is None:
            partition = self.partition()

        # initial refinement
        key_fn_left = lambda a,b:a
        key_fn_both = lambda a,b:a+(b,)
        P = partition.refinement(half_part, key_fn=key_fn_left)
        if len(P) == 0:
            raise ValueError("Inequality {} does not intersect P "
                    "(={})".format(half.inequalities()[0], partition))
        level = 1
        ans = {}
        P = self(P, key_fn=lambda a,b:(a,))
        while len(P):
            P_returned = P.refinement(half_part, key_fn=key_fn_left)
            if P_returned:
                ans[level] = P_returned
            # for what is remaining we do:
            P = P.refinement(other_half_part, key_fn=key_fn_left)
            P = P.refinement(partition, key_fn=key_fn_both)
            P = self(P, key_fn=key_fn_left)
            level += 1
        return ans

    def induced_in_partition(self, ieq, partition=None):
        r"""
        Returns the partition of the induced transformation on the domain.
        given by an inequality.

        INPUT:

        - ``ieq`` -- list, an inequality. An entry equal to "[-1,7,3,4]"
          represents the inequality 7x_1+3x_2+4x_3>= 1.
        - ``partition`` -- polyhedron partition (default:``None``), if
          None, it uses the domain partition of the transformation

        OUTPUT:

            dict of polyhedron partitions with keys giving the return time

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})

        ::

            sage: from slabbe import PolyhedronExchangeTransformation as PET
            sage: base = identity_matrix(2)
            sage: translation = vector((1/3, 0))
            sage: u = PET.toral_translation(base, translation)
            sage: ieq = [h, -1, 0]   # x0 <= h
            sage: u.induced_in_partition(ieq, P)
            {3: Polyhedron partition of 4 atoms with 4 letters}

        ::

            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq2 = [1/2, -1, 0]   # x0 <= 1/2
            sage: d = u.induced_in_partition(ieq2, P)
            sage: d
            {1: Polyhedron partition of 2 atoms with 2 letters,
             2: Polyhedron partition of 3 atoms with 3 letters,
             3: Polyhedron partition of 4 atoms with 4 letters}
        """
        out_partition = self.induced_out_partition(ieq, partition)
        in_partition = {}
        self_inv = self.inverse()
        for i,P in out_partition.items():
            for _ in range(i):
                P = self_inv(P)
                #P = self_inv(P, key_fn=lambda a,b:a)
            in_partition[i] = P
        return in_partition

    def induced_partition(self, ieq, partition=None, substitution_type='dict'):
        r"""
        Returns the partition of the induced transformation on the domain.

        INPUT:

        - ``ieq`` -- list, an inequality. An entry equal to "[-1,7,3,4]"
          represents the inequality 7x_1+3x_2+4x_3>= 1.
        - ``partition`` -- polyhedron partition (default:``None``), if
          None, it uses the domain partition of the transformation
        - ``substitution_type`` -- string (default:``'dict'``), if
          ``'column'`` or ``'row'``, it returns a substitution2d, otherwise
          it returns a dict.

        OUTPUT:

            - a polyhedron partition
            - a substitution2d or a dict

        EXAMPLES::

            sage: from slabbe import PolyhedronExchangeTransformation as PET
            sage: base = identity_matrix(2)
            sage: translation = vector((1/3, 0))
            sage: u = PET.toral_translation(base, translation)

        We compute the induced partition of a polyhedron exchange
        transformation on a subdomain given by an inequality::

            sage: ieq = [1/3, -1, 0]   # x0 <= 1/3
            sage: u.induced_partition(ieq)
            (Polyhedron partition of 1 atoms with 1 letters,
             {0: [0, 0, 1]})
            sage: ieq = [1/2, -1, 0]   # x0 <= 1/2
            sage: u.induced_partition(ieq)
            (Polyhedron partition of 3 atoms with 3 letters,
             {0: [0], 1: [0, 1], 2: [0, 0, 1]})

        The second output can be turned into a column or a row
        Substitution2d if desired::

            sage: u.induced_partition(ieq, substitution_type='row')
            (Polyhedron partition of 3 atoms with 3 letters,
             Substitution 2d: {0: [[0]], 1: [[0], [1]], 2: [[0], [0], [1]]})

        Now we construct a another coding partition::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})

        We use this other partition to compute the induced partition::

            sage: ieq = [h, -1, 0]   # x0 <= h
            sage: Q,sub = u.induced_partition(ieq, P)
            sage: Q
            Polyhedron partition of 4 atoms with 4 letters
            sage: sub
            {0: [0, 2, 2], 1: [1, 2, 2], 2: [1, 2, 3], 3: [1, 3, 3]}

        ::

            sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})
            sage: ieq2 = [1/2, -1, 0]   # x0 <= 1/2
            sage: Q,sub = u.induced_partition(ieq2, P)
            sage: Q
            Polyhedron partition of 9 atoms with 9 letters
            sage: sub
            {0: [0],
             1: [1],
             2: [2, 2],
             3: [2, 3],
             4: [3, 3],
             5: [0, 2, 2],
             6: [1, 2, 2],
             7: [1, 2, 3],
             8: [1, 3, 3]}

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
            sage: base = identity_matrix(2)
            sage: translation = vector((1/phi, 0))
            sage: u = PET.toral_translation(base, translation)
            sage: ieq = [h, -1, 0]   # x0 <= h
            sage: P1,sub01 = u.induced_partition(ieq, P)
            sage: P1
            Polyhedron partition of 7 atoms with 7 letters
            sage: sub01
            {0: [0, 2],
             1: [1, 2],
             2: [1, 3],
             3: [0, 2, 2],
             4: [1, 2, 2],
             5: [1, 3, 2],
             6: [1, 3, 3]}

        We do the induction on a smaller domain::

            sage: ieq2 = [1/phi^3, -1, 0]   # x0 <= h
            sage: P2,sub02 = u.induced_partition(ieq2, P)
            sage: P2
            Polyhedron partition of 10 atoms with 10 letters
            sage: sub02
            {0: [0, 2, 2],
             1: [1, 2, 2],
             2: [1, 3, 2],
             3: [1, 3, 3],
             4: [0, 2, 0, 2, 2],
             5: [0, 2, 1, 2, 2],
             6: [1, 2, 1, 2, 2],
             7: [1, 2, 1, 3, 2],
             8: [1, 3, 1, 3, 2],
             9: [1, 3, 1, 3, 3]}

        We check that inductions commute::

            sage: base = diagonal_matrix((phi^-2,1))
            sage: translation = vector((phi^-3, 0))
            sage: u1 = PET.toral_translation(base, translation)
            sage: P2_alt,sub12 = u1.induced_partition(ieq2, P1)
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
            WordMorphism: 0->022, 1->122, 2->132, 3->133, 4->02022, 5->02122, 6->12122, 7->12132, 8->13132, 9->13133
            sage: s01*s12 == s02
            True

        By chance, the above is true, but in general, we have::

            sage: perm = WordMorphism(P2.keys_permutation(P2_alt))
            sage: perm
            WordMorphism: 0->0, 1->1, 2->2, 3->3, 4->4, 5->5, 6->6, 7->7, 8->8, 9->9
            sage: s01*s12*perm == s02
            True
        """
        # Default partition
        if partition is None:
            partition = self.partition()

        in_partition = self.induced_in_partition(ieq, partition)

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
                    w.append(partition.code(p))
                    p = self(p)
                d[tuple(w)].append(p_copy)

        # We construct the list of (key, atom) and the substitution
        L = []
        sub = {}
        from slabbe.finite_word import sort_word_by_length_lex_key
        sorted_keys = sorted(d.keys(), key=sort_word_by_length_lex_key)
        for key,w in enumerate(sorted_keys):
            atoms = d[w]
            for atom in atoms:
                L.append((key,atom))
                sub[key] = list(w)

        # Build a substitution2d if desired
        if substitution_type == 'dict':
            pass
        elif substitution_type in ['column', 'row']:
            from slabbe import Substitution2d
            if substitution_type == 'column':
                sub = Substitution2d.from_1d_column_substitution(sub)
            elif substitution_type == 'row':
                sub = Substitution2d.from_1d_row_substitution(sub)
        else:
            raise ValueError('Unknown value for substitution_type'
                    ' (={})'.format(substitution_type))

        return PolyhedronPartition(L), sub

    def induced_transformation(self, ieq):
        r"""
        Return the induced transformation on the domain.

        INPUT:

        - ``ieq`` -- list, an inequality. An entry equal to "[-1,7,3,4]"
          represents the inequality 7x_1+3x_2+4x_3>= 1.

        OUTPUT:

            - a polyhedron exchange transformation on the subdomain
            - a substitution (dict)

        EXAMPLES::

            sage: from slabbe import PolyhedronExchangeTransformation as PET
            sage: base = identity_matrix(2)
            sage: translation = vector((1/3, 0))
            sage: u = PET.toral_translation(base, translation)

        We compute the induced transformation of a polyhedron exchange
        transformation on a subdomain given by an inequality::

            sage: ieq = [1/2, -1, 0]   # x0 <= 1/2
            sage: T,sub = u.induced_transformation(ieq)
            sage: T
            Polyhedron Exchange Transformation of
            Polyhedron partition of 3 atoms with 3 letters
            with translations {0: (1/3, 0), 1: (-1/3, 0), 2: (0, 0)}
            sage: sub
            {0: (0,), 1: (0, 1), 2: (0, 0, 1)}

        """
        out_partition = self.induced_out_partition(ieq, partition=None)

        atoms_dict = {}
        trans_dict = {}
        sub = {}
        natural_numbers = itertools.count()
        for return_time,Q in out_partition.items():
            for w,atom in Q:
                i = next(natural_numbers)
                trans_dict[i] = t = sum(self._translations[a] for a in w)
                atoms_dict[i] = atom - t
                sub[i] = w

        P = PolyhedronPartition(atoms_dict)
        return PolyhedronExchangeTransformation(P, trans_dict), sub

    def cylinder(self, word, partition=None, key_fn=None):
        r"""
        Return the region associated to the coding word.

        INPUT:

        - ``word`` -- list
        - ``partition`` -- polyhedron partition (default:``None``), if
          None, it uses the domain partition of the transformation
        - ``key_fn`` -- function (default:``lambda a,b:(a,b)``), the
          concatenation function 

        OUTPUT:

            polyhedron partition

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition
            sage: h = 1/2
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (1,1), (1,h), (h,0)])
            sage: r = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([p,q,r])

        ::

            sage: from slabbe import PolyhedronExchangeTransformation as PET
            sage: base = identity_matrix(2)
            sage: translation = vector((1/3, 0))
            sage: u = PET.toral_translation(base, translation)
            sage: c = u.cylinder([2,2], P); c
            Polyhedron partition of 1 atoms with 1 letters
            sage: c.alphabet()
            {(2, 2)}

        ::

            sage: u.cylinder([1,1], P)
            Polyhedron partition of 2 atoms with 1 letters
            sage: u.cylinder([1], P)
            Polyhedron partition of 1 atoms with 1 letters

        Cylinders of words of length 0::

            sage: u.cylinder([], P).volume()
            1

        Cylinders of words of length 1::

            sage: C1 = [u.cylinder([a], P).volume() for a in range(3)]
            sage: C1
            [1/8, 3/4, 1/8]
            sage: sum(C1)
            1

        Cylinders of words of length 2::

            sage: import itertools
            sage: L2 = itertools.product(range(3),repeat=2)
            sage: C2 = [u.cylinder([a,b], P).volume() for (a,b) in L2]
            sage: C2
            [1/72, 1/9, 0, 1/9, 19/36, 1/9, 0, 1/9, 1/72]
            sage: sum(C2)
            1

        Cylinders of words of length 3::

            sage: L3 = itertools.product(range(3),repeat=3)
            sage: C3 = [u.cylinder([a,b,c], P).volume() for (a,b,c) in L3]
            sage: sum(C3)
            1

        TESTS::

            sage: u.cylinder([0,0,0], P)
            Polyhedron partition of 0 atoms with 0 letters
            sage: u.cylinder([2,3], P)
            Polyhedron partition of 0 atoms with 0 letters
            sage: u.cylinder([2,1], P)
            Polyhedron partition of 1 atoms with 1 letters
            sage: u.cylinder([], P)
            Polyhedron partition of 3 atoms with 3 letters

        """
        # Default partition
        if partition is None:
            partition = self.partition()
        # Default key_fn
        if key_fn is None:
            key_fn = lambda a,b:(a,b)

        if not word:
            return partition
        trans_inv = self.inverse()
        reversed_word = reversed(word)
        P = partition[next(reversed_word)]
        for a in reversed_word:
            P = trans_inv(P, key_fn=lambda a,b:a)
            P = partition[a].refinement(P, key_fn=key_fn)
        return P

    def cylinders(self, size, partition=None, key_fn=None):
        r"""
        Return the cylinders of given size.

        INPUT:

        - ``size`` -- nonnegative integer
        - ``partition`` -- polyhedron partition (default:``None``), if
          None, it uses the domain partition of the transformation
        - ``key_fn`` -- function (default:``lambda a,b:a+b`` and every key
          of atoms of the partition is changed into a singleton tuple),
          the concatenation function 

        OUTPUT:

            polyhedron partition

        EXAMPLES::

            sage: from slabbe import PolyhedronExchangeTransformation as PET
            sage: base = identity_matrix(2)
            sage: translation = vector((1/3, 0))
            sage: u = PET.toral_translation(base, translation)
            sage: [u.cylinders(i) for i in range(5)]
            [Polyhedron partition of 1 atoms with 1 letters,
             Polyhedron partition of 2 atoms with 2 letters,
             Polyhedron partition of 3 atoms with 3 letters,
             Polyhedron partition of 3 atoms with 3 letters,
             Polyhedron partition of 3 atoms with 3 letters]
            sage: [u.cylinders(i).alphabet() for i in range(5)]
            [{()},
             {(0,), (1,)},
             {(0, 0), (0, 1), (1, 0)},
             {(0, 0, 1), (0, 1, 0), (1, 0, 0)},
             {(0, 0, 1, 0), (0, 1, 0, 0), (1, 0, 0, 1)}]

        """
        # Default partition
        if partition is None:
            partition = self.partition()
        # Default key_fn
        if key_fn is None:
            key_fn = lambda a,b:a+b
            partition = PolyhedronPartition([((a,),p) for a,p in partition])

        if size == 0:
            return PolyhedronPartition([(tuple(), self.domain())])

        P = partition
        trans_inv = self.inverse()
        for i in range(size-1):
            P = trans_inv(P, key_fn=lambda a,b:a)
            P = partition.refinement(P, key_fn=key_fn)
        return P

