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

def center_insphere_polytope(polytope, solver=None):
    r"""
    Return the center and radius of maximal inscribed sphere

    INPUT:

    - ``polytope`` -- polytope

    OUTPUT:

    a 2-tuple (center, radius)

    EXAMPLES::

        sage: from slabbe.polyhedron_partition import center_insphere_polytope
        sage: P = polytopes.associahedron(['A',3])
        sage: center_insphere_polytope(P)
        ([0.03553390593273766, 0.5355339059327378, 0.03553390593273766],
         1.4644660940672622)
        sage: center_insphere_polytope(P + vector((10,10,10)))
        ([10.035533905932738, 10.535533905932738, 10.035533905932738],
         1.4644660940672622)

    ::

        sage: P = Polyhedron([(0,0), (1,0), (0,10)])
        sage: center_insphere_polytope(P)
        ([0.47506218943955486, 0.47506218943955486], 0.47506218943955486)

    """
    from math import sqrt
    from sage.numerical.mip import MixedIntegerLinearProgram
    p = MixedIntegerLinearProgram(solver=solver)
    x = p.new_variable()
    d = p.new_variable()

    for ineq in polytope.inequalities_list():
        constant = ineq[0]
        coeffs = ineq[1:]
        norm = sqrt(sum(a**2 for a in coeffs))
        S = p.sum(a/norm*x[i] for (i,a) in enumerate(coeffs))
        p.add_constraint(constant/norm + S >= d[0])
    p.set_objective(d[0])

    radius = p.solve()
    soln = p.get_values(x)
    center = [soln[k] for k in sorted(x.keys())]
    return center, radius

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
        Returns the partition scaled by some factor.

        INPUT:

        - ``factor`` -- real number

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
            sage: P = P.translate((1,0))
            sage: P
            Polyhedron partition of 4 atoms with 4 letters

        """
        L = [(key, M*p) for key,p in self]
        return PolyhedronPartition(L)

    def translate(self, displacement):
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
            sage: P.translate((1,1))
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

    def plot(self, label_pos='insphere_center'):
        r"""
        INPUT:

        - ``label_pos`` -- string (default:``'insphere_center'``) or
          ``'center'`` for the center of the polytope

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
            if label_pos == 'insphere_center':
                center,radius = center_insphere_polytope(P)
            elif label_pos == 'center':
                center = P.center()
            else:
                raise ValueError('unknown value for label_pos(={})'.format(label_pos))
            G += text(key, center)
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

    def tikz(self, fontsize=r'\normalsize', scale=1, label_format = r'{}',
            label_pos='insphere_center', extra_code=''):
        r"""
        INPUT:

        - ``fontsize`` -- string (default: ``r'\normalsize'``
        - ``scale`` -- number (default: ``1``)
        - ``label_format`` -- string (default: ``r'{}'``) to be called with
          ``label_format.format(key)``
        - ``label_pos`` -- string (default:``'insphere_center'``) or
          ``'center'`` for the center of the polytope
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
        from sage.misc.functional import numerical_approx
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
            if label_pos == 'insphere_center':
                center,radius = center_insphere_polytope(P)
                center = tuple(numerical_approx(a,digits=5) for a in center)
            elif label_pos == 'center':
                center = P.center()
            else:
                raise ValueError('unknown value for label_pos(={})'.format(label_pos))
            lines.append(node_format.format(fontsize, center, label))
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
        the same image under the dictionnary and such that their union is
        convex.

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

    def refinement(self, other, certificate=False):
        r"""
        Return the polyhedron partition obtained by the intersection of the
        atoms of self with the atoms of other.

        Only atoms of positive volume are kept.

        INPUT:

        - ``other`` -- a polyhedron partition
        - ``certificate`` -- boolean (default:``False``), return a
          dictionnary for i:(p,q) if atom number i was obtained as the
          intersection of atoms p in self and q in other

        OUTPUT:

            a polyhedron partition
            
        or

            (polyhedron partition, dict)

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
        d = {}
        for ((a,p),(b,q)) in itertools.product(self, other):
            p_q = p.intersection(q)
            if p_q.is_full_dimensional():
                if certificate:
                    d[len(L)] = (a,b)
                L.append(p_q)

        P = PolyhedronPartition(L)
        if certificate:
            return P, d
        else:
            return P

