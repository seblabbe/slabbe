# -*- coding: utf-8 -*-
r"""
Polyhedron exchange transformations and induced transformations

EXAMPLES:

A polyhedron partition::

    sage: from slabbe import PolyhedronPartition
    sage: h = 1/3
    sage: p = Polyhedron([(0,h),(0,1),(h,1)])
    sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
    sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
    sage: s = Polyhedron([(h,0), (1,0), (1,h)])
    sage: P = PolyhedronPartition({0:p, 1:q, 2:r, 3:s})

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
from slabbe import PolyhedronPartition

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

        - Do we want to merge atoms mapped by the same translation?

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

    def __call__(self, p, niterations=1):
        r"""
        Apply the transformation.

        INPUT:

        - ``p`` -- vector or polyhedron or partition
        - ``niterations`` -- nonnegative integer (default: ``1``)

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

        Doing many iterations::

            sage: F((1/10, 1/10), niterations=5)
            (13/30, 1/10)
            sage: F(p, niterations=5)
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
            sage: F(P, niterations=5)
            Polyhedron partition of 3 atoms with 2 letters

        TESTS::

            sage: F((1/10, 1/10), niterations=0)
            (1/10, 1/10)

        ::

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
            Polyhedron partition of 6 atoms with 4 letters
            sage: uuP.volume()
            1

        """
        from sage.structure.element import Vector
        from sage.geometry.polyhedron.base import Polyhedron_base

        if not niterations >= 0:
            raise NotImplementedError("niterations(={}) supported only"
                    " when >= 0".format(niterations))

        if isinstance(p, (tuple, Vector)):
            p = self.ambient_space()(p)
            for _ in range(niterations):
                a = self._partition.code(Polyhedron([p]))
                t = self._translations[a]
                p = p + t
            return p

        elif isinstance(p, Polyhedron_base):
            for j in range(niterations):
                S = set(i for i,atom in self._partition if p <= atom)
                if len(S) == 1:
                    a = next(iter(S))
                    t = self._translations[a]
                    p = p + t
                elif len(S) > 1:
                    raise ValueError('During {}-th iteration, image of {} is not' 
                    ' well-defined as it belongs to many distinct atoms(={})' 
                    ' of the partition'.format(j,p,S))
                else:
                    raise ValueError('During {}-th iteration, image of polyhedron' 
                    ' (={}) is not defined as it overlaps distinct'
                    ' atoms of the partition'.format(j,p))
            return p

        elif isinstance(p, PolyhedronPartition):
            for _ in range(niterations):
                p,d = p.refinement(self._partition, certificate=True)
                L = []
                for key,atom in p:
                    good_key,pet_key = d[key]
                    a = self._partition.code(atom)
                    assert a == pet_key, "I think this is true (to check)"
                    t = self._translations[a]
                    L.append((good_key, atom + t))
                p = PolyhedronPartition(L)
            return p

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

        ::

            sage: h = 1/3
            sage: p = Polyhedron([(0,h),(0,1),(h,1)])
            sage: q = Polyhedron([(0,0), (0,h), (h,1), (h,0)])
            sage: r = Polyhedron([(h,1), (1,1), (1,h), (h,0)])
            sage: s = Polyhedron([(h,0), (1,0), (1,h)])
            sage: P = PolyhedronPartition([(0,p), (0,q), (1,r), (1,s)])
            sage: T = {0:(1-h,0), 1:(-h,0)}
            sage: F = PolyhedronExchangeTransformation(P, T)
            sage: F.image_partition()
            Polyhedron partition of 4 atoms with 2 letters

        """
        return PolyhedronPartition([(a,p+self._translations[a])
                                    for (a,p) in self._partition])

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

    def __mul__(self, other):
        r"""
        Return the product of polyhedron exchange transformations.

        INPUT:

        - ``other`` -- polyhedron exchange transformation

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
            with translations {0: (2/5, 0), 1: (-3/5, 0), 2: (-3/5, 0)}
            sage: F * F * F
            Polyhedron Exchange Transformation of
            Polyhedron partition of 4 atoms with 4 letters
            with translations {0: (3/5, 0), 1: (-2/5, 0), 2: (-2/5, 0), 3: (-2/5, 0)}

        """
        if not isinstance(other, PolyhedronExchangeTransformation):
            return NotImplemented
        R,d = self.image_partition().refinement(other.partition(),
                                                certificate=True)

        atoms_dict = {}
        trans_dict = {}
        for key,atom in R:
            (a,b) = d[key]
            atoms_dict[key] = atom - self._translations[a]
            trans_dict[key] = self._translations[a] + other._translations[b]

        P = PolyhedronPartition(atoms_dict)
        return PolyhedronExchangeTransformation(P, trans_dict)

    def __rmul__(self, factor):
        r"""
        Returns the PET scaled by some factor.

        INPUT:

        - ``factor`` -- real number

        OUTPUT:

        - polyhedron exchange transformation

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 4/5
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: T = {0:(1-h,0), 1:(-h,0)}
            sage: F = PolyhedronExchangeTransformation(P, T)
            sage: 2 * F
            Polyhedron Exchange Transformation of
            Polyhedron partition of 2 atoms with 2 letters
            with translations {0: (2/5, 0), 1: (-8/5, 0)}
            sage: -2 * F
            Polyhedron Exchange Transformation of
            Polyhedron partition of 2 atoms with 2 letters
            with translations {0: (-2/5, 0), 1: (8/5, 0)}

        """
        if isinstance(factor, PolyhedronExchangeTransformation):
            return NotImplemented
        P = factor * self.partition()
        trans_dict = {key:factor*t for (key,t) in self.translations().items()}
        return PolyhedronExchangeTransformation(P, trans_dict)

    def translate_domain(self, displacement):
        """
        Return the PET on a domain translated by some displacement.

        INPUT:

        - ``displacement`` -- a displacement vector or a list/tuple of
          coordinates that determines a displacement vector.

        OUTPUT:

        The translated PET

        EXAMPLES::

            sage: from slabbe import PolyhedronPartition, PolyhedronExchangeTransformation
            sage: h = 4/5
            sage: p = Polyhedron([(0,0),(h,0),(h,1),(0,1)])
            sage: q = Polyhedron([(1,0),(h,0),(h,1),(1,1)])
            sage: P = PolyhedronPartition({0:p, 1:q})
            sage: T = {0:(1-h,0), 1:(-h,0)}
            sage: F = PolyhedronExchangeTransformation(P, T)
            sage: Ft = F.translate_domain((3,1))
            sage: Ft
            Polyhedron Exchange Transformation of
            Polyhedron partition of 2 atoms with 2 letters
            with translations {0: (1/5, 0), 1: (-4/5, 0)}
            sage: Ft.domain().vertices()
            (A vertex at (3, 1),
             A vertex at (3, 2),
             A vertex at (4, 1),
             A vertex at (4, 2))

        """
        P = self.partition().translate(displacement)
        trans_dict = {key:t for (key,t) in self.translations().items()}
        return PolyhedronExchangeTransformation(P, trans_dict)

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

        # good side of the hyperplane
        half_polyhedron = Polyhedron(ieqs=[ieq])
        half = PolyhedronPartition([half_polyhedron])

        # the other side of the hyperplane
        other_half_polyhedron = Polyhedron(ieqs=[[-a for a in ieq]])
        other_half = PolyhedronPartition([other_half_polyhedron])

        # the window is the intersection of the domain with the half space
        W = self.domain().intersection(half_polyhedron)
        if W.volume() == 0:
            raise ValueError("Inequality {} does not intersect partition "
                    "(={})".format(half_polyhedron.inequalities()[0], partition))

        # Compute the induced partition and associated return words
        Q = PolyhedronPartition([(tuple(), W)])
        S = []
        self_inv = self.inverse()
        while len(Q):
            Q = self_inv(Q)
            # Compute the refinement of P and Q (concatenate the labels)
            PQ,d = partition.refinement(Q, certificate=True)
            Q = PolyhedronPartition([((d[i][0],)+d[i][1], q) for (i,q) in PQ])
            # Take what has returned to the window (keep labels from Q only)
            Q_returned,d = Q.refinement(half, certificate=True)
            S.extend((d[i][0], q) for (i,q) in Q_returned)
            # Continue with what is left (keep labels from Q only)
            Q,d = Q.refinement(other_half, certificate=True)
            Q = PolyhedronPartition([(d[i][0], q) for (i,q) in Q])

        # We sort the keys and relabel them with nonnegative integers
        from slabbe.finite_word import sort_word_by_length_lex_key
        return_words = set(w for (w,q) in S)
        sorted_return_words = sorted(return_words, key=sort_word_by_length_lex_key)
        key_to_word = {key:list(w) for (key,w) in enumerate(sorted_return_words)}
        word_to_key = {w:key for (key,w) in enumerate(sorted_return_words)}
        induced_partition = PolyhedronPartition([(word_to_key[w],q) for (w,q) in S])

        # Build a substitution2d if desired
        if substitution_type == 'dict':
            sub = key_to_word
        elif substitution_type == 'column':
            from slabbe import Substitution2d
            sub = Substitution2d.from_1d_column_substitution(key_to_word)
        elif substitution_type == 'row':
            from slabbe import Substitution2d
            sub = Substitution2d.from_1d_row_substitution(key_to_word)
        else:
            raise ValueError('Unknown value for substitution_type'
                    ' (={})'.format(substitution_type))

        return induced_partition, sub

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
            {0: [0], 1: [0, 1], 2: [0, 0, 1]}

        """
        newP, sub = self.induced_partition(ieq)
        T = self.translations()
        newT = {a:sum(T[b] for b in sub[a]) for a in sub}
        return PolyhedronExchangeTransformation(newP, newT), sub

    def cylinder(self, word, partition=None):
        r"""
        Return the region associated to the coding word.

        INPUT:

        - ``word`` -- list
        - ``partition`` -- polyhedron partition (default:``None``), if
          None, it uses the domain partition of the transformation

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
            {0}

        ::

            sage: u.cylinder([1,1], P)
            Polyhedron partition of 2 atoms with 2 letters
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
        if not word:
            return partition
        trans_inv = self.inverse()
        reversed_word = reversed(word)
        P = partition[next(reversed_word)]
        for a in reversed_word:
            P = trans_inv(P)
            P = partition[a].refinement(P)
        return P

    def cylinders(self, size, partition=None):
        r"""
        Return the cylinders of given size.

        INPUT:

        - ``size`` -- nonnegative integer
        - ``partition`` -- polyhedron partition (default:``None``), if
          None, it uses the domain partition of the transformation

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
            [{()}, {0, 1}, {0, 1, 2}, {0, 1, 2}, {0, 1, 2}]

        """
        # Default partition
        if partition is None:
            partition = self.partition()

        if size == 0:
            return PolyhedronPartition([(tuple(), self.domain())])

        P = partition
        trans_inv = self.inverse()
        for i in range(size-1):
            P = trans_inv(P)
            P = partition.refinement(P)
        return P

