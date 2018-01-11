# -*- coding: utf-8 -*-
r"""
Digital geometry primitives

Subsets of ZZ^d with the edge relation +e_i and -e_i.

EXAMPLES::

    sage: from slabbe import DiscreteSubset
    sage: DiscreteSubset(dimension=2)
    Subset of ZZ^2
    sage: DiscreteSubset(dimension=4)
    Subset of ZZ^4

A subset from an iterable::

    sage: L = [(0,0,0,0), (1,0,0,0), (2,0,0,0), (3,0,0,0)]
    sage: s = DiscreteSubset.from_subset(L)
    sage: s
    Subset of ZZ^4

A discrete 2d disk::

    sage: D = DiscreteSubset(dimension=2, predicate=lambda (x,y) : x^2 + y^2 < 4)
    sage: D.list()
    [(0, 0), (0, 1), (0, -1), (1, 0), (-1, 0), (-1, 1), (1, -1), (1, 1), (-1, -1)]
    sage: D
    Subset of ZZ^2

A discrete 3d ball::

    sage: predicate = lambda (x,y,z) : x^2 + y^2 + z^2 <= 4
    sage: D = DiscreteSubset(dimension=3, predicate=predicate)
    sage: D
    Subset of ZZ^3
    sage: (0,0,0) in D
    True
    sage: (10,10,10) in D
    False
    sage: len(D.list())
    33
    sage: D.plot()    # optional long

A discrete 4d hyperplane::

    sage: predicate = lambda (x,y,z,w) : 0 <= 2*x + 3*y + 4*z + 5*w < 14
    sage: D = DiscreteSubset(dimension=4, predicate=predicate)
    sage: D
    Subset of ZZ^4
    sage: D.an_element()
    (0, 0, 0, 0)

A 2d discrete box::

    sage: from slabbe import DiscreteBox
    sage: b = DiscreteBox([-5,5], [-5,5])
    sage: b
    Box: [-5, 5] x [-5, 5]
    sage: b.plot()       # optional long

A 3d discrete box::

    sage: b = DiscreteBox([-2,2], [-5,5], [-5,5])
    sage: b
    Box: [-2, 2] x [-5, 5] x [-5, 5]
    sage: b.plot()       # optional long

The intersection of two discrete objects of the same dimension::

    sage: circ = DiscreteSubset(dimension=2, predicate=lambda p: p[0]^2+p[1]^2<=100)
    sage: b = DiscreteBox([0,10], [0,10])
    sage: I = circ & b
    sage: I
    Intersection of the following objects:
    Subset of ZZ^2
    [0, 10] x [0, 10]
    sage: I.an_element()
    (0, 0)
    sage: I.plot()      # optional long

A discrete tube (preimage of a discrete box by a matrix)::

    sage: from slabbe import M3to2
    sage: M3to2
    [-0.866025403784439  0.866025403784439  0.000000000000000]
    [-0.500000000000000 -0.500000000000000   1.00000000000000]
    sage: from slabbe import DiscreteTube
    sage: tube = DiscreteTube([-5,5],[-5,5], projmat=M3to2)
    sage: tube
    DiscreteTube: Preimage of [-5, 5] x [-5, 5] by a 2 by 3 matrix
    sage: it = iter(tube)
    sage: [next(it) for _ in range(4)]
    [(0, 0, 0), (1, 0, 0), (0, 0, 1), (0, 0, -1)]

TODO:

    - Code Complement
    - The method projection_matrix should be outside of the class?
    - DiscreteTube should have a method projection_matrix
    - Their should be an input saying whether the object is connected or
      not and what kind of neighbor connectedness
"""
#*****************************************************************************
#       Copyright (C) 2010-2016 Sébastien Labbé <slabqc@gmail.com>
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
from sage.functions.other import sqrt
from sage.matrix.constructor import matrix
from sage.rings.real_mpfr import RR
from sage.modules.free_module_element import vector
from sage.structure.sage_object import SageObject
from sage.modules.free_module import FreeModule
from sage.rings.integer_ring import ZZ
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
from sage.misc.decorators import options
from sage.misc.cachefunc import cached_method
from sage.misc.functional import round
from sage.misc.latex import LatexExpr
from sage.plot.graphics import Graphics
from sage.plot.point import point
from sage.plot.circle import circle
from sage.plot.plot3d.shapes2 import text3d
from sage.plot.line import line
from sage.plot.text import text
from sage.plot.plot3d.platonic import cube
from slabbe.tikz_picture import TikzPicture
from slabbe.matrices import M3to2

################################################
# Convex boundary of a set of 2d points
################################################
def convex_boundary(L):
    r"""
    EXAMPLES::

        sage: from slabbe.discrete_subset import convex_boundary
        sage: convex_boundary([(3,4), (1,2), (3,5)])
        [(3, 5), (1, 2), (3, 4)]

    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    P = Polyhedron(L)
    from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d
    return map(tuple, cyclic_sort_vertices_2d(P.vertices()))
    #if P.ambient_dim() != 2:
    #    raise ValueError("polyhedron ambient dim (=%s) must be <=2" % P.ambient_dim())
    #good = P.plot()._objects[-1]
    #assert isinstance(good, sage.plot.polygon.Polygon), "good not polygon"
    #return zip(good.xdata, good.ydata)
################################################
# Discret subsets of ZZ^d
################################################
class DiscreteSubset(SageObject):
    r"""
    A subset of ZZ^d.

    INPUT:

    - ``dimension`` -- integer, dimension of the space
    - ``predicate`` -- function ZZ^d -> {False, True} (default: ``None``)
    - ``edge_predicate`` -- function ZZ^d,ZZ^d -> {False, True} (default:
      ``None``)
    - ``iterator`` -- function (default: ``None``) returning an iterator of
      points, it must be consistent with the predicate
    - ``roots`` -- list (default: ``None``) of some elements in self. If
      ``iterator`` is not provided, it is used to iterate the elements throught
      connectedness.

    EXAMPLES::

        sage: from slabbe import DiscreteSubset
        sage: DiscreteSubset(dimension=3)
        Subset of ZZ^3

    ::

        sage: p = DiscreteSubset(dimension=3, predicate=lambda x:True)
        sage: p
        Subset of ZZ^3

    ::

        sage: fn = lambda p : p[0]+p[1]<p[2]
        sage: p = DiscreteSubset(dimension=3, predicate=fn, roots=[(0,0,1)])
        sage: p
        Subset of ZZ^3

    ::

        sage: F = lambda p: Integers(7)(2*p[0]+5*p[1])
        sage: edge_predicate = lambda p,s: F(s) < F(s)
        sage: D = DiscreteSubset(dimension=3, edge_predicate=edge_predicate)
        sage: D
        Subset of ZZ^3

    From a list::

        sage: L = [(0,0,0), (1,0,0), (2,0,0), (3,0,0)]
        sage: s = DiscreteSubset.from_subset(L)
        sage: s
        Subset of ZZ^3

    Providing a root may be necessary if zero (the origin) is not inside::

        sage: predicate = lambda (x,y) : 4 < x^2 + y^2 < 25
        sage: D = DiscreteSubset(dimension=2, predicate=predicate, roots=[(3,0)])

    TESTS:

    No edges go outside of the box::

        sage: from slabbe import DiscreteBox
        sage: B = DiscreteBox([-1,1],[-1,1])
        sage: len(list(B.edges_iterator()))
        12
        sage: sorted(B.edges_iterator())
        [((-1, -1), (-1, 0)), ((-1, -1), (0, -1)), ((-1, 0), (-1, 1)),
        ((-1, 0), (0, 0)), ((-1, 1), (0, 1)), ((0, -1), (0, 0)), ((0, -1),
        (1, -1)), ((0, 0), (0, 1)), ((0, 0), (1, 0)), ((0, 1), (1, 1)),
        ((1, -1), (1, 0)), ((1, 0), (1, 1))]
    """
    def __init__(self, dimension=3, predicate=None, edge_predicate=None,
            iterator=None, roots=None):
        r"""
        Constructor.

        See class:`DiscreteSubset` for complete documentation.

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: DiscreteSubset(dimension=3)
            Subset of ZZ^3
        """
        self._space = FreeModule(ZZ, dimension)
        if predicate is None:
            self._predicate = lambda p: len(p) == dimension
        else:
            self._predicate = predicate
        if edge_predicate is None:
            self._edge_predicate = lambda p,s: p in self and s in self
        else:
            self._edge_predicate = edge_predicate
        self._iterator = iterator
        if roots is None:
            self._roots = None
        else:
            self._roots = map(self._space, roots)
            for p in self._roots: p.set_immutable()

    @cached_method
    def roots(self):
        r"""
        Return the roots, i.e., a list of elements in self.

        It also makes sure the roots are in self and raises an error otherwise.

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: s = DiscreteSubset.from_subset([(0,0,0)])
            sage: s.roots()
            [(0, 0, 0)]

        ::

            sage: predicate = lambda (x,y) : 4 < x^2 + y^2 < 25
            sage: D = DiscreteSubset(dimension=2, predicate=predicate, roots=[(3,0)])
            sage: D.roots()
            [(3, 0)]

        TESTS::

            sage: predicate = lambda (x,y) : 4 < x^2 + y^2 < 25
            sage: D = DiscreteSubset(dimension=2, predicate=predicate, roots=[(2,0)])
            sage: D.roots()
            Traceback (most recent call last):
            ...
            ValueError: root element (=(2, 0)) provided at initialisation is not in self

        An error is raised if the roots are inconsistent::

            sage: s = DiscreteSubset.from_subset([])
            sage: s.roots()
            Traceback (most recent call last):
            ...
            ValueError: default element (=(0, 0, 0)) is not in self, please provide one at initialisation
        """
        if self._roots is None:
            element = self._space(0)
            element.set_immutable()
            if not element in self: 
                raise ValueError("default element (={}) is not in self, please provide "
                        "one at initialisation".format(element))
            self._roots = [element]
        else:
            for p in self._roots:
                if not p in self: 
                    raise ValueError("root element (={}) provided at"
                            " initialisation is not in self".format(p))
        return self._roots

    @classmethod
    def from_subset(self, subset):
        r"""
        Constructor from a finite subset.

        INPUT:

        - ``subset`` -- iterable of integer coordinate points

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: L = [(0,0,0), (1,0,0), (2,0,0), (3,0,0)]
            sage: s = DiscreteSubset.from_subset(L)
            sage: s
            Subset of ZZ^3
            sage: all(p in s for p in s)
            True

        Note that tuple or mutable vectors work fine::

            sage: (0,0,0) in s
            True
            sage: vector((0,0,0)) in s
            True

        TESTS::

            sage: DiscreteSubset.from_subset([])
            Subset of ZZ^3
        """
        try:
            dimension = len(next(iter(subset)))
        except StopIteration:
            dimension = 3
        space = FreeModule(ZZ, dimension)
        space_subset = set()
        for p in subset:
            p = space(p)
            p.set_immutable()
            space_subset.add(p)
        def predicate(p):
            p = space(p)
            p.set_immutable()
            return p in space_subset
        iterator = space_subset.__iter__
        return DiscreteSubset(dimension=dimension, predicate=predicate,
                iterator=iterator)

    def _repr_(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: DiscreteSubset(dimension=3)
            Subset of ZZ^3

        ::

            sage: DiscreteSubset(dimension=3, predicate=lambda x:True)
            Subset of ZZ^3

        """
        return "Subset of ZZ^{}".format(self.dimension())

    def dimension(self):
        r"""
        Returns the dimension of the ambiant space.

        OUTPUT:

            integer

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: d = DiscreteSubset(dimension=3)
            sage: d.dimension()
            3

        ::

            sage: from slabbe import DiscreteBox
            sage: p = DiscreteBox([0,3], [0,3], [0,3], [0,3])
            sage: p.dimension()
            4
        """
        return self._space.dimension()

    def __contains__(self, p):
        r"""
        Returns True if the point p is in self.

        INPUT:

        - ``p`` - point in the space

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: d = DiscreteSubset(dimension=3)
            sage: (0,0,0) in d
            True
            sage: (0,0,0,0) in d
            False

        ::

            sage: fn = lambda p : p[0]+p[1]<p[2]
            sage: p = DiscreteSubset(dimension=3, predicate=fn)
            sage: vector((1,2,4)) in p
            True
            sage: vector((1,2,2)) in p
            False
        """
        return self._predicate(p)

    def has_edge(self, p, s):
        r"""
        Returns whether it has the edge (p, s) where s-p is a canonical
        vector.

        INPUT:

        - ``p`` - point in the space
        - ``s`` - point in the space

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: F = lambda p: Integers(7)(2*p[0]+5*p[1])
            sage: edge_predicate = lambda p,s: F(p) < F(s)
            sage: D = DiscreteSubset(dimension=3, edge_predicate=edge_predicate)
            sage: D.has_edge(vector((0,0)),vector((1,0)))
            True
            sage: D.has_edge(vector((0,0)),vector((-1,0)))
            True
            sage: D.has_edge(vector((-1,1)),vector((1,0)))
            False
        """
        return self._edge_predicate(p, s)

    def an_element(self):
        r"""
        Returns an immutable element in self.

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: p = DiscreteSubset(dimension=3)
            sage: p.an_element()
            (0, 0, 0)
            sage: p.an_element().is_immutable()
            True

        ::

            sage: predicate = lambda (x,y) : 4 < x^2 + y^2 < 25
            sage: D = DiscreteSubset(dimension=2, predicate=predicate, roots=[(3,0)])
            sage: D.an_element()
            (3, 0)
        """
        return self.roots()[0]

    def __and__(self, other):
        r"""
        Return the intersection of self and other.

        INPUT:

        - ``other`` - Discrete object of dimension equal to s.

        EXAMPLES::

            sage: from slabbe import DiscreteSubset, DiscretePlane
            sage: d3 = DiscreteSubset(dimension=3)
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: p & d3
            Intersection of the following objects:
            Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8
            Subset of ZZ^3

        TESTS::

            sage: d3 = DiscreteSubset(dimension=3)
            sage: d3 & 4
            Traceback (most recent call last):
            ...
            TypeError: Impossible to construct an intersection containing the object : 4

        ::

            sage: d3 = DiscreteSubset(dimension=3)
            sage: d5 = DiscreteSubset(dimension=5)
            sage: d3 & d5
            Traceback (most recent call last):
            ...
            ValueError: Intersection not defined for objects not of the same dimension

        Intersection of intersection do not stack up::

            sage: d3 = DiscreteSubset(dimension=3)
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: I = p & d3
            sage: p & I
            Intersection of the following objects:
            Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8
            Subset of ZZ^3
            Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8
        """
        if isinstance(other, Intersection):
            return other.__and__(self)
        return Intersection((self, other))

    def __le__(self, other):
        r"""
        Return whether self is a subset of other.

        Runs into an infinite loop if self is infinite.

        EXAMPLES::

            sage: from slabbe import DiscretePlane, DiscreteTube
            sage: P = DiscretePlane([4,6,7], 17, mu=0)
            sage: tube = DiscreteTube([-6.4, 6.4], [-5.2, 5.2])
            sage: I = tube & P
            sage: I <= tube
            True
            sage: I <= P
            True
        """
        return all(p in other for p in self)

    def __ge__(self, other):
        r"""
        Return whether self is a superset of other.

        Runs into an infinite loop if other is infinite.

        EXAMPLES::

            sage: from slabbe import DiscretePlane, DiscreteTube
            sage: P = DiscretePlane([4,6,7], 17, mu=0)
            sage: tube = DiscreteTube([-6.4, 6.4], [-5.2, 5.2])
            sage: I = tube & P
            sage: tube >= I
            True
            sage: P >= I
            True
        """
        return all(p in self for p in other)

    @cached_method
    def base_edges(self):
        r"""
        Return a list of positive canonical vectors.

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: d = DiscreteSubset(dimension=2)
            sage: d.base_edges()
            [(1, 0), (0, 1)]

        ::

            sage: from slabbe import DiscretePlane
            sage: P = DiscretePlane([3,4,5], 12)
            sage: P.base_edges()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

        """
        d = self.dimension()
        canonical_base = [copy(self._space.zero()) for _ in range(d)]
        for i in range(d):
            canonical_base[i][i] = 1
        return canonical_base
        #return canonical_base + [-v for v in canonical_base]

    def children(self, p):
        r"""

        EXAMPLES::

            sage: from slabbe import DiscretePlane
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: list(p.children(vector((0,0,0))))
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        """
        for v in self.base_edges():
            u = p+v
            u.set_immutable()
            if u in self: yield u
            u = p-v
            u.set_immutable()
            if u in self: yield u

    def d_neighbors(self, p, d=2):
        r"""
        Retourne le voisinage du point p, i.e. les points parmi les 3^d
        possible qui appartiennent a l'objet discret.

        INPUT:

        - ``p`` - un point discret
        - ``d`` - integer (optional, default:2), 

        OUTPUT:

        liste de points

        EXAMPLES::

            sage: from slabbe import DiscretePlane
            sage: p = DiscretePlane([1,3,7], 10)
            sage: p.d_neighbors((0,0,0))
            [(-1, -1, 1), (-1, 0, 1), (-1, 1, 0), (-1, 1, 1), (0, -1, 1),
            (0, 0, 0), (0, 0, 1), (0, 1, 0), (1, -1, 1), (1, 0, 0), (1, 0,
            1), (1, 1, 0)]
        """
        if d != 2:
            raise NotImplementedError("implemented only for d=2")
        p = self._space(p)
        v = [-1, 0, 1]
        L = []
        for a in itertools.product(*[v]*self.dimension()):
            b = p + self._space(a)
            if b in self:
                L.append(b)
        return L

    def connected_component_iterator(self, roots=None):
        r"""
        Return an iterator over the connected component of the root.

        INPUT:

        - ``roots`` - list of some elements immutable in self

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: p = DiscreteSubset(dimension=3, roots=[(0,0,0)])
            sage: it = p.connected_component_iterator()
            sage: [next(it) for _ in range(5)]
            [(0, 0, 0), (1, 0, 0), (0, 0, 1), (0, 0, -1), (0, -1, 0)]

        ::

            sage: from slabbe import DiscretePlane
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: root = vector((0,0,0))
            sage: root.set_immutable()
            sage: it = p.connected_component_iterator(roots=[root])
            sage: [next(it) for _ in range(5)]
            [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 1, 0)]
        """
        roots = roots if roots else self.roots()
        if not all(root in self for root in roots):
            raise ValueError("roots (=%s) must all be in self(=%s)" % (roots, self))
        #for root in roots:
        #    root.set_immutable()
        C = RecursivelyEnumeratedSet(seeds=roots, successors=self.children, structure='symmetric')
        return C.breadth_first_search_iterator()

    def __iter__(self):
        r"""
        Return an iterator over self.

        EXAMPLES::

            sage: from slabbe import DiscretePlane
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: it = iter(p)
            sage: [next(it) for _ in range(5)]
            [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 1, 0)]

        ::

            sage: from slabbe import DiscreteSubset
            sage: L = [(0,0,0), (1,0,0), (2,0,0), (3,0,0)]
            sage: s = DiscreteSubset.from_subset(L)
            sage: sorted(s)
            [(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)]
            sage: all(p in s for p in s)
            True
        """
        if self._iterator is None:
            return self.connected_component_iterator()
        else:
            return self._iterator()

    def list(self):
        r"""
        Return the list of elements in self.

        EXAMPLES::

            sage: from slabbe import DiscretePlane, DiscreteTube
            sage: P = DiscretePlane([3,4,5], 12, mu=20)
            sage: tube = DiscreteTube([0,2],[0,2])
            sage: I = P & tube
            sage: sorted(I.list())
            [(-3, -1, -1), (-3, -1, 0), (-2, -2, -1), (-2, -2, 0), (-2, -1,
            -1), (-2, -1, 0), (-2, 0, -1), (-1, -1, -1)]
        """
        return list(self)

    def level_iterator(self):
        r"""
        This returns an iterator of the levels according to the given roots.

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: p = DiscreteSubset(dimension=3, roots=[(0,0,0)])
            sage: it = p.level_iterator()
            sage: sorted(next(it))
            [(0, 0, 0)]
            sage: sorted(next(it))
            [(-1, 0, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (0, 1, 0), (1, 0, 0)]
            sage: sorted(next(it))
            [(-2, 0, 0),
             (-1, -1, 0),
             (-1, 0, -1),
             (-1, 0, 1),
             (-1, 1, 0),
             (0, -2, 0),
             (0, -1, -1),
             (0, -1, 1),
             (0, 0, -2),
             (0, 0, 2),
             (0, 1, -1),
             (0, 1, 1),
             (0, 2, 0),
             (1, -1, 0),
             (1, 0, -1),
             (1, 0, 1),
             (1, 1, 0),
             (2, 0, 0)]

        ::

            sage: from slabbe import DiscretePlane
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: it = p.level_iterator()
            sage: sorted(next(it))
            [(0, 0, 0)]
            sage: sorted(next(it))
            [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
            sage: sorted(next(it))
            [(-1, 0, 1),
             (-1, 1, 0),
             (0, -1, 1),
             (0, 1, 1),
             (0, 2, 0),
             (1, 0, 1),
             (1, 1, 0),
             (2, 0, 0)]
        """
        roots = self.roots()
        t = RecursivelyEnumeratedSet(seeds=roots, successors=self.children, structure='symmetric')
        return t.graded_component_iterator()

    def edges_iterator(self):
        r"""
        Returns an iterator over the pair of points in self that are
        adjacents, i.e. their difference is a canonical vector.

        It considers only points that are connected to the given roots.

        INPUT:

        EXAMPLES::

            sage: from slabbe import DiscretePlane
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: it = p.edges_iterator()
            sage: next(it)
            ((0, 0, 0), (1, 0, 0))
            sage: next(it)
            ((0, 0, 0), (0, 1, 0))
            sage: next(it)
            ((0, 0, 0), (0, 0, 1))
            sage: next(it)
            ((-1, 1, 0), (0, 1, 0))
            sage: next(it)
            ((-2, 1, 0), (-1, 1, 0))
        """
        base_edges = self.base_edges()
        for p in self:
            # astuce pour eviter les doublons pour cette connexite
            if sum(p)%2 != 0: continue
            for e in base_edges:
                s = p+e
                if self.has_edge(p, s):
                    yield (p, s)
                s = p-e
                if self.has_edge(s, p):
                    yield (s, p)

    def projection_matrix(self, m='isometric', oblique=None):
        r"""
        Return a projection matrix.

        INPUT:

        - ``m`` -- projection matrix (default: ``'isometric'``), it can be one of
          the following:
        
          - ``'isometric'`` - the isometric projection is used by default
          - matrix - a 2 x 3 matrix
          - ``'belle'`` - shortcut for ``matrix(2, [1/3.0, 1, 0, 2/3.0, 0, 1])``
          - vector - defines the projection on the plane orthogonal to the
            vector.

        - ``oblique`` -- vector (default: ``None``), vector perpendicular
          to the range space

        EXAMPLES::

            sage: from slabbe import DiscreteSubset
            sage: d = DiscreteSubset(dimension=3)
            sage: d.projection_matrix(vector((2,3,4))) # tolerance 0.00001
            [  1.00000000000000  0.000000000000000 -0.500000000000000]
            [ 0.000000000000000   1.00000000000000 -0.750000000000000]
            sage: d.projection_matrix((2,3,4)) # tolerance 0.00001
            [  1.00000000000000  0.000000000000000 -0.500000000000000]
            [ 0.000000000000000   1.00000000000000 -0.750000000000000]
            sage: d.projection_matrix()         # tolerance 0.00001
            [-0.866025403784  0.866025403784             0.0]
            [           -0.5            -0.5             1.0]
            sage: d.projection_matrix(_) # tolerance 0.00001
            [-0.866025403784439  0.866025403784439  0.000000000000000]   
            [-0.500000000000000 -0.500000000000000   1.00000000000000]   
            sage: d.projection_matrix('belle')    # tolerance 0.00001
            [0.333333333333            1.0            0.0]
            [0.666666666667            0.0            1.0]
        """
        if m == 'isometric':
            return matrix(2, [-1.7320508075688772*0.5, 1.7320508075688772*0.5, 0, 
                                           -0.5, -0.5, 1])
        elif m == 'belle':
            return matrix(2, [1/3.0, 1, 0, 2/3.0, 0, 1])
        elif m == 'plusbelle':
            return matrix(2, [1, 0, -0.4, 0, 1, -0.4])
        else:
            try:
                v = matrix(RR, m)
            except:
                raise
                #TypeError, "unknown input for m(=%s)" % m
            dim = v.ncols(), v.nrows()
            if dim == (3,2):
                return v
            elif dim == (3,1):
                return v.right_kernel().basis_matrix()
            else:
                raise ValueError, "uncorrect dimension (=%s) " % (dim,)

    def plot_points(self, color='blue', m=None):
        r"""
        Returns a 2d or 3d graphics object of the points of self.

        INPUT:

        - ``color`` -- string (default: ``'blue'``), the color of the points
        - ``m`` -- projection matrix (default: ``None``), it can be one of
          the following:
        
          - ``None`` - no projection is done
          - ``'isometric'`` - the isometric projection
          - matrix - a 2 x n projection matrix
          - ``'belle'`` - shortcut for ``matrix(2, [1/3.0, 1, 0, 2/3.0, 0, 1])``
          - vector - defines the projection on the plane orthogonal to the
            vector.

        EXAMPLES:

        A 2d plot of a 2d object::

            sage: from slabbe import DiscreteSubset, DiscreteBox
            sage: D = DiscreteSubset(dimension=2)
            sage: box = DiscreteBox([-5,5],[-5,5])
            sage: I = D & box
            sage: I.plot_points(color='green')      # optional long

        A 3d plot of a 3d object::

            sage: D = DiscreteSubset(dimension=3)
            sage: box = DiscreteBox([-5,5],[-5,5],[-5,5])
            sage: I = D & box
            sage: I.plot_points(color='green')      # optional long

        A 2d plot of a 3d object::

            sage: D = DiscreteSubset(dimension=3)
            sage: box = DiscreteBox([-5,5],[-5,5],[-5,5])
            sage: I = D & box
            sage: I.plot_points(color='green', m='isometric')      # optional long
        """
        if m is None:
            L = self.list()
        else:
            m = self.projection_matrix(m=m)
            L = [m*p for p in self]
        if not all(len(p) in [2,3] for p in L[:5]):
            raise ValueError("dimension of points after projection must be 2 or 3")
        G = point(L, color=color)
        #G.set_aspect_ratio(1) # for 2d only
        #G.axes(False)         # for 2d only
        return G

    def plot_points_at_distance(self, k, color='blue', projmat=None):
        r"""
        Plot points at distance k from the roots.

        INPUT:

        - ``k`` - integer

        EXAMPLES::

            sage: alpha = solve(x+x**2+x**3==1, x)[2].right()
            sage: vv = vector((alpha, alpha+alpha**2, 1))
            sage: omega = (1+alpha)**2 / 2
            sage: from slabbe import DiscretePlane
            sage: Pr = DiscretePlane(vv, omega, mu=pi, prec=200)
            sage: Pr.plot_points_at_distance(200)               # optional long
            sage: Pr.plot_points_at_distance(200, projmat='isometric') # optional long
        """
        it = self.level_iterator()
        G = Graphics()
        for i in range(k):
            B = next(it)
        if projmat is None:
            L = [p for p in B]
        else:
            m = self.projection_matrix(m=projmat)
            L = [m * p for p in B]
            roots = [m * root for root in self.roots()]
        G += point(L, color='black')
        for root in roots:
            G += circle(root, 0.5, color='red', thickness=4)
        return G

    def plot_edges(self, color='blue', m=None):
        r"""
        Returns the mesh of the plane. The mesh is the union of segments
        joining two adjacents points.

        INPUT:

        - ``color`` -- string (default: ``'blue'``), the color of the edges
        - ``m`` -- projection matrix (default: ``None``), it can be one of
          the following:
        
          - ``None`` - no projection is done
          - ``'isometric'`` - the isometric projection
          - matrix - a 2 x 3 matrix
          - ``'belle'`` - shortcut for ``matrix(2, [1/3.0, 1, 0, 2/3.0, 0, 1])``
          - vector - defines the projection on the plane orthogonal to the
            vector.

        EXAMPLES:

        A 2d plot of a 2d object::

            sage: from slabbe import DiscreteSubset, DiscreteBox
            sage: D = DiscreteSubset(dimension=2)
            sage: box = DiscreteBox([-5,5],[-5,5])
            sage: I = D & box
            sage: I.plot_edges(color='green') # optional long

        A 3d plot of a 3d object::

            sage: D = DiscreteSubset(dimension=3)
            sage: box = DiscreteBox([-3,3],[-3,3],[-3,3])
            sage: I = D & box
            sage: I.plot_edges(color='green') # optional long

        A 2d plot of a 3d object::

            sage: D = DiscreteSubset(dimension=3)
            sage: box = DiscreteBox([-3,3],[-3,3],[-3,3])
            sage: I = D & box
            sage: I.plot_edges(color='green', m='isometric') # optional long
        """
        if m is None:
            L = self.edges_iterator()
        else:
            m = self.projection_matrix(m=m)
            L = [(m*p,m*q) for (p,q) in self.edges_iterator()]
        G = Graphics()
        for edge in L:
                G += line(edge, color=color)
        # G.axes(False) # for 2d only
        return G
    def plot_cubes(self, **kwds):
        r"""
        Returns the discrete object as cubes in 3d.

        EXAMPLES::

            sage: from slabbe import DiscretePlane, DiscreteTube
            sage: P = DiscretePlane([3,4,5], 12, mu=20)
            sage: tube = DiscreteTube([-5,5],[-5,5])
            sage: I = P & tube
            sage: I.plot_cubes(color='red', frame_thickness=1 # optional long)

        TESTS::

            sage: from slabbe import DiscreteBox
            sage: box = DiscreteBox([-5,5],[-5,5])
            sage: box.plot_cubes() # optional long
            Traceback (most recent call last):
            ...
            ValueError: this method is currently implemented only for objects living in 3 dimensions

        """
        if self.dimension() != 3:
            raise ValueError("this method is currently implemented only for "
                             "objects living in 3 dimensions")
        p = sum(cube(p, **kwds) for p in self)
        return p

    def plot(self, frame=False, edgecolor="blue", pointcolor="blue"):
        r"""
        Return a plot (2d or 3d) of the points and edges of self.

        INPUT:

        - ``frame`` - (default: False) if True, draw a bounding frame with
          labels
        - ``edgecolor`` -- string (default: ``'blue'``), the color of the edges
        - ``pointcolor`` -- string (default: ``'blue'``), the color of the points

        EXAMPLES:

        2d example::

            sage: from slabbe import DiscreteBox
            sage: box = DiscreteBox([-5,5],[-5,5])
            sage: box.plot() # optional long

        3d example::

            sage: from slabbe import DiscretePlane, DiscreteTube
            sage: P = DiscretePlane([1,3,7], 11)
            sage: tube = DiscreteTube([-5,5],[-5,5])
            sage: I = P & tube
            sage: I.plot() # optional long
        """
        G = self.plot_edges(color=edgecolor) 
        G += self.plot_points(color=pointcolor)
        if self.dimension() == 2:
            G += text(repr(self), (0.5,1.1), axis_coords=True)
            G.axes(False)
            G.set_aspect_ratio(1)
        elif self.dimension() == 3:
            G += text3d(repr(self), (-2,-2,-2))
        else:
            raise NotImplementedError("implemented only for dimension(=%s) 2 or 3" % self.dimension())
        return G

    def tikz_projection_scale(self, projmat='isometric', scale=1, extra=""):
        r"""
        INPUT:

        - ``projmat`` -- (default: ``'isometric'``) It can be one of the following:

          - ``'isometric'`` - the isometric projection is used by default
          - matrix - a 2 x 3 matrix
          - ``'belle'`` - shortcut for ``matrix(2, [1/3.0, 1, 0, 2/3.0, 0, 1])``
          - vector - defines the projection on the plane orthogonal to the
            vector.

        - ``scale`` -- real number (default: 1), scaling constant for the
          whole figure

        - ``extra`` -- string (default: ``''``)

        EXAMPLES::

            sage: from slabbe import DiscretePlane
            sage: p = DiscretePlane([1,3,7], 11)
            sage: p.tikz_projection_scale()
            [x={(-0.866025cm,-0.500000cm)}, y={(0.866025cm,-0.500000cm)},
            z={(0.000000cm,1.000000cm)}, scale=1]
            sage: p.tikz_projection_scale(extra="xshift=4cm")
            [x={(-0.866025cm,-0.500000cm)}, y={(0.866025cm,-0.500000cm)},
            z={(0.000000cm,1.000000cm)}, scale=1,xshift=4cm]

        """
        if self.dimension() != 3:
            raise NotImplementedError("implemented only for dimension (=%s) 3" % self.dimension())

        projmat = self.projection_matrix(projmat)

        e1 = projmat*vector([1,0,0])
        e2 = projmat*vector([0,1,0])
        e3 = projmat*vector([0,0,1])

        t = (e1[0], e1[1], e2[0], e2[1], e3[0], e3[1], scale)
        s = 'x={(%fcm,%fcm)}, y={(%fcm,%fcm)}, \nz={(%fcm,%fcm)}, scale=%s' % t
        if extra:
            s += "," + extra
        return LatexExpr("[%s]\n" % s)

    def tikz_axes(self, xshift=0, yshift=0, label="e", projmat='isometric'):
        r"""
        Return the tikz code for drawing axes.

        INPUT:

        - ``xshift`` - integer (default: ``0``), x shift
        - ``yshift`` - integer (default: ``0``), y shift
        - ``label`` - string (default: ``"e"``), label for base vectors
        - ``projmat`` - matrix (default: ``'isometric'``), projection matrix

        OUTPUT:

        string

        EXAMPLES:

        2d example::

            sage: from slabbe import DiscreteSubset
            sage: d = DiscreteSubset(dimension=2)
            sage: d.tikz_axes()
            %the axes
            \begin{scope}[xshift=0cm,yshift=0cm]
            \draw[->,>=latex, very thick, blue] (0,0) -- (1, 0);
            \draw[->,>=latex, very thick, blue] (0,0) -- (0, 1);
            \node at (1.40000000000000,0) {$e_1$};
            \node at (0,1.40000000000000) {$e_2$};
            \end{scope}

        3d example::

            sage: d = DiscreteSubset(dimension=3)
            sage: d.tikz_axes(projmat='isometric')
            %the axes
            \begin{scope}
            [x={(-0.866025cm,-0.500000cm)}, y={(0.866025cm,-0.500000cm)}, 
            z={(0.000000cm,1.000000cm)}, scale=1,xshift=0,yshift=0]
            \draw[fill=white] (2,0,0) rectangle (-1.8,.1,1);
            \draw[->,>=latex, very thick, blue] (0,0,0) -- (1, 0, 0);
            \draw[->,>=latex, very thick, blue] (0,0,0) -- (0, 1, 0);
            \draw[->,>=latex, very thick, blue] (0,0,0) -- (0, 0, 1);
            \node at (1.40000000000000,0,0) {$e_1$};
            \node at (0,1.40000000000000,0) {$e_2$};
            \node at (0,0,1.40000000000000) {$e_3$};
            \end{scope}

        """
        if self.dimension() == 2:
            s = "%the axes\n"
            s += "\\begin{scope}[xshift=%scm,yshift=%scm]\n" % (xshift, yshift)
            s += "\\draw[->,>=latex, very thick, blue] (0,0) -- (1, 0);\n"
            s += "\\draw[->,>=latex, very thick, blue] (0,0) -- (0, 1);\n"
            s += "\\node at (1.40000000000000,0) {$%s_1$};\n" % label
            s += "\\node at (0,1.40000000000000) {$%s_2$};\n" % label
            s += "\\end{scope}\n"
            return LatexExpr(s)
        elif self.dimension() == 3:
            s = "%the axes\n"
            s += "\\begin{scope}\n"
            extra = "xshift=%s,yshift=%s" % (xshift, yshift)
            s += self.tikz_projection_scale(projmat, scale=1, extra=extra)
            s += "\\draw[fill=white] (2,0,0) rectangle (-1.8,.1,1);\n"
            s += "\\draw[->,>=latex, very thick, blue] (0,0,0) -- (1, 0, 0);\n"
            s += "\\draw[->,>=latex, very thick, blue] (0,0,0) -- (0, 1, 0);\n"
            s += "\\draw[->,>=latex, very thick, blue] (0,0,0) -- (0, 0, 1);\n"
            s += "\\node at (1.40000000000000,0,0) {$%s_1$};\n" % label
            s += "\\node at (0,1.40000000000000,0) {$%s_2$};\n" % label
            s += "\\node at (0,0,1.40000000000000) {$%s_3$};\n" % label
            s += "\\end{scope}\n"
            return LatexExpr(s)
        else:
            raise ValueError("this method is currently implemented only for "
                             "objects living in 2 or 3 dimensions")

    def tikz_edges(self, style='very thick', color='blue', projmat=None):
        r"""
        Returns the mesh of the object. The mesh is the union of segments
        joining two adjacents points.

        INPUT:

        - ``style`` - string (default: ``'dashed, very thick'``)
        - ``color`` - string or callable (default: ``'blue'``), the color
          of all edges or a function : (u,v) -> color of the edge (u,v)
        - ``projmat`` - matrix (default: ``None``), projection matrix, if
          None, no projection is done.

        EXAMPLES::

            sage: from slabbe import DiscretePlane
            sage: p = DiscretePlane([2,3,5], 4)
            sage: p.tikz_edges()
            \draw[very thick, blue] (0, 0, 0) -- (1, 0, 0);
            \draw[very thick, blue] (0, 0, 0) -- (0, 1, 0);
            \draw[very thick, blue] (-1, 1, 0) -- (0, 1, 0);

        ::

            sage: p.tikz_edges(color='orange')
            \draw[very thick, orange] (0, 0, 0) -- (1, 0, 0);
            \draw[very thick, orange] (0, 0, 0) -- (0, 1, 0);
            \draw[very thick, orange] (-1, 1, 0) -- (0, 1, 0);

        ::

            sage: c = lambda u,v: 'red' if u == 0 else 'blue'
            sage: p.tikz_edges(color=c)
            \draw[very thick, red] (0, 0, 0) -- (1, 0, 0);
            \draw[very thick, red] (0, 0, 0) -- (0, 1, 0);
            \draw[very thick, blue] (-1, 1, 0) -- (0, 1, 0);

        ::

            sage: from slabbe.discrete_subset import M3to2
            sage: p.tikz_edges(projmat=M3to2)
            \draw[very thick, blue] (0.00000, 0.00000) -- (-0.86603, -0.50000);
            \draw[very thick, blue] (0.00000, 0.00000) -- (0.86603, -0.50000);
            \draw[very thick, blue] (1.73205, 0.00000) -- (0.86603, -0.50000);

        """
        s = ''
        if isinstance(color, str):
            color_f = lambda u,v: color
        else:
            color_f = color
        for (p, q) in self.edges_iterator():
            c = color_f(p,q)
            if projmat:
                p = projmat * p
                q = projmat * q
                p = "(" + ", ".join("%.5f" % a for a in p) + ")"
                q = "(" + ", ".join("%.5f" % a for a in q) + ")"
            s += "\\draw[%s, %s] %s -- %s;\n" % (style, c,  p, q)
        return LatexExpr(s)

    def tikz_points(self, size='0.8mm', label=None, label_pos='right',
            fill='black', options="", filter=None,
            projmat=None):
        r"""
        INPUT:

        - ``size`` - string (default: ``'0.8mm'``), size of the
          points
        - ``label`` - function (default: ``None``), print some label next
          to the point
        - ``label_pos`` - function (default: ``'right'``), tikz label
          position
        - ``fill`` - string (default: ``'black'``), fill color
        - ``options`` - string (default: ``''``), author tikz node circle options
        - ``filter`` - boolean function, if filter(p) is False, the point p
          is not drawn
        - ``projmat`` - matrix (default: ``None``), projection matrix, if
          None, no projection is done.

        EXAMPLES::

            sage: from slabbe import DiscreteBox
            sage: p = DiscreteBox([0,3], [0,3], [0,3])
            sage: s = p.tikz_points()
            sage: lines = s.splitlines()
            sage: lines[0]
            '\\node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (...) {};'

        ::

            sage: from slabbe import DiscretePlane, DiscreteTube
            sage: p = DiscretePlane([1,3,7], 11)
            sage: d = DiscreteTube([-1,1],[-1,1])
            sage: I = p & d
            sage: I.tikz_points()
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (0, 0, 0) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (1, 0, 0) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (0, 1, 0) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (0, 0, 1) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (0, 1, 1) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (1, 1, 0) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (1, 0, 1) {};


        Using a filter on the points::

            sage: I.tikz_points(filter=lambda x:sum(x)==1)
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (1, 0, 0) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (0, 1, 0) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (0, 0, 1) {};

        Of a finite subset::

            sage: from slabbe import DiscreteSubset
            sage: V = [(0,0,0), (1,1,0), (1,-1,1), (-2,1,0), (2,0,1), (-1,2,0),
            ....:      (-1,0,1), (0,1,1)]
            sage: s = DiscreteSubset.from_subset(V)
            sage: s.tikz_points()
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (-2, 1, 0) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (1, 1, 0) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (0, 1, 1) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (-1, 2, 0) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (-1, 0, 1) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (0, 0, 0) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (1, -1, 1) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (2, 0, 1) {};
        """
        it = iter(self)
        if filter:
            it = itertools.ifilter(filter, it)
        s = ''
        if projmat:
            for pt in it:
                pt2d = projmat * pt
                pt2d = "(%.5f, %.5f)" % tuple(pt2d)
                s += '\\node[circle,fill=%s,draw=black,' % fill
                s += 'minimum size=%s,inner sep=0pt,%s] at %s {};\n' % (size, options, pt2d)
                if label:
                    s += '\\node[%s] at %s {$%s$};\n' % (label_pos, pt2d, label(pt))
        else:
            for pt in it:
                s += '\\node[circle,fill=%s,draw=black,' % fill
                s += 'minimum size=%s,inner sep=0pt,%s] at %s {};\n' % (size, options, pt)
                if label:
                    s += '\\node[%s] at %s {$%s$};\n' % (label_pos, pt, label(pt))
        return LatexExpr(s)

    def tikz_noprojection(self, projmat=None, scale=1, clip=[], edges=True,
            points=True, axes=False, point_kwds={}, edge_kwds={}, axes_kwds={}, extra_code=''):
        r"""
        Return the tikz code of self.
        
        In this version, the points are not projected. If the points are in
        3d, the tikz 3d picture is used.

        INPUT:

        - ``projmat`` -- (default: None) 2*3 projection matrix
          for drawing unit faces, the isometric projection is used by
          default
        - ``scale`` -- real number (default: 1), scaling constant for the
          whole figure
        - ``clip`` - list (default: ``[]``), list of points describing a
          cliping path once projected. Works only if ``self.dimension()``
          is 2.
        - ``edges`` - bool (default: ``True``), whether to draw edges
        - ``points`` - bool (default: ``True``), whether to draw points
        - ``axes`` - bool (default: ``False``), whether to draw axes
        - ``point_kwds`` - dict (default: ``{}``)
        - ``edge_kwds`` - dict (default: ``{}``)
        - ``axes_kwds`` - dict (default: ``{}``)
        - ``extra_code`` -- string (default: ``''``), extra tikz code to add

        EXAMPLES:

        Object in 2d::

            sage: from slabbe import DiscreteLine, DiscreteBox
            sage: L = DiscreteLine([2,5], 2+5, mu=0)
            sage: b = DiscreteBox([-5,5],[-5,5])
            sage: I = L & b
            sage: point_kwds = {'label':lambda p:2*p[0]+5*p[1],'label_pos':'above right'}
            sage: tikz = I.tikz_noprojection(scale=0.5,point_kwds=point_kwds)
            sage: tikz
            \documentclass[tikz]{standalone}
            \usepackage{amsmath}
            \begin{document}
            \begin{tikzpicture}
            [scale=0.500000000000000]
            \draw[very thick, blue] (0, 0) -- (1, 0);
            \draw[very thick, blue] (0, 0) -- (0, 1);
            \draw[very thick, blue] (2, 0) -- (3, 0);
            ...
            ... 40 lines not printed (2659 characters in total) ...
            ...
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (-5, 2) {};
            \node[above right] at (-5, 2) {$0$};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (-5, 3) {};
            \node[above right] at (-5, 3) {$5$};
            \end{tikzpicture}
            \end{document}

        Object in 3d::

            sage: from slabbe import DiscretePlane, DiscreteTube
            sage: p = DiscretePlane([1,3,7], 11)
            sage: d = DiscreteTube([-5,5],[-5,5])
            sage: I = p & d
            sage: s = I.tikz_noprojection()
            sage: s
            \documentclass[tikz]{standalone}
            \usepackage{amsmath}
            \begin{document}
            \begin{tikzpicture}
             [x={(-0.866025cm,-0.500000cm)}, y={(0.866025cm,-0.500000cm)},
            z={(0.000000cm,1.000000cm)}, scale=1]
            \draw[very thick, blue] (0, 0, 0) -- (1, 0, 0);
            \draw[very thick, blue] (0, 0, 0) -- (0, 1, 0);
            ...
            ... 311 lines not printed (20339 characters in total) ...
            ...
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (1, -4, 3) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (4, 4, -1) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (5, 3, -1) {};
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (6, 2, -1) {};
            \end{tikzpicture}
            \end{document}
        """
        s = '\\begin{tikzpicture}\n'
        if self.dimension() == 2:
            s += "[scale=%s]\n" % scale
        elif self.dimension() == 3:
            if projmat is None:
                projmat = 'isometric'
            s += self.tikz_projection_scale(projmat, scale)
        else:
            raise NotImplementedError("implemented only for dimension(=%s) 2 or 3" % self.dimension())
        if axes:
            s += self.tikz_axes(**axes_kwds)
        if clip and self.dimension() == 2:
            ss = " -- ".join(str(pt) for pt in clip)
            #s += "\\draw %s -- cycle;" % ss
            s += "\\clip %s -- cycle;" % ss
        if edges:
            s += self.tikz_edges(**edge_kwds)
        if points:
            s += self.tikz_points(**point_kwds)
        s += extra_code
        s += '\\end{tikzpicture}\n'
        return TikzPicture(s)

    def tikz(self, projmat=M3to2, scale=1, clip=[], contour=[],
            edges=True, points=True, axes=False, point_kwds={},
            edge_kwds={}, axes_kwds={}, extra_code=''):
        r"""
        INPUT:

        - ``projmat`` -- (default: M3to2) 2 x dim projection
          matrix where dim is the dimensoin of self, the isometric
          projection is used by default
        - ``scale`` -- real number (default: 1), scaling constant for the
          whole figure
        - ``clip`` - list (default: ``[]``), list of points whose convex
          hull describes a cliping path
        - ``contour`` - list (default: ``[]``), list of points describing a
          contour path to be drawn
        - ``edges`` - bool (default: ``True``), whether to draw edges
        - ``points`` - bool (default: ``True``), whether to draw points
        - ``axes`` - bool (default: ``False``), whether to draw axes
        - ``point_kwds`` - dict (default: ``{}``)
        - ``edge_kwds`` - dict (default: ``{}``)
        - ``axes_kwds`` - dict (default: ``{}``)
        - ``extra_code`` -- string (default: ``''``), extra tikz code to add

        EXAMPLES::

            sage: from slabbe import DiscretePlane
            sage: p = DiscretePlane([2,3,5], 10)
            sage: p.tikz(points=False, edges=False)
            \documentclass[tikz]{standalone}
            \usepackage{amsmath}
            \begin{document}
            \begin{tikzpicture}
            [scale=1]
            \end{tikzpicture}
            \end{document}
        """
        s = '\\begin{tikzpicture}\n'
        s += '[scale=%s]\n' % scale
        #clip
        if clip:
            if len(clip[0]) == 2:
                pass
            elif len(clip[0]) == projmat.ncols() and projmat.nrows() == 2:
                clip = [projmat * c for c in clip]
            else:
                raise ValueError("incorrect space dimension of clip data (=%s)" %len(clip[0]))
            clip_convex = convex_boundary(clip)
            clip_s = " -- ".join(str(pt) for pt in clip_convex)
            s += "\\clip %s -- cycle;\n" % clip_s
        # edges
        if edges:
            s += self.tikz_edges(projmat=projmat,**edge_kwds)
        # points
        if points:
            s += self.tikz_points(projmat=projmat,**point_kwds)
        # contour
        if contour:
            if 2 != len(contour[0]) == projmat.ncols() and projmat.nrows() == 2:
                contour = [projmat * c for c in contour]
            else:
                raise ValueError("incorrect space dimension of contour data (=%s)" %len(contour[0]))
            contour_convex = convex_boundary(contour)
            contour_s = " -- ".join(str(pt) for pt in contour_convex)
            if clip:
                s += "\\filldraw[fill=white,very thick,dotted,opacity=0.5,even odd rule]\n"
                s += "  %s -- cycle\n" % clip_s
                s += "  %s -- cycle;\n" % contour_s
            else:
                s += "\\draw %s -- cycle;\n" % contour_s
        s += extra_code
        # axes
        if axes:
            s += self.tikz_axes(projmat=projmat,**axes_kwds)
        s += '\\end{tikzpicture}\n'
        return TikzPicture(s)

class Intersection(DiscreteSubset):
    r"""
    Intersection

    todo:

    - Rendre l'heritage 3d automatique

    INPUT:

    - ``objets`` - un tuple d'objets discrets

    EXAMPLES:

    Intersection de deux plans::

        sage: from slabbe import DiscretePlane, Intersection
        sage: p = DiscretePlane([1,3,7],11)
        sage: q = DiscretePlane([1,3,5],9)
        sage: Intersection((p,q))
        Intersection of the following objects:
        Set of points x in ZZ^3 satisfying: 0 <= (1, 3, 7) . x + 0 < 11
        Set of points x in ZZ^3 satisfying: 0 <= (1, 3, 5) . x + 0 < 9

    Shortcut::

        sage: p & q
        Intersection of the following objects:
        Set of points x in ZZ^3 satisfying: 0 <= (1, 3, 7) . x + 0 < 11
        Set of points x in ZZ^3 satisfying: 0 <= (1, 3, 5) . x + 0 < 9

    Intersection of a plane and a tube::

        sage: from slabbe import DiscreteTube
        sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
        sage: d = DiscreteTube([-5,5],[-5,5])
        sage: I = p & d
        sage: I
        Intersection of the following objects:
        Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8
        DiscreteTube: Preimage of [-5, 5] x [-5, 5] by a 2 by 3 matrix
        sage: len(list(I))
        115

    Intersection of a line and a box::

        sage: from slabbe import DiscreteLine, DiscreteBox
        sage: L = DiscreteLine([2,5], 2+5, mu=0)
        sage: b = DiscreteBox([-5,5],[-5,5])
        sage: I = L & b
        sage: I
        Intersection of the following objects:
        Set of points x in ZZ^2 satisfying: 0 <= (2, 5) . x + 0 < 7
        [-5, 5] x [-5, 5]

    TESTS:

    Intersected objects must be of the same dimension::

        sage: box = DiscreteBox([-5,5],[-5,5])
        sage: p = DiscretePlane([1,pi,7], 1+pi+7)
        sage: p & box
        Traceback (most recent call last):
        ...
        ValueError: Intersection not defined for objects not of the same dimension

    """
    def __init__(self, objets):
        r"""
        Constructeur d'intersection d'objets.

        EXAMPLES::

            sage: from slabbe import DiscretePlane, DiscreteSubset
            sage: d3 = DiscreteSubset(dimension=3)
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: I = p & d3
            sage: I
            Intersection of the following objects:
            Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8
            Subset of ZZ^3

        TESTS::

            sage: from slabbe import Intersection
            sage: type(Intersection([DiscreteSubset(dimension=2)]))
            <class 'slabbe.discrete_subset.Intersection'>
            sage: type(Intersection([DiscreteSubset(dimension=3)]))
            <class 'slabbe.discrete_subset.Intersection'>
            sage: type(Intersection([DiscreteSubset(dimension=4)]))
            <class 'slabbe.discrete_subset.Intersection'>
        """
        for o in objets:
            if not isinstance(o, DiscreteSubset):
                raise TypeError("Impossible to construct an intersection "
                        "containing the object : {}".format(o))
        dimension = objets[0].dimension()
        if not all(o.dimension() == dimension for o in objets):
            raise ValueError("Intersection not defined for objects not of the same dimension")
        self._objets = objets
        DiscreteSubset.__init__(self, dimension=dimension)
    @cached_method
    def roots(self):
        r"""
        EXAMPLES::

            sage: from slabbe import DiscreteBox, DiscreteSubset
            sage: d3 = DiscreteSubset(dimension=3, roots=[(0,0,0), (1,1,1)])
            sage: box = DiscreteBox([-5,5],[-5,5],[-5,5])
            sage: I = d3 & box
            sage: sorted(d3.roots())
            [(0, 0, 0), (1, 1, 1)]
            sage: box.roots()
            [(0, 0, 0)]
            sage: sorted(I.roots())
            [(0, 0, 0), (1, 1, 1)]
        """
        s = set()
        for o in self._objets:
            s.update(o.roots())
        return [p for p in s if p in self]

    def _repr_(self):
        r"""
        Retourne la représentation de l'objet en chaîne de caractères.

        EXAMPLES::

            sage: from slabbe import DiscretePlane, Intersection
            sage: p = DiscretePlane([1,3,7], 11)
            sage: q = DiscretePlane([1,3,5], 9)
            sage: Intersection((p,q))
            Intersection of the following objects:
            Set of points x in ZZ^3 satisfying: 0 <= (1, 3, 7) . x + 0 < 11
            Set of points x in ZZ^3 satisfying: 0 <= (1, 3, 5) . x + 0 < 9
        """
        s = "Intersection of the following objects:\n"
        s += '\n'.join(map(str, self._objets))
        return s

    def __contains__(self, p):
        r"""
        Returns True if the point p is in self.

        INPUT:

        - ``p`` - point in the space

        EXAMPLES::

            sage: from slabbe import DiscretePlane, DiscreteSubset
            sage: d3 = DiscreteSubset(dimension=3)
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: I = p & d3
            sage: vector((0,0,0)) in I
            True
        """
        return all(p in o for o in self._objets)

    def has_edge(self, p, s):
        r"""
        Returns whether it has the edge (p, s) where s-p is a canonical
        vector.

        INPUT:

        - ``p`` - point in the space
        - ``s`` - point in the space

        EXAMPLES::

            sage: from slabbe import DiscretePlane, DiscreteSubset
            sage: d3 = DiscreteSubset(dimension=3)
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: I = p & d3
            sage: I.has_edge(vector((0,0,0)),vector((0,0,1)))
            True
            sage: I.has_edge(vector((0,0,0)),vector((0,0,-1)))
            False

        TESTS::

            sage: from slabbe import DiscreteBox
            sage: F = lambda p: (2*p[0]+5*p[1]) % 7
            sage: edge_predicate = lambda p,s: F(p) < F(s)
            sage: D = DiscreteSubset(dimension=2, edge_predicate=edge_predicate)
            sage: b = DiscreteBox([-5,5],[-5,5])
            sage: I = D & b
            sage: all(I.has_edge(a,b) for a,b in I.edges_iterator())
            True
            sage: all(D.has_edge(a,b) for a,b in I.edges_iterator())
            True

        ::

            sage: from slabbe import ChristoffelGraph
            sage: C = ChristoffelGraph((2,5))
            sage: b = DiscreteBox([-5,5],[-5,5])
            sage: I = C & b
            sage: all(I.has_edge(a,b) for a,b in I.edges_iterator())
            True
            sage: all(C.has_edge(a,b) for a,b in I.edges_iterator())
            True

        """
        return all(o.has_edge(p,s) for o in self._objets)

    def __and__(self, other):
        r"""
        Return the intersection of self and other.

        INPUT:

        - ``other`` - Discrete object of dimension equal to s.

        EXAMPLES::

            sage: from slabbe import DiscretePlane, DiscreteSubset
            sage: d3 = DiscreteSubset(dimension=3)
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: I = p & d3
            sage: I & p
            Intersection of the following objects:
            Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8
            Subset of ZZ^3
            Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8

        ::

            sage: d3 = DiscreteSubset(dimension=3)
            sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
            sage: I = p & d3
            sage: p & I
            Intersection of the following objects:
            Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8
            Subset of ZZ^3
            Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8
        """
        if not isinstance(other, DiscreteSubset):
            raise TypeError("other(=%s) must be an instance of DiscreteSubset" % other)
        L = list(self._objets)
        if isinstance(other, Intersection):
            L.extend(other._objets)
        else:
            L.append(other)
        return Intersection(L)

    def an_element(self):
        r"""
        Returns an element in self.

        EXAMPLES::

            sage: from slabbe import DiscretePlane, DiscreteTube
            sage: P = DiscretePlane([4,6,7], 17, mu=0)
            sage: tube = DiscreteTube([-6.4, 6.4], [-5.2, 5.2])
            sage: I = tube & P
            sage: I.an_element()
            (0, 0, 0)
            sage: I.an_element() in I
            True

        TESTS::

            sage: P = DiscretePlane([4,6,7], 17, mu=0)
            sage: def contain(p): return 0 < P._v.dot_product(p) + P._mu <= P._omega
            sage: P._predicate = contain
            sage: tube = DiscreteTube([-6.4, 6.4], [-5.2, 5.2])
            sage: I = tube & P
            sage: I.an_element()
            (0, 0, 0) not in the plane
            trying similar points
            (0, 0, 1)

        """
        for o in self._objets:
            e = o.an_element()
            if e in self:
                return e
        else:
            raise ValueError("unable to find an element in this intersection")


class DiscreteBox(DiscreteSubset):
    r"""
    Cartesian product of intervals.

    INPUT:

    - ``*args`` - intervals, lists of size two : [min, max]

    EXAMPLES::

        sage: from slabbe import DiscreteBox
        sage: DiscreteBox([-5,5],[-5,5])
        Box: [-5, 5] x [-5, 5]

    ::

        sage: D = DiscreteBox([-3,3],[-3,3],[-3,3],[-3,3])
        sage: next(iter(D))
        (0, 0, 0, 0)

    TESTS::

        sage: d = DiscreteBox([-5,5], [-5,5], [-4,4])
        sage: d.edges_iterator().next()
        ((0, 0, 0), (1, 0, 0))

    """
    def __init__(self, *args):
        r"""
        Constructor.

        EXAMPLES::

            sage: from slabbe import DiscreteBox
            sage: b = DiscreteBox([2,10],[3,4])
            sage: b
            Box: [2, 10] x [3, 4]
            sage: b.roots()
            [(6, 4)]
        """
        self._intervals = args
        def predicate(p):
            return all(xmin <= x <= xmax for (x,(xmin,xmax)) in
                    itertools.izip(p,self._intervals))
        dim = len(self._intervals)
        root = tuple(round((xmax+xmin)/2) for (xmin,xmax) in self._intervals)
        DiscreteSubset.__init__(self, dimension=dim,
                predicate=predicate, roots=[root])

    def __str__(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from slabbe import DiscreteBox
            sage: str(DiscreteBox([2,10],[3,4]))
            '[2, 10] x [3, 4]'
            sage: str(DiscreteBox([-5,5], [-5,5], [-1,3]))
            '[-5, 5] x [-5, 5] x [-1, 3]'
        """
        L = ( "[%s, %s]" % (a, b) for (a,b) in self._intervals)
        return " x ".join(L)

    def _repr_(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from slabbe import DiscreteBox
            sage: DiscreteBox([-5,5])
            Box: [-5, 5]
            sage: DiscreteBox([-5,5], [-5,5])
            Box: [-5, 5] x [-5, 5]
            sage: DiscreteBox([-5,5], [-5,5], [-4,4])
            Box: [-5, 5] x [-5, 5] x [-4, 4]
        """
        return "Box: %s" % self

    def clip(self, space=1):
        r"""
        Return a good clip rectangle for this box.

        INPUT:

        - ``space`` -- number (default: ``1``), inner space within the box

        EXAMPLES::

            sage: from slabbe import DiscreteBox
            sage: box = DiscreteBox([-6,6],[-6,6])
            sage: box
            Box: [-6, 6] x [-6, 6]
            sage: box.clip()
            [(-5, -5), (5, -5), (5, 5), (-5, 5), (-5, -5)]

        ::

            sage: box = DiscreteBox([-6,6],[-4,3])
            sage: box.clip()
            [(-5, -3), (5, -3), (5, 2), (-5, 2), (-5, -3)]
        """
        if self.dimension() != 2:
            raise NotImplementedError("implemented only for dimension (=%s) 2" % self.dimension())
        xmin,xmax = self._intervals[0]
        ymin,ymax = self._intervals[1]
        xmin = xmin + space
        xmax = xmax - space
        ymin = ymin + space
        ymax = ymax - space
        return [(xmin,ymin), (xmax,ymin), (xmax, ymax), (xmin,ymax),
                (xmin,ymin)]

class DiscreteTube(DiscreteSubset):
    r"""
    Discrete Tube (preimage of a box by a projection matrix)

    Subset of a discrete object such that its projection by a matrix is
    inside a certain box.

    INPUT:

    - ``*args`` - intervals, lists of size two : [min, max]
    - ``projmat`` - matrix (default: ``M3to2``), projection matrix

    EXAMPLES::

        sage: from slabbe import DiscreteTube
        sage: DiscreteTube([-5,5],[-5,5])
        DiscreteTube: Preimage of [-5, 5] x [-5, 5] by a 2 by 3 matrix

    ::

        sage: m = matrix(3,4,range(12))
        sage: DiscreteTube([2,10],[3,4],[6,7], projmat=m)
        DiscreteTube: Preimage of [2, 10] x [3, 4] x [6, 7] by a 3 by 4 matrix

    EXAMPLES::

        sage: from slabbe import DiscretePlane, DiscreteTube
        sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
        sage: tube = DiscreteTube([-5,5],[-5,5])
        sage: I = p & tube
        sage: I
        Intersection of the following objects:
        Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8
        DiscreteTube: Preimage of [-5, 5] x [-5, 5] by a 2 by 3 matrix
        sage: len(list(I))
        115
    """
    @options(projmat=M3to2)
    def __init__(self, *args, **kwds):
        r"""
        Constructor.

        Return the subset of self such that the projection by the matrix is
        inside a certain box.

        EXAMPLES::

            sage: from slabbe import DiscreteTube
            sage: DiscreteTube([2,10],[3,4])
            DiscreteTube: Preimage of [2, 10] x [3, 4] by a 2 by 3 matrix

        """
        self._box = DiscreteBox(*args)
        projmat = kwds['projmat']
        def predicate(p):
            return projmat * p in self._box
        root = projmat.pseudoinverse() * self._box.an_element()
        root = vector(map(round, root))
        DiscreteSubset.__init__(self, dimension=projmat.ncols(),
                predicate=predicate,roots=[root])
        self._projmat = projmat

    def _repr_(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from slabbe import DiscreteTube
            sage: DiscreteTube([2,10],[3,4])
            DiscreteTube: Preimage of [2, 10] x [3, 4] by a 2 by 3 matrix
        """
        s = "DiscreteTube: Preimage of "
        s += str(self._box)
        s += " by a %s by %s matrix" % (self._projmat.nrows(), self._projmat.ncols())
        return s

    def clip(self, space=1):
        r"""
        Return a good clip rectangle for this box.

        INPUT:

        - ``space`` -- number (default: ``1``), inner space within the box

        EXAMPLES::

            sage: from slabbe import DiscreteTube
            sage: tube = DiscreteTube([-6,6],[-4,3])
            sage: tube.clip()
            [(-5, -3), (5, -3), (5, 2), (-5, 2), (-5, -3)]
        """
        return self._box.clip(space=space)

