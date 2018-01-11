# -*- coding: utf-8 -*-
r"""
Billiard words

EXAMPLES::

    sage: from slabbe import BilliardCube
    sage: b = BilliardCube((1,pi,sqrt(2)))
    sage: b
    Cubic billiard of direction (1, pi, sqrt(2))

TODO:

    - Rewrite some parts in cython because it is slow
    - Should handle any direction
    - Should use Forest structure for enumeration
    - Should use +e_i only for children
    - Fix documentation of class
    - Fix issue with the assertion error in the step iterator
    - not robust for non integral start point

"""

#*****************************************************************************
#       Copyright (C) 2010-2014 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
from copy import copy
from sage.modules.free_module_element import vector
from sage.combinat.words.word import Word
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
from slabbe.discrete_plane import DiscretePlane 
from slabbe.discrete_subset import DiscreteSubset, Intersection
################################################
# Discrete Line
################################################
class BilliardCube(Intersection):
    r"""
    This is the set of point `p` such that

        `0 \leq  p \cdot v - mu < \omega`  #fix me
        `0 \leq  p \cdot v - mu < \omega`  #fix me
        `0 \leq  p \cdot v - mu < \omega`  #fix me

    INPUT:

    - ``v`` - directive vector
    - ``start`` - initial point (default = (0,0,0))

    EXAMPLES::

        sage: from slabbe import BilliardCube
        sage: b = BilliardCube((1,pi,sqrt(2)))
        sage: b
        Cubic billiard of direction (1, pi, sqrt(2))

    ::

        sage: b = BilliardCube((1,pi,sqrt(2)))
        sage: it = iter(b)
        sage: [next(it) for _ in range(20)]
        [(0, 0, 0),
         (0, 1, 0),
         (0, 1, 1),
         (0, 2, 1),
         (1, 2, 1),
         (1, 3, 1),
         (1, 3, 2),
         (1, 4, 2),
         (1, 5, 2),
         (2, 5, 2),
         (2, 6, 2),
         (2, 6, 3),
         (2, 7, 3),
         (2, 8, 3),
         (2, 8, 4),
         (3, 8, 4),
         (3, 9, 4),
         (3, 10, 4),
         (3, 10, 5),
         (3, 11, 5)]

     ::

         sage: b = BilliardCube((1,sqrt(2),pi), start=(11,13,14))
         sage: b.to_word()
         word: 3231323313233213323132331233321332313233...

    """
    def __init__(self, v, start=(0,0,0)):
        r"""
        EXAMPLES::

            sage: from slabbe import BilliardCube
            sage: b = BilliardCube((1,pi,sqrt(2)))
            sage: b
            Cubic billiard of direction (1, pi, sqrt(2))

        TESTS::

            sage: vector((0,0,0)) in b
            True
            sage: vector((0,0,1)) in b
            False
            sage: vector((0,1,0)) in b
            True
            sage: vector((1,0,0)) in b
            False
            sage: vector((0,-1,0)) in b
            True

        """
        a,b,c = self._v = vector(v)
        sx,sy,sz = self._start = vector(start)
        px = DiscretePlane([0,c,-b], b+c, mu=(b+c)/2 - sy*c + sz*b)
        py = DiscretePlane([c,0,-a], a+c, mu=(a+c)/2 - sx*c + sz*a)
        pz = DiscretePlane([b,-a,0], a+b, mu=(a+b)/2 - sx*b + sy*a)
        Intersection.__init__(self, (px,py,pz))

    def _repr_(self):
        r"""
        EXAMPLES::
        
            sage: from slabbe import BilliardCube
            sage: b = BilliardCube((1,pi,sqrt(2)))
            sage: b
            Cubic billiard of direction (1, pi, sqrt(2))
        """
        return "Cubic billiard of direction %s" % self._v

    def an_element(self):
        r"""
        Returns an element in self.

        EXAMPLES::

            sage: from slabbe import BilliardCube
            sage: b = BilliardCube((1,pi,sqrt(2)))
            sage: b.an_element()
            (0, 0, 0)
        """
        v = copy(self._start)
        v.set_immutable()
        return v

    def children(self, p):
        r"""
        Return the children of a point.

        This method overwrites the methods
        :meth:`slabbe.discrete_subset.DiscreteSubset.children`, because for
        billiard words, we go only in one direction in each axis.

        EXAMPLES::

            sage: from slabbe import BilliardCube
            sage: b = BilliardCube((1,pi,sqrt(2)))
            sage: list(b.children(vector((0,0,0))))
            [(0, 1, 0)]
        """
        for v in self.base_edges():
            u = p+v
            u.set_immutable()
            if u in self: yield u

    def connected_component_iterator(self, roots=None):
        r"""
        Return an iterator over the connected component of the root.

        This method overwrites the methods
        :meth:`slabbe.discrete_subset.DiscreteSubset.connected_component_iterator`,
        because for billiard words, we go only in one direction in each
        axis which allows to use a forest structure for the enumeration.

        INPUT:

        - ``roots`` - list of some elements immutable in self

        EXAMPLES::

            sage: from slabbe import BilliardCube
            sage: p = BilliardCube([1,pi,sqrt(7)])
            sage: root = vector((0,0,0))
            sage: root.set_immutable()
            sage: it = p.connected_component_iterator(roots=[root])
            sage: [next(it) for _ in range(5)]
            [(0, 0, 0), (0, 1, 0), (0, 1, 1), (0, 2, 1), (1, 2, 1)]

        ::

            sage: p = BilliardCube([1,pi,7.45], start=(10.2,20.4,30.8))
            sage: it = p.connected_component_iterator()
            sage: [next(it) for _ in range(5)]
            [(10.2000000000000, 20.4000000000000, 30.8000000000000),
             (10.2000000000000, 20.4000000000000, 31.8000000000000),
             (10.2000000000000, 21.4000000000000, 31.8000000000000),
             (10.2000000000000, 21.4000000000000, 32.8000000000000),
             (10.2000000000000, 21.4000000000000, 33.8000000000000)]
        """
        roots = roots if roots else [self.an_element()]
        if not all(root in self for root in roots):
            raise ValueError("roots (=%s) must all be in self(=%s)" % (roots, self))
        #for root in roots:
        #    root.set_immutable()
        C = RecursivelyEnumeratedSet(seeds=roots, successors=self.children, structure='forest')
        return C.breadth_first_search_iterator()

    def step_iterator(self):
        r"""
        Return an iterator coding the steps of the discrete line.

        EXAMPLES::

            sage: from slabbe import BilliardCube
            sage: b = BilliardCube((1,pi,sqrt(2)))
            sage: it = b.step_iterator()
            sage: [next(it) for _ in range(5)]
            [(0, 1, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 1, 0)]

        TESTS:

        Fix this::

            sage: from slabbe import BilliardCube
            sage: B = BilliardCube((1.1,2.2,3.3))
            sage: B.to_word()
            Traceback (most recent call last):
            ...
            AssertionError: step(=(-1, 0, 1)) is not a canonical basis
            vector.
        """
        possible_steps = map(vector, ((1,0,0), (0,1,0), (0,0,1)))
        i = iter(self)
        j = iter(self)
        j.next()
        while True:
            step = j.next() - i.next()
            step.set_immutable()
            assert step in possible_steps, ("step(=%s) is not a " % step +
                      "canonical basis vector.")
            yield step

    def to_word(self, alphabet=[1,2,3]):
        r"""
        Return the billiard word.

        INPUT:

        - ``alphabet`` - list

        EXAMPLES::

            sage: from slabbe import BilliardCube
            sage: b = BilliardCube((1,pi,sqrt(2)))
            sage: b.to_word()
            word: 2321232212322312232123221322231223212322...

        ::

            sage: B = BilliardCube((sqrt(3),sqrt(5),sqrt(7)))
            sage: B.to_word()
            word: 3213213231232133213213231232132313231232...

        """
        steps = map(vector, ((1,0,0), (0,1,0), (0,0,1)))
        for step in steps: step.set_immutable()
        coding = dict(zip(steps, alphabet)) 
        it = (coding[step] for step in self.step_iterator())
        #WP = WordPaths(alphabet, [(1,0,0),(0,1,0),(0,0,1)])
        return Word(it, alphabet=alphabet)

