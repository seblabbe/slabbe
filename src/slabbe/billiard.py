# -*- coding: utf-8 -*-
r"""
Billiard words

EXAMPLES:

    ...

TODO:

    

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
from copy import copy
from sage.modules.free_module_element import vector
from sage.combinat.words.word import Word
from slabbe.discrete_plane import Plan 
from slabbe.discrete_subset import DiscreteSubset, Intersection
################################################
# Discrete Line
################################################
class BillardCube(Intersection):
    r"""
    This is the set of point `p` such that

        `0 \leq  p \cdot v - mu < \omega`  #fix me
        `0 \leq  p \cdot v - mu < \omega`  #fix me
        `0 \leq  p \cdot v - mu < \omega`  #fix me

    INPUT:

    - ``v`` - vecteur directeur
    - ``start`` - point initial (default = (0,0,0))

    EXAMPLES::

        sage: b = BillardCube((1,pi,sqrt(2)))
        sage: b
        Billard cubique de direction (1, pi, sqrt(2))

    ::

        sage: b = BillardCube((1,pi,sqrt(2)))
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

    """
    def __init__(self, v, start=(0,0,0)):
        r"""
        EXAMPLES::

            sage: b = BillardCube((1,pi,sqrt(2)))
            sage: b
            Billard cubique de direction (1, pi, sqrt(2))

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
            False

        """
        a,b,c = v
        self._v = vector(v)
        self._start = vector(start)
        assert a>=0 and b>=0 and c>=0, "We assume positive entries for now"
        px = Plan([0,c,-b], b+c, mu=(b+c)/2)
        py = Plan([c,0,-a], a+c, mu=(a+c)/2)
        pz = Plan([b,-a,0], a+b, mu=(a+b)/2)
        def is_positive(p):
            return all(p[i]>=self._start[i] for i in range(len(p)))
        pos = DiscreteSubset(dimension=3, predicate=is_positive)
        Intersection.__init__(self, (px,py,pz,pos))

    def _repr_(self):
        r"""
        EXAMPLES::
        
            sage: b = BillardCube((1,pi,sqrt(2)))
            sage: b
            Billard cubique de direction (1, pi, sqrt(2))
        """
        return "Billard cubique de direction %s" % self._v

    def an_element(self):
        r"""
        Returns an element in self.

        EXAMPLES::

            sage: b = BillardCube((1,pi,sqrt(2)))
            sage: b.an_element()
            (0, 0, 0)
        """
        v = copy(self._start)
        v.set_immutable()
        return v

    def step_iterator(self):
        r"""
        Return an iterator coding the steps of the discrete line.

        EXAMPLES::

            sage: b = BillardCube((1,pi,sqrt(2)))
            sage: it = b.step_iterator()
            sage: [next(it) for _ in range(5)]
            [(0, 1, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 1, 0)]
        """
        i = iter(self)
        j = iter(self)
        j.next()
        while True:
            step = j.next() - i.next()
            step.set_immutable()
            yield step

    def to_word(self, alphabet=[1,2,3]):
        r"""
        Return the billiard word.

        INPUT:

        - ``alphabet`` - list

        EXAMPLES::

            sage: b = BillardCube((1,pi,sqrt(2)))
            sage: b.to_word()
            word: 2321232212322312232123221322231223212322...

        """
        steps = map(vector, ((1,0,0), (0,1,0), (0,0,1)))
        for step in steps: step.set_immutable()
        coding = dict(zip(steps, alphabet)) 
        it = (coding[step] for step in self.step_iterator())
        #WP = WordPaths(alphabet, [(1,0,0),(0,1,0),(0,0,1)])
        return Word(it, alphabet=alphabet)

