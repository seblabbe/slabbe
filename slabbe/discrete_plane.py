# -*- coding: utf-8 -*-
r"""
Discrete Hyperplanes

Intersection of a plane and a tube::

    sage: from slabbe import DiscretePlane, DiscreteTube
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
    sage: L = DiscreteLine([pi,sqrt(2)], pi+sqrt(2), mu=0)
    sage: b = DiscreteBox([-5,5],[-5,5])
    sage: I = L & b
    sage: I
    Intersection of the following objects:
    Set of points x in ZZ^2 satisfying: 0 <= (pi, sqrt(2)) . x + 0 < pi + sqrt(2)
    [-5, 5] x [-5, 5]


TODO:

    - do some dimension checking for DiscreteLine and DiscretePlane

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
from sage.rings.real_mpfr import RealField
from sage.functions.other import ceil, sqrt
from sage.modules.free_module_element import vector
from sage.misc.cachefunc import cached_method
from slabbe.discrete_subset import DiscreteSubset
################################################
# Discrete plane and hyperplanes
################################################
class DiscreteHyperplane(DiscreteSubset):
    r"""
    This is the set of point `p` such that

        `0 \leq  p \cdot v - mu < \omega`

    INPUT:

    - ``v`` - normal vector
    - ``omega`` - width
    - ``mu`` - intercept (optional, default: 0)

    EXAMPLES::

        sage: from slabbe import DiscreteLine
        sage: L = DiscreteLine([pi,sqrt(2)], pi+sqrt(2), mu=10)
        sage: L
        Set of points x in ZZ^2 satisfying: 0 <= (pi, sqrt(2)) . x + 10 < pi + sqrt(2)

    ::

        sage: from slabbe import DiscretePlane
        sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
        sage: p
        Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8

    ::

        sage: from slabbe import DiscreteHyperplane
        sage: p = DiscreteHyperplane([1,3,7,9], 20, mu=13)
        sage: p
        Set of points x in ZZ^4 satisfying: 0 <= (1, 3, 7, 9) . x + 13 < 20

    TESTS::

        sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=20)
        sage: vector((0,0,0)) in p
        False
        sage: p = DiscretePlane([1,pi,7], 1+pi+7, mu=0)
        sage: vector((0,0,0)) in p
        True

    ::

        sage: p = DiscreteHyperplane((2,3,4,5), 10)
        sage: p.dimension()
        4

    ::

        sage: L = DiscreteLine([1,pi], 1+pi, mu=20)
        sage: vector((0,0)) in L
        False
        sage: L = DiscreteLine([1,pi], 1+pi, mu=0)
        sage: vector((0,0)) in L
        True
    """
    def __init__(self, v, omega, mu=0, prec=None):
        r"""
        EXAMPLES::

            sage: from slabbe import DiscreteHyperplane
            sage: p = DiscreteHyperplane([1,pi,7], 1+pi+7, mu=0)
            sage: p
            Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 0 < pi + 8

        ::

            sage: p = DiscreteHyperplane([1,pi,7], 1+pi+7, mu=20)
            sage: vector((0,0,0)) in p
            False
            sage: p = DiscreteHyperplane([1,pi,7], 1+pi+7, mu=0)
            sage: vector((0,0,0)) in p
            True
        """
        if prec is None:
            self._v = vector(v)
            self._omega = omega
            self._mu = mu
        else:
            RF = RealField(prec=prec)
            self._v = vector(RF, v)
            self._omega = RF(omega)
            self._mu = RF(mu)
        def contain(p):
            #print("est-ce proche : ", self._v.dot_product(p) + self._mu)
            return  0 <= self._v.dot_product(p) + self._mu < self._omega
        DiscreteSubset.__init__(self, dimension=len(self._v), predicate=contain)

    @cached_method
    def roots(self):
        r"""
        Return the roots, i.e., a list of elements in self.

        It also makes sure the roots are in self and raises an error otherwise.

        EXAMPLES::

            sage: from slabbe import DiscretePlane
            sage: P = DiscretePlane([3,4,5], 12, mu=20)
            sage: P.roots()
            [(-1, -1, -1)]
            sage: all(p in P for p in P.roots())
            True
        """
        p = self.an_element()
        p = self._space(p)
        p.set_immutable()
        if not p in self: 
            raise ValueError("root element (={}) provided at"
                    " initialisation is not in self".format(p))
        self._roots = [p]
        return self._roots

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from slabbe import DiscreteHyperplane
            sage: p = DiscreteHyperplane([1,pi,7], 1+pi+7, mu=10)
            sage: p
            Set of points x in ZZ^3 satisfying: 0 <= (1, pi, 7) . x + 10 < pi + 8
        """
        s = "Set of points x in ZZ^{} satisfying: ".format(self.dimension())
        s +=  "0 <= {} . x + {} < {}".format(self._v,self._mu, self._omega)
        return s

    def an_element(self, x=0, y=0):
        r"""
        Returns an element in self.

        EXAMPLES::

            sage: from slabbe import DiscreteHyperplane
            sage: p = DiscreteHyperplane([1,pi,7], 1+pi+7, mu=10)
            sage: p.an_element()
            (0, 0, 0)

        ::

            sage: from slabbe import DiscreteLine
            sage: L = DiscreteLine([pi,sqrt(2)], pi+sqrt(2), mu=10)
            sage: L.an_element()
            (-2, -2)

        ::

            sage: L = DiscreteLine([pi,sqrt(2)], pi+sqrt(2), mu=0)
            sage: L.an_element()
            (0, 0)

        """
        dimension = len(self._v)
        if dimension == 2:
            return self._an_element_2d(x=x,y=y)
        elif dimension == 3:
            return self._an_element_3d(x=x,y=y)
        else:
            raise NotImplementedError("implemented only for dimension(=%s) 2 or 3" % dimension)

    def _an_element_3d(self, x=0, y=0):
        r"""
        Returns an element in self.

        EXAMPLES::

            sage: from slabbe import DiscreteHyperplane
            sage: p = DiscreteHyperplane([1,pi,7], 1+pi+7, mu=10)
            sage: p._an_element_3d()
            (0, 0, 0)
        """
        a,b,c = self._v
        x_sqrt3 = ceil(x / sqrt(3))
        left  = ((a+b) * y + (a-b) * x_sqrt3 - self._mu) / (a+b+c)
        right = ((a+b) * y + (a-b) * x_sqrt3 - self._mu + self._omega) / (a+b+c)
        #print("left, right = ", left, right)
        #print("left, right = ", ceil(left), ceil(right)-1)
        # left <= z <= right
        znew = ceil(left)
        xnew = znew - y - x_sqrt3
        ynew = znew - y + x_sqrt3
        znew = ceil(right)-1
        #print(xnew, ynew, znew)
        #print(vector((xnew, ynew, znew)) in self)
        #print(vector((x,y,ceil(right)-1)) in self)
        v = vector((xnew, ynew, znew))
        if v in self:
            v.set_immutable()
            return v
        else:
            print("%s not in the plane" % v)
            print("trying similar points")
            v = vector((xnew, ynew, znew-1))
            if v in self:
                v.set_immutable()
                return v
            v = vector((xnew, ynew, znew+1))
            if v in self:
                v.set_immutable()
                return v
        raise ValueError("%s not in the plane" % v)

        # minimum = - floor(self._mu)
        # maximum = ceil(self._omega - self._mu) - 1
        # # print(minimum, maximum)
        # # minimum <= z <= maximum
        # assert vector((0,0,maximum)) in self
        # return vector((0,0,maximum))
    def _an_element_2d(self, x=0, y=0):
        r"""
        Returns an element in self.

        EXAMPLES::

            sage: from slabbe import DiscreteLine
            sage: L = DiscreteLine([pi,sqrt(2)], pi+sqrt(2), mu=10)
            sage: L._an_element_2d()
            (-2, -2)

        """
        assert self.dimension() == 2, "dimension must be 2"
        p = vector((x,y))
        while True:
            if p in self:
                p.set_immutable()
                return p
            elif -p in self:
                mp = -p
                mp.set_immutable()
                return mp
            else:
                p += vector((1,1))


    def level_value(self, p):
        r"""
        Return the level value of a point p.

        INPUT:

        - ``p`` - point in the space

        EXAMPLES::

            sage: from slabbe import DiscreteHyperplane
            sage: H = DiscreteHyperplane([1,3,7,9], 20, mu=13)
            sage: p = H._space((1,2,3,4))
            sage: H.level_value(p)
            64
        """
        return self._v.dot_product(p)

DiscretePlane = DiscreteHyperplane
DiscreteLine = DiscreteHyperplane
