# -*- coding: utf-8 -*-
r"""
Inflation rule

See Definition 5.17 in [BG13]_.

EXAMPLES:

This examples is taken from [BFG19]_::

    sage: z = polygen(QQ, 'z')
    sage: K = NumberField(z**2-z-1, 'tau', embedding=RR(1.6))
    sage: tau = K.gen()
    sage: from slabbe import InflationRule
    sage: import itertools
    sage: d = {(i,j):[] for i,j in itertools.product(range(4),repeat=2)}
    sage: d[(0,3)] = [vector(K, (tau,tau))]
    sage: d[(1,2)] = d[(1,3)] = [vector(K, (0,tau))]
    sage: d[(3,0)] = d[(3,1)] = d[(3,2)] = d[(3,3)] = [vector(K, (0,0))]
    sage: d[(2,1)] = d[(2,3)] = [vector(K, (tau,0))]
    sage: s = InflationRule(tau, 4, d)
    sage: s
    Inflation Rule (multiplier=tau) defined on 4 sets

TODO:

    - code Affine IFS, see for example:
      http://larryriddle.agnesscott.org/ifs/ifs.htm
      https://demonstrations.wolfram.com/AttractorsOfIteratedAffineTransformSystems/
      https://www.arxiv-vanity.com/papers/math/0604547/
      https://encyclopediaofmath.org/wiki/Iterated_function_system

    - or rename AffineIteratedTransformSystem to include both inflation and
      contractions ...

      code affine graph-directed ifs, see:

      http://jolivet.org/timo/docs/undfrac.pdf or recent [BV20]_

REFERENCES:

.. [BFG19] Michael Baake, Natalie Priebe Frank, Uwe Grimm. Three variations on a
   theme by Fibonacci. http://arxiv.org/abs/1910.00988

.. [BG13] Michael Baake, Uwe Grimm. Aperiodic order. Vol. 1. Vol. 149.
   Encyclopedia of Mathematics and its Applications. Cambridge University
   Press, Cambridge, 2013. http://www.ams.org/mathscinet-getitem?mr=3136260.

.. [BV20] Michael Barnsley, Andrew Vince. Tilings from Graph Directed
   Iterated Function Systems. Geometriae Dedicata, 9 août 2020.
   https://doi.org/10.1007/s10711-020-00560-4

"""
#*****************************************************************************
#       Copyright (C) 2020 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
import itertools

from sage.modules.free_module_element import vector

class InflationRule(object):
    r"""
    INPUT:

    - ``multiplier`` -- real number, inflation multiplier
    - ``number_of_sets`` -- integer
    - ``d`` -- dict, the displacement matrix, where each key (i,j) is
      mapped to a list of translations

    EXAMPLES:

    This examples is taken from [BFG19]_::

        sage: z = polygen(QQ, 'z')
        sage: K = NumberField(z**2-z-1, 'tau', embedding=RR(1.6))
        sage: tau = K.gen()
        sage: from slabbe import InflationRule
        sage: import itertools
        sage: d = {(i,j):[] for i,j in itertools.product(range(4),repeat=2)}
        sage: d[(0,3)] = [vector(K, (tau,tau))]
        sage: d[(1,2)] = d[(1,3)] = [vector(K, (0,tau))]
        sage: d[(3,0)] = d[(3,1)] = d[(3,2)] = d[(3,3)] = [vector(K, (0,0))]
        sage: d[(2,1)] = d[(2,3)] = [vector(K, (tau,0))]
        sage: s = InflationRule(tau, 4, d)
        sage: s
        Inflation Rule (multiplier=tau) defined on 4 sets

    """
    def __init__(self, multiplier, number_of_sets, d):
        r"""
        See class documentation.

        EXAMPLES::

            sage: from slabbe import inflation_rules
            sage: s = inflation_rules.direct_product_two_Fibonacci_chains()
        """
        self._multiplier = multiplier
        self._number_of_sets = number_of_sets
        self._d = d

    @classmethod
    def from_one_dimensional_substitution(cls, m):
        r"""
        Return the inflation rule defined by a unidimensional primitive
        substitution

        INPUT:

        - ``m`` -- WordMorphism, primitive substitution

        EXAMPLES::

            sage: from slabbe import InflationRule
            sage: m = WordMorphism('a->ab,b->a')
            sage: m
            WordMorphism: a->ab, b->a
            sage: s = InflationRule.from_one_dimensional_substitution(m)
            sage: s
            Inflation Rule (multiplier=root) defined on 2 sets

        """
        from slabbe.matrices import perron_right_eigenvector_in_number_field
        M = m.incidence_matrix()
        root, perron_left = perron_right_eigenvector_in_number_field(M.T)
        K = root.parent()
        alphabet = m.domain().alphabet()
        number_of_sets = alphabet.cardinality()
        d = {(i,j):[] for i,j in itertools.product(range(number_of_sets),repeat=2)}
        for i,a in enumerate(alphabet):
            m_a = m(a)
            pos = vector(K,(0,0))
            for b in m_a:
                j = alphabet.index(b)
                d[(i,j)].append(pos)
                pos += vector(K,(perron_left[j],0))  # horizontal vector
        return InflationRule(root, number_of_sets, d)

    @classmethod
    def from_two_dimensional_substitution(cls, m):
        r"""
        Return the inflation rule defined by a 2-dimensional primitive
        substitution

        INPUT:

        - ``s`` -- Substitution2d, primitive substitution

        EXAMPLES::

            sage: from slabbe import InflationRule, Substitution2d
            sage: d = {}                         # not tested
            sage: m = Substitution2d(d)          # not tested
            sage: s = InflationRule.from_one_dimensional_substitution(s)  # not tested
            sage: s                              # not tested
            Inflation Rule (multiplier=root) defined on 2 sets

        """
        raise NotImplementedError

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import inflation_rules
            sage: inflation_rules.direct_product_two_Fibonacci_chains()
            Inflation Rule (multiplier=tau) defined on 4 sets
        """
        return ("Inflation Rule (multiplier={})"
                " defined on {} sets".format(self._multiplier, 
                                             self._number_of_sets))

    def __eq__(self, other):
        r"""
        INPUT:

        - ``other`` -- inflation rule

        EXAMPLES::

            sage: from slabbe import inflation_rules
            sage: s = inflation_rules.direct_product_two_Fibonacci_chains()
            sage: t = inflation_rules.direct_product_two_Fibonacci_chains()
            sage: s == t
            True

        """
        return (isinstance(other, InflationRule) and 
                self._multiplier == other._multiplier and
                self._d == other._d)

    def number_of_sets(self):
        return self._number_of_sets

    def ambient_dimension(self):
        r"""
        Return the ambient dimension of the space

        EXAMPLES::

            sage: from slabbe import inflation_rules
            sage: s = inflation_rules.direct_product_two_Fibonacci_chains()
            sage: s.ambient_dimension()
            2
        """
        for v in self._d.values():
            if v:
                return len(v[0])
        else:
            return None

    def galois_conjugate(self):
        r"""
        EXAMPLES:

        Fibonacci substitution::

            sage: from slabbe import InflationRule
            sage: m = WordMorphism('a->ab,b->a')
            sage: s = InflationRule.from_one_dimensional_substitution(m)
            sage: s.galois_conjugate()
            Inflation Rule (multiplier=-root + 1) defined on 2 sets

        Direct Product of 2 Fibonacci::

            sage: from slabbe import inflation_rules
            sage: s = inflation_rules.direct_product_two_Fibonacci_chains()
            sage: s.galois_conjugate()
            Inflation Rule (multiplier=-tau + 1) defined on 4 sets

        """
        d = {k:[vector(a.galois_conjugate() for a in t) for t in v] 
             for (k,v) in self._d.items()}
        return InflationRule(self._multiplier.galois_conjugate(),
                             self._number_of_sets, d)

    def __call__(self, S):
        r"""
        Return the image of the list of list of points.

        INPUT:

        - ``S`` -- list of list of points

        EXAMPLES::
            
            sage: from slabbe import inflation_rules
            sage: s = inflation_rules.direct_product_two_Fibonacci_chains()
            sage: s([[],[],[],[]])
            [[], [], [], []]

        ::

            sage: zero = vector((0,0))
            sage: s([[],[],[],[zero]])
            [[(tau, tau)], [(0, tau)], [(tau, 0)], [(0, 0)]]
            sage: s(_)
            [[(tau, tau)],
             [(tau + 1, tau), (0, tau)],
             [(tau, tau + 1), (tau, 0)],
             [(tau + 1, tau + 1), (0, tau + 1), (tau + 1, 0), (0, 0)]]
            sage: s(_)
            [[(3*tau + 1, 3*tau + 1), (tau, 3*tau + 1), (3*tau + 1, tau), (tau, tau)],
             [(tau + 1, 3*tau + 1),
              (tau + 1, tau),
              (2*tau + 1, 3*tau + 1),
              (0, 3*tau + 1),
              (2*tau + 1, tau),
              (0, tau)],
             [(3*tau + 1, tau + 1),
              (tau, tau + 1),
              (3*tau + 1, 2*tau + 1),
              (tau, 2*tau + 1),
              (3*tau + 1, 0),
              (tau, 0)],
             [(tau + 1, tau + 1),
              (2*tau + 1, tau + 1),
              (0, tau + 1),
              (tau + 1, 2*tau + 1),
              (tau + 1, 0),
              (2*tau + 1, 2*tau + 1),
              (0, 2*tau + 1),
              (2*tau + 1, 0),
              (0, 0)]]

        TESTS::

            sage: s([[],[],[]])
            Traceback (most recent call last):
            ...
            ValueError: size of input (=3) must match the number of sets of
            this inflation rule (=4)

        """
        if not len(S) == self._number_of_sets:
            raise ValueError("size of input (={}) must match the number of"
                    " sets of this inflation rule (={})".format(len(S),
                        self._number_of_sets))
        S_image = []
        for i in range(self._number_of_sets):
            S_image_i = []
            for j in range(self._number_of_sets):
                expanded = [self._multiplier * v for v in S[j]]
                translations = self._d[(i,j)]
                S_image_i.extend(v + t for v in expanded 
                                       for t in translations)
            S_image.append(S_image_i)
        return S_image

    def ifs(self, n_iterations=3, S=None):
        r"""
        Return the IFS after few iterations

        INPUT:

        - ``n_iterations`` -- integer (default: ``3``)
        - ``S`` -- list of list of points (default: ``None``)

        OUTPUT:

        list of list of points

        EXAMPLES::

            sage: from slabbe import inflation_rules
            sage: s = inflation_rules.direct_product_two_Fibonacci_chains()
            sage: s.ifs(1)
            [[], [], [], [(0, 0)]]
            sage: s.ifs(0)
            [[(0, 0)], [], [], []]
            sage: s.ifs(1)
            [[], [], [], [(0, 0)]]
            sage: s.ifs(2)
            [[(tau, tau)], [(0, tau)], [(tau, 0)], [(0, 0)]]
            sage: s.ifs(3)
            [[(tau, tau)],
             [(tau + 1, tau), (0, tau)],
             [(tau, tau + 1), (tau, 0)],
             [(tau + 1, tau + 1), (0, tau + 1), (tau + 1, 0), (0, 0)]]
        
        """
        if S is None:
            S = [[] for _ in range(self._number_of_sets)]
            zero = vector([0]*self.ambient_dimension())
            S[0].append(zero)
        for _ in range(n_iterations):
            S = self(S)
        return S

    def plot_ifs(self, n_iterations=3):
        r"""
        Return a graphic image of the IFS after few iterations

        INPUT:

        - ``n_iterations`` -- integer (default: ``3``)
        - ``S`` -- list of list of points (default: ``None``)

        OUTPUT:

        Graphics object

        EXAMPLES:
        
        The usual Fibonacci chain::

            sage: from slabbe import InflationRule
            sage: m = WordMorphism('a->ab,b->a')
            sage: s = InflationRule.from_one_dimensional_substitution(m)
            sage: G = s.plot_ifs(7)

        Its contracting IFS::

            sage: sc = s.galois_conjugate()
            sage: G = sc.plot_ifs(10)

        The direct product of two Fibonacci chains::

            sage: from slabbe import inflation_rules
            sage: s = inflation_rules.direct_product_two_Fibonacci_chains()
            sage: G = s.plot_ifs(7)

        This inflation rule is related to a contracting IFS whose unique
        solution is given in formula (4.5) of [BFG19]_::

            sage: sc = s.galois_conjugate()
            sage: G = sc.plot_ifs(10)
        """
        from sage.plot.colors import rainbow
        from sage.plot.point import points
        from sage.plot.graphics import Graphics
        G = Graphics()
        colors = rainbow(self._number_of_sets)
        ifs = self.ifs(n_iterations=n_iterations)
        for i,P in enumerate(ifs):
            G += points(P, color=colors[i], legend_label=str(i))
        return G

class inflation_rules_generator(object):
    def direct_product_two_Fibonacci_chains(self):
        r"""
        EXAMPLES::

            sage: from slabbe import inflation_rules
            sage: inflation_rules.direct_product_two_Fibonacci_chains()
            Inflation Rule (multiplier=tau) defined on 4 sets

        """
        from sage.rings.polynomial.polynomial_ring import polygen
        from sage.rings.real_mpfr import RR
        from sage.rings.rational_field import QQ
        from sage.rings.number_field.number_field import NumberField
        z = polygen(QQ, 'z')
        K = NumberField(z**2-z-1, 'tau', embedding=RR(1.6))
        tau = K.gen()
        d = {(i,j):[] for i,j in itertools.product(range(4),repeat=2)}
        d[(0,3)] = [vector(K, (tau,tau))]
        d[(1,2)] = d[(1,3)] = [vector(K, (0,tau))]
        d[(3,0)] = d[(3,1)] = d[(3,2)] = d[(3,3)] = [vector(K, (0,0))]
        d[(2,1)] = d[(2,3)] = [vector(K, (tau,0))]
        return InflationRule(tau, 4, d)

inflation_rules = inflation_rules_generator()

