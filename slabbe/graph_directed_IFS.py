# -*- coding: utf-8 -*-
r"""
Graph-directed iterated function system (GIFS)

See [JK14] or [BV20] or

- http://larryriddle.agnesscott.org/ifs/ifs.htm
- https://encyclopediaofmath.org/wiki/Iterated_function_system

We allow the functions to be contracting or not. When the functions are
inflations, it allows to represent inflation rules as in Definition 5.17 of
[BG13]_.

EXAMPLES:

The Cantor set::

    sage: F = AffineGroup(1, QQ)
    sage: f1 = F.linear(1/3)
    sage: f2 = F(1/3, vector([2/3]))
    sage: f1
    x |-> [1/3] x + [0]
    sage: f2
    x |-> [1/3] x + [2/3]
    sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
    sage: cantor_IFS = GIFS(QQ^1, [(0,0,f1),(0,0,f2)])
    sage: cantor_IFS
    GIFS defined by 2 maps on 
    Vector space of dimension 1 over Rational Field

REFERENCES:

.. [JK14] Jolivet, Timo, et Jarkko Kari. « Undecidable Properties of Self-Affine
   Sets and Multi-Tape Automata ». In Mathematical Foundations of Computer
   Science 2014, édité par Erzsébet Csuhaj-Varjú, Martin Dietzfelbinger,
   et Zoltán Ésik, 8634:352‑64. Berlin, Heidelberg: Springer Berlin
   Heidelberg, 2014. https://doi.org/10.1007/978-3-662-44522-8_30.

.. [BV20] Michael Barnsley, Andrew Vince. Tilings from Graph Directed
   Iterated Function Systems. Geometriae Dedicata, 9 août 2020.
   https://doi.org/10.1007/s10711-020-00560-4

.. [BG13] Michael Baake, Uwe Grimm. Aperiodic order. Vol. 1. Vol. 149.
   Encyclopedia of Mathematics and its Applications. Cambridge University
   Press, Cambridge, 2013. http://www.ams.org/mathscinet-getitem?mr=3136260.

.. [BFG19] Michael Baake, Natalie Priebe Frank, Uwe Grimm. Three variations on a
   theme by Fibonacci. http://arxiv.org/abs/1910.00988

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

class GraphDirectedIteratedFunctionSystem(object):
    r"""
    INPUT:

    - ``module`` -- the module on which the functions are defined
    - ``edges`` -- list, list of triples (u,v,f) where f is a function
      associated to the directed edge (u,v).

    EXAMPLES:

    The Cantor set::

        sage: F = AffineGroup(1, QQ)
        sage: f1 = F.linear(1/3)
        sage: f2 = F(1/3, vector([2/3]))
        sage: f1
        x |-> [1/3] x + [0]
        sage: f2
        x |-> [1/3] x + [2/3]
        sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
        sage: GIFS(QQ^1, [(0,0,f1),(0,0,f2)])
        GIFS defined by 2 maps on 
        Vector space of dimension 1 over Rational Field

    """
    def __init__(self, module, edges):
        r"""
        See class documentation.

        EXAMPLES::

            sage: F = AffineGroup(1, QQ)
            sage: f1 = F.linear(1/3)
            sage: f2 = F(1/3, vector([2/3]))
            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: ifs = GIFS(QQ^1, [(0,0,f1),(0,0,f2)])
        """
        self._module = module
        self._edges = edges

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: F = AffineGroup(1, QQ)
            sage: f1 = F.linear(1/3)
            sage: f2 = F(1/3, vector([2/3]))
            sage: GIFS(QQ^1, [(0,0,f1),(0,0,f2)])
            GIFS defined by 2 maps on 
            Vector space of dimension 1 over Rational Field

        """
        return ("GIFS defined by {} maps on {}".format(len(self._edges),
                self._module))

    @classmethod
    def from_one_dimensional_substitution(cls, m):
        r"""
        Return the GIFS defined by a unidimensional primitive
        substitution

        INPUT:

        - ``m`` -- WordMorphism, primitive substitution

        EXAMPLES::

            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: m = WordMorphism('a->ab,b->a')
            sage: g = GIFS.from_one_dimensional_substitution(m)
            sage: g
            GIFS defined by 3 maps on
            Vector space of dimension 1 over
            Number Field in root with defining polynomial x^2 - x - 1 with
            root = 1.618033988749895?

        """
        from slabbe.matrices import perron_right_eigenvector_in_number_field
        M = m.incidence_matrix()
        root, perron_left = perron_right_eigenvector_in_number_field(M.T)
        K = root.parent()
        alphabet = m.domain().alphabet()
        size = alphabet.cardinality()
        module = K**1
        d = {(i,j):[] for i,j in itertools.product(range(size),repeat=2)}
        for i,a in enumerate(alphabet):
            m_a = m(a)
            pos = module.zero()
            for b in m_a:
                j = alphabet.index(b)
                d[(i,j)].append(pos)
                pos += module([perron_left[j]])
        return cls.from_inflation_rule(module, root, d)

    @classmethod
    def from_two_dimensional_substitution(cls, m):
        r"""
        Return the GIFS defined by a 2-dimensional primitive
        substitution

        INPUT:

        - ``s`` -- Substitution2d, primitive substitution

        EXAMPLES::

            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: from slabbe import Substitution2d
            sage: d = {}                         # not tested
            sage: m = Substitution2d(d)          # not tested
            sage: s = InflationRule.from_one_dimensional_substitution(s)  # not tested
            sage: s                              # not tested
            Inflation Rule (multiplier=root) defined on 2 sets

        """
        raise NotImplementedError

    @classmethod
    def from_inflation_rule(cls, module, multiplier, displacement_matrix):
        r"""
        Return the GIFS defined by a 2-dimensional primitive
        substitution

        INPUT:

        - ``module`` -- module or vector space
        - ``multiplier`` -- real number, inflation multiplier
        - ``d`` -- dict, the displacement matrix, where each key (i,j) is
          mapped to a list of translations

        EXAMPLES:

        This examples is taken from [BFG19]_::

            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: z = polygen(QQ, 'z')
            sage: K = NumberField(z**2-z-1, 'tau', embedding=RR(1.6))
            sage: tau = K.gen()
            sage: import itertools
            sage: d = {(i,j):[] for i,j in itertools.product(range(4),repeat=2)}
            sage: d[(0,3)] = [vector(K, (tau,tau))]
            sage: d[(1,2)] = d[(1,3)] = [vector(K, (0,tau))]
            sage: d[(2,1)] = d[(2,3)] = [vector(K, (tau,0))]
            sage: d[(3,0)] = d[(3,1)] = d[(3,2)] = d[(3,3)] = [vector(K, (0,0))]
            sage: GIFS.from_inflation_rule(K^2, tau, d)
            GIFS defined by 9 maps on Vector space of dimension 2 over
            Number Field in tau with defining polynomial z^2 - z - 1
            with tau = 1.618033988749895?

        """
        from sage.groups.affine_gps.affine_group import AffineGroup
        from sage.matrix.special import identity_matrix
        dimension = module.dimension()
        ring = module.base_ring()
        F = AffineGroup(dimension, ring)
        M = multiplier * identity_matrix(dimension)

        edges = [(j,i,F(M, translation)) 
                 for (i,j),L in displacement_matrix.items() 
                 for translation in L]

        return GraphDirectedIteratedFunctionSystem(module, edges)

    def to_digraph(self):
        r"""
        EXAMPLES::

            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: F = AffineGroup(1, QQ)
            sage: f1 = F.linear(1/3)
            sage: f2 = F(1/3, vector([2/3]))
            sage: cantor_ifs = GIFS(QQ^1, [(0,0,f1),(0,0,f2)])
            sage: cantor_ifs.to_digraph()
            Looped multi-digraph on 1 vertex

        """
        from sage.graphs.digraph import DiGraph
        edges = [(u,v) for (u,v,f) in self._edges]
        return DiGraph(edges, format='list_of_edges', loops=True,
                multiedges=True)

    def vertices(self):
        r"""
        EXAMPLES::

            sage: F = AffineGroup(1, QQ)
            sage: f1 = F.linear(1/3)
            sage: f2 = F(1/3, vector([2/3]))
            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: cantor_ifs = GIFS(QQ^1, [(0,0,f1),(0,0,f2)])
            sage: cantor_ifs.vertices()
            [0]

        """
        U = [u for (u,v,f) in self._edges]
        V = [v for (u,v,f) in self._edges]
        return sorted(set(U)|set(V))

    def num_vertices(self):
        r"""
        EXAMPLES::

            sage: F = AffineGroup(1, QQ)
            sage: f1 = F.linear(1/3)
            sage: f2 = F(1/3, vector([2/3]))
            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: cantor_ifs = GIFS(QQ^1, [(0,0,f1),(0,0,f2)])
            sage: cantor_ifs.num_vertices()
            1

        """
        return len(self.vertices())

    def galois_conjugate(self):
        r"""
        Return the element-wise Galois conjugate of this GIFS

        INPUT:

        - ``self`` -- an Affine GIFS, defined on a ring where elements have
          a method ``.galois_conjugate`` (e.g., quadratic number field elements)

        EXAMPLES:

        Fibonacci substitution::

            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: m = WordMorphism('a->ab,b->a')
            sage: s = GIFS.from_one_dimensional_substitution(m)
            sage: s.galois_conjugate()
            GIFS defined by 3 maps on Vector space of dimension 1 over
            Number Field in root with defining polynomial x^2 - x - 1 with
            root = 1.618033988749895?

        Direct Product of 2 Fibonacci::

            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: from slabbe import Substitution2d
            sage: m = Substitution2d(...) # not tested
            sage: s = GIFS.from_two_dimensional_substitution(m)
            sage: s.galois_conjugate()
            Inflation Rule (multiplier=-tau + 1) defined on 4 sets

        ::

            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: z = polygen(QQ, 'z')
            sage: K = NumberField(z**2-z-1, 'tau', embedding=RR(1.6))
            sage: tau = K.gen()
            sage: import itertools
            sage: d = {(i,j):[] for i,j in itertools.product(range(4),repeat=2)}
            sage: d[(0,3)] = [vector(K, (tau,tau))]
            sage: d[(1,2)] = d[(1,3)] = [vector(K, (0,tau))]
            sage: d[(2,1)] = d[(2,3)] = [vector(K, (tau,0))]
            sage: d[(3,0)] = d[(3,1)] = d[(3,2)] = d[(3,3)] = [vector(K, (0,0))]
            sage: ifs = GIFS.from_inflation_rule(K^2, tau, d)
            sage: ifs.galois_conjugate()
            GIFS defined by 9 maps on Vector space of dimension 2 over
            Number Field in tau with defining polynomial z^2 - z - 1 with
            tau = 1.618033988749895?

        """
        edges = [(u,v,galois_conjugate(f)) for (u,v,f) in self._edges]
        return GraphDirectedIteratedFunctionSystem(self._module, edges)

    def __call__(self, S=None, n_iterations=1):
        r"""
        Return the image of the list of list of points.

        INPUT:

        - ``S`` -- list or dict, list of list of points or dictionary
          associating a list of points to each vertex. If a list is used,
          we assume the vertices are integers 0,1,...,n-1.
        - ``n_iterations`` -- integer (default: ``1``)

        EXAMPLES::

            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: F = AffineGroup(1, QQ)
            sage: f1 = F.linear(1/3)
            sage: f2 = F(1/3, vector([2/3]))
            sage: cantor_ifs = GIFS(QQ^1, [(0,0,f1),(0,0,f2)])
            sage: cantor_ifs({0:[vector([0])]})
            {0: [(0), (2/3)]}
            sage: cantor_ifs(_)
            {0: [(0), (2/9), (2/3), (8/9)]}
            sage: cantor_ifs(_)
            {0: [(0), (2/27), (2/9), (8/27), (2/3), (20/27), (8/9), (26/27)]}
            sage: cantor_ifs(_)
            {0: [(0),
              (2/81),
              (2/27),
              (8/81),
              (2/9),
              (20/81),
              (8/27),
              (26/81),
              (2/3),
              (56/81),
              (20/27),
              (62/81),
              (8/9),
              (74/81),
              (26/27),
              (80/81)]}

        ::

            sage: cantor_ifs([[vector([0])]], 2)
            {0: [(0), (2/9), (2/3), (8/9)]}
            sage: cantor_ifs([[vector([0])]], 3)
            {0: [(0), (2/27), (2/9), (8/27), (2/3), (20/27), (8/9), (26/27)]}

        ::

            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: z = polygen(QQ, 'z')
            sage: K = NumberField(z**2-z-1, 'tau', embedding=RR(1.6))
            sage: tau = K.gen()
            sage: import itertools
            sage: d = {(i,j):[] for i,j in itertools.product(range(4),repeat=2)}
            sage: d[(0,3)] = [vector(K, (tau,tau))]
            sage: d[(1,2)] = d[(1,3)] = [vector(K, (0,tau))]
            sage: d[(2,1)] = d[(2,3)] = [vector(K, (tau,0))]
            sage: d[(3,0)] = d[(3,1)] = d[(3,2)] = d[(3,3)] = [vector(K, (0,0))]
            sage: ifs = GIFS.from_inflation_rule(K^2, tau, d)
            sage: ifs(n_iterations=1)
            {0: [], 1: [], 2: [], 3: [(0, 0)]}
            sage: ifs(n_iterations=2)
            {0: [(tau, tau)], 1: [(0, tau)], 2: [(tau, 0)], 3: [(0, 0)]}
            sage: ifs(n_iterations=3)
            {0: [(tau, tau)],
             1: [(tau + 1, tau), (0, tau)],
             2: [(tau, tau + 1), (tau, 0)],
             3: [(tau + 1, tau + 1), (0, tau + 1), (tau + 1, 0), (0, 0)]}
            sage: ifs(n_iterations=4)
            {0: [(3*tau + 1, 3*tau + 1), (tau, 3*tau + 1), (3*tau + 1, tau), (tau, tau)],
             1: [(tau + 1, 3*tau + 1),
              (tau + 1, tau),
              (2*tau + 1, 3*tau + 1),
              (0, 3*tau + 1),
              (2*tau + 1, tau),
              (0, tau)],
             2: [(3*tau + 1, tau + 1),
              (tau, tau + 1),
              (3*tau + 1, 2*tau + 1),
              (tau, 2*tau + 1),
              (3*tau + 1, 0),
              (tau, 0)],
             3: [(tau + 1, tau + 1),
              (2*tau + 1, tau + 1),
              (0, tau + 1),
              (tau + 1, 2*tau + 1),
              (tau + 1, 0),
              (2*tau + 1, 2*tau + 1),
              (0, 2*tau + 1),
              (2*tau + 1, 0),
              (0, 0)]}

        TESTS::

            sage: cantor_ifs([[vector([0])],[vector([0])]])
            Traceback (most recent call last):
            ...
            ValueError: size of input (=2) must match the number of vertices of this GIFS (=1)

        """
        # input S
        if S is None:
            zero = self._module.zero()
            S = {0:[zero]}
        elif isinstance(S, list):
            if not len(S) == self.num_vertices():
                raise ValueError("size of input (={}) must match the number of"
                        " vertices of this GIFS (={})".format(len(S),
                            self.num_vertices()))
            S = dict(enumerate(S))

        # one iteration
        if n_iterations == 1:
            S_image = {}
            for v in self.vertices():
                Ev = [(u,v_,f) for (u,v_,f) in self._edges if v_ == v]
                S_image[v] = [f(p) for (u,_,f) in Ev for p in S.get(u,[])]
            return S_image

        # many iterations
        if n_iterations > 1:
            for _ in range(n_iterations):
                S = self(S, n_iterations=1)
            return S

        # invalid number of iterations
        if n_iterations < 1:
            raise ValueError('n_iterations(={}) must be larger or equal to'
                    ' 1'.format(n_iterations))

    def plot(self, S=None, n_iterations=1, projection=None):
        r"""
        Return a graphic image of the IFS after few iterations

        INPUT:

        - ``S`` -- list or dict, list of list of points or dictionary
          associating a list of points to each vertex. If a list is used,
          we assume the vertices are integers 0,1,...,n-1.
        - ``n_iterations`` -- integer (default: ``1``)
        - ``projection`` -- matrix (default: ``None``), projection matrix
          to 2-dimensional space

        OUTPUT:

        Graphics object

        EXAMPLES:

        The Cantor set::

            sage: F = AffineGroup(1, QQ)
            sage: f1 = F.linear(1/3)
            sage: f2 = F(1/3, vector([2/3]))
            sage: from slabbe import GraphDirectedIteratedFunctionSystem as GIFS
            sage: cantor_ifs = GIFS(QQ^1, [(0,0,f1),(0,0,f2)])
            sage: G = cantor_ifs.plot(n_iterations=7)

        Projection on the vertical y-axis::

            sage: G = cantor_ifs.plot(n_iterations=7, projection=matrix(2,[0,1]))
        
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
        from sage.matrix.constructor import matrix
        from sage.plot.colors import rainbow
        from sage.plot.graphics import Graphics
        from sage.plot.point import points

        if self._module.dimension() == 1 and projection is None:
            # default projection on the x-axis
            projection = matrix([[1],[0]])
        elif self._module.dimension() != 2 and projection is None:
            raise ValueError('a projection matrix must be provided'
                    ' when the dimension of the GIFS (={}) is not'
                    ' 2'.format(self._module.dimension()))

        G = Graphics()
        vertex_to_color = dict(zip(self.vertices(), rainbow(self.num_vertices())))
        ifs = self(S=S, n_iterations=n_iterations)
        for v,P in ifs.items():
            if not self._module.dimension() == 2:
                P = [projection*p for p in P]
            G += points(P, color=vertex_to_color[v], legend_label=str(v))
        return G

def galois_conjugate(f):
    r"""
    Return the element-wise Galois conjugate of an element of an affine
    group 

    INPUT:

    - ``f`` -- affine group element

    EXAMPLES::

        sage: from slabbe.graph_directed_IFS import galois_conjugate
        sage: z = polygen(QQ, 'z')
        sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
        sage: phi = K.gen()
        sage: F = AffineGroup(2, K)
        sage: f = F(phi*identity_matrix(2), (phi,0))
        sage: galois_conjugate(f)
              [-phi + 1        0]     [-phi + 1]
        x |-> [       0 -phi + 1] x + [       0]

    """
    from sage.matrix.constructor import matrix
    F = f.parent()
    dim = F.degree() + 1
    M = matrix(dim,[a.galois_conjugate() for a in f.matrix().list()])
    return F(M)


