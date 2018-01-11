# -*- coding: utf-8 -*-
r"""
Combinatorics methods

EXEMPLES::

    sage: from slabbe.combinat import random_composition
    sage: c = random_composition(24,9)
    sage: c        # random
    [1, 4, 2, 2, 3, 5, 4, 2, 1]

::

    sage: from slabbe.combinat import random_simplex_point
    sage: random_simplex_point(3)      # random
    [0.2493321790694003, 0.5353600549544871, 0.21530776597611256]

::

    sage: from slabbe.combinat import random_interior_point
    sage: p = polytopes.hypercube(3)
    sage: p = p + vector([20,0,0])
    sage: random_interior_point(p)          # random
    (19.33174562788114, -0.5428002756082744, -0.3568284089832092)
"""
#*****************************************************************************
#       Copyright (C) 2016 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
from random import randint, random, randrange

def non_uniform_randint(L):
    r"""
    Return a random integer from 0 to len(L)-1 with probabilities
    proportional to the (integral) entries of L.

    INPUT:

    - ``L`` -- list of integers

    EXEMPLES::

        sage: from slabbe.combinat import non_uniform_randint
        sage: non_uniform_randint([2,3,5])    # random
        1
        sage: non_uniform_randint([2,3,5])    # random
        2
        sage: non_uniform_randint([2,3,5])    # random
        2
        sage: from collections import Counter
        sage: L = [non_uniform_randint([2,3,5]) for _ in range(100000)]
        sage: Counter(L)            # random
        Counter({2: 49805, 1: 30228, 0: 19967})
    """
    sum_L = sum(L)
    R = randrange(sum_L)
    for i,a in enumerate(L):
        R -= a
        if R < 0:
            return i
    else:
        assert False, "problem, algorithm should not reach this point"

def random_composition(n, length):
    r"""
    EXEMPLES::

        sage: from slabbe.combinat import random_composition
        sage: random_composition(4,2)   # random
        [1, 3]
        sage: random_composition(4,2)   # random
        [2, 2]
        sage: random_composition(4,2)   # random
        [3, 1]
        sage: c = random_composition(24,9)
        sage: c  # random
        [1, 4, 2, 2, 3, 5, 4, 2, 1]
        sage: sum(c)
        24

    Because this is very slow!!::

        sage: C = Compositions(24, length=9)
        sage: %time C.random_element()     # not tested
        CPU times: user 43.3 s, sys: 31.9 ms, total: 43.3 s
        Wall time: 43.2 s
        [2, 2, 5, 2, 8, 1, 1, 2, 1]
    """
    L = [randint(0,n-length) for _ in range(length-1)]
    L.sort()
    L = [a+i+1 for i,a in enumerate(L)]
    L.insert(0,0)
    L.append(n)
    return [L[i+1]-L[i] for i in range(length)]

def random_simplex_point(d):
    r"""
    Return a random vector of d positive real numbers summing to 1.

    INPUT:

    - ``d`` -- integer

    EXEMPLES::

        sage: from slabbe.combinat import random_simplex_point
        sage: random_simplex_point(7)          # random
        [0.06137280030263492,
         0.08066113584919432,
         0.09019666554921013,
         0.24473802319989957,
         0.41761622259683495,
         0.10043545384643937,
         0.004979698655786735]
        sage: random_simplex_point(2)          # random
        [0.5677654878488222, 0.4322345121511778]
        sage: random_simplex_point(3)          # random
        [0.2493321790694003, 0.5353600549544871, 0.21530776597611256]

    TESTS::

        sage: sum(random_simplex_point(4))
        1.0
        sage: sum(random_simplex_point(7))
        1.0
    """
    L = [random() for _ in range(d-1)]
    L.append(0)
    L.append(1)
    L.sort()
    return [L[i+1]-L[i] for i in range(d)]

def random_interior_point(self, a=10, integer=False):
    r"""
    Return a random interior point of a polytope.

    INPUT:

    - ``a`` -- number, amplitude of random deplacement in the direction of each
      ray.
    - ``integer`` -- bool, whether the output must be with integer
      coordinates

    EXEMPLES::

        sage: from slabbe.combinat import random_interior_point
        sage: p = polytopes.hypercube(3)
        sage: p = p + vector([20,0,0])
        sage: p.center()
        (20, 0, 0)
        sage: random_interior_point(p)     # random
        (19.33174562788114, -0.5428002756082744, -0.3568284089832092)
        sage: random_interior_point(p)     # random
        (20.039169786976075, -0.4121594862234468, -0.05623023234688396)
        sage: random_interior_point(p, integer=True)     # random
        (21, 0, 0)
    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    L = list(self.vertices())
    L.extend(a*random()*ray.vector() for ray in self.rays())
    P = Polyhedron(L)
    return random_interior_point_compact_polytope(P)

def random_interior_point_compact_polytope(self, uniform='simplex', integer=False):
    r"""
    Return a random interior point of a compact polytope.

    INPUT:

    - ``uniform`` -- ``'points'`` (slow) or ``'simplex'`` (fast), whether
      to take the probability uniformly with respect to the set of integral
      points or with respect to the simplexes.
    - ``integer`` -- bool, whether the output must be with integer
      coordinates

    EXEMPLES::

        sage: from slabbe.combinat import random_interior_point_compact_polytope
        sage: p = polytopes.hypercube(3)
        sage: p = p + vector([30,20,10])
        sage: p.center()
        (30, 20, 10)
        sage: random_interior_point_compact_polytope(p)     # random
        (19.33174562788114, -0.5428002756082744, -0.3568284089832092)
        sage: random_interior_point_compact_polytope(p)     # random
        (20.039169786976075, -0.4121594862234468, -0.05623023234688396)
        sage: random_interior_point_compact_polytope(p, integer=True) # random
        (30, 19, 9)
    """
    assert self.is_compact(), "input is not compact"
    from sage.geometry.polyhedron.constructor import Polyhedron
    triangulation = self.triangulate()
    T = []
    for simplex_indices in self.triangulate():
        simplex_vertices = [self.Vrepresentation(i) for i in simplex_indices]
        simplex = Polyhedron(simplex_vertices)
        T.append(simplex)
    if uniform == 'points':
        v = [p.integral_points_count() for p in T] # slow
        i = non_uniform_randint(v)
    elif uniform == 'simplex':
        i = randrange(len(T))
    else:
        raise ValueError('unknown value for uniform(={})'.format(uniform))
    return random_interior_point_simplex(T[i], integer=integer)
def random_interior_point_simplex(self, integer=False):
    r"""
    Return a random interior point of a simplex.

    This method was based on the code ``P.center()`` of sage.

    INPUT:

    - ``integer`` -- bool, whether the output must be with integer
      coordinates

    EXEMPLES::

        sage: from slabbe.combinat import random_interior_point_simplex
        sage: P = 10 * polytopes.simplex(3)
        sage: random_interior_point_simplex(P)          # random
        (2.8787864522849462, 5.302173919578364, 1.7059355910006113, 0.11310403713607808)
        sage: a = random_interior_point_simplex(P, integer=True)
        sage: a             # random
        (0, 7, 1, 2)
        sage: a in P
        True
    """
    from sage.modules.free_module_element import vector
    from sage.misc.functional import round
    from sage.functions.generalized import sign
    assert self.is_simplex(), "input is not a simplex"
    d = self.n_vertices()
    R = random_simplex_point(d)
    vertex_sum = vector(self.base_ring(), [0]*self.ambient_dim())
    for v,r in zip(self.vertex_generator(), R):
        vertex_sum += r*v.vector()
    if integer:
        v = vector(map(round, vertex_sum))
        r = round(sum(vertex_sum) - sum(v))
        s = sign(r)
        for _ in range(abs(r)):
            v[randrange(len(v))] += s
        v.set_immutable()
        return v
    else:
        vertex_sum.set_immutable()
        return vertex_sum

def integral_points_count_union_of_polytopes(L):
    r"""
    Return the cardinality of an union of polytopes.

    See https://en.wikipedia.org/wiki/Inclusion–exclusion_principle

    INPUT:

    - ``L`` -- list of polytopes

    EXEMPLES::

        sage: P = Polyhedron(ieqs=[[0,1,0],[0,0,1],[9,-1,0],[9,0,-1]])
        sage: Q = Polyhedron(ieqs=[[-5,1,0],[-5,0,1],[14,-1,0],[14,0,-1]])
        sage: P.integral_points_count()   # optional -- latte_int
        100
        sage: Q.integral_points_count()   # optional -- latte_int
        100
        sage: from slabbe.combinat import integral_points_count_union_of_polytopes
        sage: integral_points_count_union_of_polytopes([P,Q])  # optional -- latte_int
        175
    """
    from sage.combinat.subset import Subsets
    n = len(L)
    range_n = range(n)
    res = 0
    for k in range_n:
        sgn = (-1)**k
        for indices in Subsets(range_n, k+1):
            I = intersection_of_polytopes(L[i] for i in indices)
            try:
                card = I.integral_points_count()
            except IndexError:
                print("IndexError when indices={}, see https://trac.sagemath.org/ticket/21491".format(indices))
                card = 0
            res += sgn * card
    return res

def intersection_of_polytopes(L):
    r"""
    Return the intersection of a list of polytopes.

    INPUT:

    - ``L`` -- list of polytopes

    EXEMPLES::

        sage: from slabbe.combinat import intersection_of_polytopes
        sage: P = Polyhedron(ieqs=[[0,1,0],[0,0,1],[9,-1,0],[9,0,-1]])
        sage: Q = Polyhedron(ieqs=[[-5,1,0],[-5,0,1],[14,-1,0],[14,0,-1]])
        sage: I = intersection_of_polytopes([P,Q])
        sage: I
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
        sage: I.integral_points_count()      # optional -- latte_int
        25

    TESTS::

        sage: intersection_of_polytopes(iter([P,Q]))
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
        sage: intersection_of_polytopes(iter([]))
        Traceback (most recent call last):
        ...
        NotImplementedError: intersection of an empty list of polytopes not defined
    """
    it = iter(L)
    try:
        rep = next(it)
    except StopIteration:
        raise NotImplementedError("intersection of an empty list of polytopes not defined")
    for p in it:
        rep = rep.intersection(p)
    return rep

