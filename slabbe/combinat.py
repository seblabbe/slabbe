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
from random import randint, random, randrange

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

    This method was based on the code ``P.center()`` of sage.

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
    from sage.modules.free_module_element import vector
    from sage.misc.functional import round
    from sage.functions.generalized import sign
    vertex_sum = vector(self.base_ring(), [0]*self.ambient_dim())
    d = self.n_vertices()
    R = random_simplex_point(d)
    for v,r in zip(self.vertex_generator(), R):
        vertex_sum += r*v.vector()
    for ray in self.rays():
        vertex_sum += a*random()*ray.vector()
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

