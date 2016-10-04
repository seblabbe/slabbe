r"""
Simultaneous diophantine approximation

EXAMPLES::

    sage: from slabbe.diophantine_approximation import simultaneous_convergents
    sage: it = simultaneous_convergents([e, pi])
    sage: [next(it) for _ in range(5)]          # not tested (60s)
    [((3, 3), 1),
     ((19, 22), 7),
     ((1843, 2130), 678),
     ((51892, 59973), 19090),
     ((113018, 130618), 41577)]

The relation with multidimensional continued fraction algorithms (the above
first 3 approximations appear in the convergents of ARP algorithm, but not the
4th neither the 5th)::

    sage: from slabbe.mult_cont_frac import ARP
    sage: algo = ARP()
    sage: algo.n_matrix((1,e,pi), 3)
    [1 1 1]
    [2 3 3]
    [2 3 4]
    sage: algo.n_matrix((1,e,pi), 6)
    [ 3 12  7]
    [ 8 33 19]
    [ 9 38 22]
    sage: algo.n_matrix((1,e,pi), 16)
    [ 678  600  919]
    [1843 1631 2498]
    [2130 1885 2887]
    sage: algo.n_matrix((1,e,pi), 20)
    [10184  4753 13659]
    [27683 12920 37129]
    [31994 14932 42911]
    sage: algo.n_matrix((1,e,pi), 21)
    [23843 18412 13659]
    [64812 50049 37129]
    [74905 57843 42911]
    sage: algo.n_matrix((1,e,pi), 22)
    [ 23843  42255  37502]
    [ 64812 114861 101941]
    [ 74905 132748 117816]
    sage: algo.n_matrix((1,e,pi), 23)
    [ 66098  42255  79757]
    [179673 114861 216802]
    [207653 132748 250564]

AUTHORS:

- Sébastien Labbé, September 22, 2016

TODO:

- Move this code to cython and see the speed improvement

"""
#*****************************************************************************
#       Copyright (C) 2016 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.functions.other import floor
from sage.misc.functional import round
from sage.modules.free_module_element import vector

def frac(x):
    r"""
    Return the fractional part of real number x.

    Not always perfect...

    EXAMPLES::

        sage: from slabbe.diophantine_approximation import frac
        sage: frac(3.2)
        0.200000000000000
        sage: frac(-3.2)
        0.800000000000000
        sage: frac(pi)
        pi - 3
        sage: frac(pi).n()
        0.141592653589793

    This looks suspicious...::

        sage: frac(pi*10**15).n()
        0.000000000000000
    """
    return x - floor(x)

def distance_to_nearest_integer(x):
    raise NotImplementedError

def simultaneous_diophantine_approximation(v, Q, start=1, verbose=False):
    r"""
    Return a simultaneous diophantine approximation of vector ``v`` at distance
    ``1/Q``.

    INPUT:

    - ``v`` -- list of real numbers
    - ``Q`` -- real number, Q>1
    - ``start`` -- integer (default: ``1``), starting value to check
    - ``verbose`` -- boolean (default: ``False``)

    OUTPUT:

    - A tuple (p, q, R) where p is a vector, q is an integer, R is a real number
      such that coordinates of vector p/q are at most 1/R from the coordinates
      of v.

    EXAMPLES::

        sage: from slabbe.diophantine_approximation import simultaneous_diophantine_approximation
        sage: simultaneous_diophantine_approximation([e,pi], 2)
        ((3, 3), 1, 3.54964677830384)
        sage: simultaneous_diophantine_approximation([e,pi], 4)
        ((19, 22), 7, 35.7490143326079)
        sage: simultaneous_diophantine_approximation([e,pi], 35)
        ((19, 22), 7, 35.7490143326079)
        sage: simultaneous_diophantine_approximation([e,pi], 36)
        ((1843, 2130), 678, 203.239442934072)
        sage: simultaneous_diophantine_approximation([e,pi], 203)   # long time (1s)
        ((1843, 2130), 678, 203.239442934072)

    We can start the next computation at step 678::

        sage: simultaneous_diophantine_approximation([e,pi], 204, start=678) # not tested (25s)
        ((51892, 59973), 19090, 266.167750949912)

    TESTS::

        sage: simultaneous_diophantine_approximation([1,e,pi], 1)
        Traceback (most recent call last):
        ...
        ValueError: argument Q(=1) must be > 1
    """
    if not Q > 1:
        raise ValueError("argument Q(={}) must be > 1".format(Q))
    v = vector(v)
    d = len(v)
    Qinv = 1. / Q
    un_moins_Qinv = 1 - Qinv
    for q in range(start, Q**d):
        q_v = q*v
        frac_q_v = map(frac, q_v)
        if verbose:
            print q,[frac(a).n() for a in q_v]
        if all(a <= Qinv or un_moins_Qinv <= a for a in frac_q_v):
            p = map(round, q_v)
            error = max((a if a < .5 else 1-a) for a in frac_q_v)
            return vector(p), q, ~error.n()
    else:
        raise RuntimeError('Did not find diophantine approximation of vector '
                'v={} with parameter Q={}'.format(v, Q))

def simultaneous_convergents(v):  
    r"""
    Return the sequence of convergents to a vector of real number according to
    Dirichlet theorem on simultaneous approximations.

    INPUT:

    - ``v`` -- list of real numbers

    OUTPUT:

    - iterator

    EXAMPLES::

        sage: from slabbe.diophantine_approximation import simultaneous_convergents
        sage: it = simultaneous_convergents([e, pi])
        sage: next(it)
        ((3, 3), 1)
        sage: next(it)
        ((19, 22), 7)
        sage: next(it)          # long time (1s)
        ((1843, 2130), 678)
        sage: next(it)          # not tested (26s)
        ((51892, 59973), 19090)
        sage: next(it)          # not tested (30s)
        ((113018, 130618), 41577)

    Correspondance with continued fraction when d=1::

        sage: it = simultaneous_convergents([e])
        sage: [next(it) for _ in range(10)]
        [((3), 1),
         ((8), 3),
         ((11), 4),
         ((19), 7),
         ((87), 32),
         ((106), 39),
         ((193), 71),
         ((1264), 465),
         ((1457), 536),
         ((2721), 1001)]
        sage: continued_fraction(e).convergents()[:11].list()
        [2, 3, 8/3, 11/4, 19/7, 87/32, 106/39, 193/71, 1264/465, 1457/536, 2721/1001]
    """
    Q = 2
    start = 1
    while True:
        p,q,Q = simultaneous_diophantine_approximation(v, Q, start)
        yield p,q
        Q = floor(Q) + 1
        start = q

