# -*- coding: utf-8 -*-
r"""
Simultaneous diophantine approximation

EXAMPLES:

The code gets between 1000x faster and 10000x faster compared to the same
function in Python code (diophantine_approximation.py)::

    sage: from slabbe.diophantine_approx_pyx import good_simultaneous_convergents_upto
    sage: from slabbe.diophantine_approx import _best_simultaneous_convergents_upto
    sage: good_simultaneous_convergents_upto([e,pi], 203)     # 493 µs
    [((1843, 2130, 678), 203.24),
     ((51892, 59973, 19090), 266.17),
     ((53735, 62103, 19768), 207.68)]
    sage: _best_simultaneous_convergents_upto([e,pi], 203)         # 905 ms
    ((1843, 2130, 678), 203.239442934072)
    sage: 905/493. * 1000
    1835.69979716024

::

    sage: good_simultaneous_convergents_upto([e,pi], 204)     # 2.25 ms
    [((51892, 59973, 19090), 266.17),
     ((53735, 62103, 19768), 207.68),
     ((113018, 130618, 41577), 279.19)]
    sage: _best_simultaneous_convergents_upto([e,pi], 204)         # 25s (not tested)
    ((51892, 59973, 19090), 266.17)
    sage: 25 / 2.25 * 1000
    11111.1111111111

AUTHORS:

- Sébastien Labbé, October 19, 2016
"""
#*****************************************************************************
#       Copyright (C) 2016 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from math import floor
from builtins import round

include "cysignals/signals.pxi"   # ctrl-c interrupt block support
include "cysignals/memory.pxi"

def good_simultaneous_convergents_upto(v, double Q, int start=1, int step=1):
    r"""
    Return a list of potential best simultaneous diophantine approximation of
    vector ``v`` at distance ``1/Q``.

    It searches for all possibilities of denominators in the interval [1, Q^d]
    starting at start by step. Argument step is used for parallel
    computations purposes.

    INPUT:

    - ``v`` -- list of real numbers to approximate
    - ``Q`` -- real number, Q>1
    - ``start`` -- integer (default: ``1``), starting value to check
    - ``step`` -- integer (default: ``1``), step

    OUTPUT:

    - A tuple (u, R) where u=(p_1, ..., p_n, q) is a vector, q is an
      integer, R is a real number such that coordinates of vector p/q are
      at most 1/R from the coordinates of v.

    EXAMPLES::

        sage: from slabbe.diophantine_approx_pyx import good_simultaneous_convergents_upto
        sage: good_simultaneous_convergents_upto([e,pi], 2)
        [((3, 3, 1), 3.55), ((5, 6, 2), 2.29), ((8, 9, 3), 2.35), ((11, 13, 4), 2.31)]
        sage: good_simultaneous_convergents_upto([e,pi], 4)
        [((19, 22, 7), 35.75), ((38, 44, 14), 17.87), ((41, 47, 15), 4.43)]
        sage: good_simultaneous_convergents_upto([e,pi], 35)
        [((19, 22, 7), 35.75), ((1843, 2130, 678), 203.24), ((1862, 2152, 685), 43.38)]
        sage: good_simultaneous_convergents_upto([e,pi], 36)
        [((1843, 2130, 678), 203.24), ((1862, 2152, 685), 43.38)]
        sage: good_simultaneous_convergents_upto([e,pi], 204)
        [((51892, 59973, 19090), 266.17),
         ((53735, 62103, 19768), 207.68),
         ((113018, 130618, 41577), 279.19)]
        sage: good_simultaneous_convergents_upto([e,pi], 203)
        [((1843, 2130, 678), 203.24),
         ((51892, 59973, 19090), 266.17),
         ((53735, 62103, 19768), 207.68)]

    We can start the next computation at step 678::

        sage: good_simultaneous_convergents_upto([e,pi], 204, start=678)
        [((51892, 59973, 19090), 266.17),
         ((53735, 62103, 19768), 207.68),
         ((113018, 130618, 41577), 279.19)]

    Use start and step for parallel computation purposes::

        sage: good_simultaneous_convergents_upto([e,pi], 2935, start=1)
        [((1174662, 1357589, 432134), 2935.31),
         ((3970510, 4588827, 1460669), 3654.3),
         ((21640489, 25010505, 7961091), 6257.09),
         ((22815151, 26368094, 8393225), 3013.59)]
        sage: good_simultaneous_convergents_upto([e,pi], 2935, start=1, step=3)
        []
        sage: good_simultaneous_convergents_upto([e,pi], 2935, start=2, step=3)
        [((1174662, 1357589, 432134), 2935.31),
         ((3970510, 4588827, 1460669), 3654.3),
         ((22815151, 26368094, 8393225), 3013.59)]
        sage: good_simultaneous_convergents_upto([e,pi], 2935, start=3, step=3)
        [((21640489, 25010505, 7961091), 6257.09)]

    TESTS::

        sage: good_simultaneous_convergents_upto([1,e,pi], 1)
        Traceback (most recent call last):
        ...
        ValueError: argument Q(=1.0) must be > 1

    If the interval is too short, then noting is found::

        sage: good_simultaneous_convergents_upto([e,pi], 204, start=42000, step=2)
        []

    Problem, no solution is found here, why? precision problems?::

        sage: good_simultaneous_convergents_upto([e,pi], 136500, start=1)  # not tested (2min)
        Traceback (most recent call last):
        ...
        ValueError: no solution found
    """
    if not Q > 1:
        raise ValueError("argument Q(={}) must be > 1".format(Q))
    cdef int d = len(v)
    cdef double Qinv = 1. / Q
    cdef double un_moins_Qinv = 1 - Qinv
    cdef double a, error
    cdef list p, q_v, frac_q_v, L
    cdef double* vdouble = <double*>check_allocarray(d, sizeof(double))
    cdef int i,q
    cdef int stop = int(floor(Q**d)) + 1
    for i in range(d):
        vdouble[i] = v[i]
    # for q in range(start, stop, step):
    # because of the use of step, the range is not translated into a C for-loop
    # http://stackoverflow.com/questions/21382180/cython-pure-c-loop-optimization
    L = []
    for q from start <= q < stop by step:
        sig_check() # Check for Keyboard interupt
        for i in range(d):
            a = q*vdouble[i]
            a = a - floor(a) 
            if not (a <= Qinv or un_moins_Qinv <= a):
                break
        else:
            q_v = [q*vdouble[i] for i in range(d)]
            p = [int(round(a)) for a in q_v]
            p.append(q)
            frac_q_v = [a-floor(a) for a in q_v]
            error = max((a if a < .5 else 1-a) for a in frac_q_v)
            L.append((tuple(p), round(1./error, 2)))
    sig_free(vdouble)
    if not L and step == 1:
        # If it was a full search, we should not reach that point
        # We assume that start was set so that there must be solution
        raise ValueError("Did not find an approximation p,q to v={} s.t.  |p-qv|<=1/Q"
            " where Q={}".format(v,Q))
    return L

