# -*- coding: utf-8 -*-
r"""
Simultaneous diophantine approximation

EXAMPLES:

The code gets between 1000x faster and 50000x faster compared to the same
function in Python code (diophantine_approximation.py)::

    sage: from slabbe.diophantine_approx_pyx import good_simultaneous_convergents_upto
    sage: from slabbe.diophantine_approx import _best_simultaneous_convergents_upto
    sage: good_simultaneous_convergents_upto([e,pi], 203)     # 493 µs
    [678, 19090, 19768]
    sage: _best_simultaneous_convergents_upto([e,pi], 203)         # 905 ms
    ((1843, 2130, 678), 203.239442934072)
    sage: 905/493. * 1000
    1835.69979716024

::

    sage: good_simultaneous_convergents_upto([e,pi], 204)    # 493 µs
    [19090, 19768, 41577]
    sage: _best_simultaneous_convergents_upto([e,pi], 204)         # 25s (not tested)
    ((51892, 59973, 19090), 266.17)
    sage: 25 / 493. * 1000000
    50709.9391480730

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
from __future__ import absolute_import, print_function

from libc.math cimport floor

from cysignals.signals cimport sig_check   # ctrl-c interrupt block support
from cysignals.memory cimport check_allocarray, sig_free

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

    - list of integers q such that entries of q*v are at most 1/Q from the
      integer lattice

    EXAMPLES::

        sage: from slabbe.diophantine_approx_pyx import good_simultaneous_convergents_upto
        sage: good_simultaneous_convergents_upto([e,pi], 2)
        [1, 2, 3, 4]
        sage: good_simultaneous_convergents_upto([e,pi], 4)
        [7, 14, 15]
        sage: good_simultaneous_convergents_upto([e,pi], 35)
        [7, 678, 685]
        sage: good_simultaneous_convergents_upto([e,pi], 36)
        [678, 685]
        sage: good_simultaneous_convergents_upto([e,pi], 203)
        [678, 19090, 19768]
        sage: good_simultaneous_convergents_upto([e,pi], 204)
        [19090, 19768, 41577]

    The previous results mean that these integers are good approximations::

        sage: v = [e,pi]
        sage: [(7 * a).n() for a in v]
        [19.0279727992133, 21.9911485751286]
        sage: [(678 * a).n() for a in v]
        [1842.99507969523, 2129.99981913388]
        sage: [(19090 * a).n() for a in v]
        [51892.0001052832, 59973.0037570292]

    Knowing the error of previous results, we can start the next computation at
    step 678::

        sage: min(1/abs(678*a - p) for a,p in zip([e,pi], (1843, 2130))).n()
        203.239442934072
        sage: good_simultaneous_convergents_upto([e,pi], 204, start=678)
        [19090, 19768, 41577]

    Use start and step for parallel computation purposes::

        sage: good_simultaneous_convergents_upto([e,pi], 2935, start=1)
        [432134, 1460669, 7961091, 8393225]
        sage: good_simultaneous_convergents_upto([e,pi], 2935, start=1, step=3)
        []
        sage: good_simultaneous_convergents_upto([e,pi], 2935, start=2, step=3)
        [432134, 1460669, 8393225]
        sage: good_simultaneous_convergents_upto([e,pi], 2935, start=3, step=3)
        [7961091]

    TESTS::

        sage: good_simultaneous_convergents_upto([1,e,pi], 1)
        Traceback (most recent call last):
        ...
        ValueError: argument Q(=1.0) must be > 1

    If the interval is too short, then noting is found::

        sage: good_simultaneous_convergents_upto([e,pi], 204, start=42000, step=2)
        []

    No solution was found when the q was `cdef int` instead of `cdef long`::

        sage: good_simultaneous_convergents_upto([e,pi], 136500, start=1)  # not tested (2min)
        [3111494861]

    This used to give a Value Error::

        sage: good_simultaneous_convergents_upto([e,pi], 102300, start=10^9, step=2) # not tested (1min)
        []
        sage: good_simultaneous_convergents_upto([e,pi], 102300, start=10^9+1, step=2) # not tested (1min)
        [3111494861, 3398281569]

    """
    if not Q > 1:
        raise ValueError("argument Q(={}) must be > 1".format(Q))
    cdef int d = len(v)
    cdef double Qinv = 1. / Q
    cdef double un_moins_Qinv = 1 - Qinv
    cdef double a
    cdef list L
    cdef double* vdouble = <double*>check_allocarray(d, sizeof(double))
    cdef int i
    cdef long q
    cdef long stop = int(floor(Q**d)) + 1
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
            L.append(q)
    sig_free(vdouble)
    return L

