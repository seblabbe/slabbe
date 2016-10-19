# -*- coding: utf-8 -*-
r"""
Simultaneous diophantine approximation

EXAMPLES:

The code gets between 1000x faster and 10000x faster compared to the same
function in Python code (diophantine_approximation.py)::

    sage: from slabbe.diophantine_approx_pyx import simultaneous_dioph_approx_pyx
    sage: from slabbe.diophantine_approx import simultaneous_dioph_approx
    sage: simultaneous_dioph_approx_pyx([e,pi], 203)     # 493 µs
    ([1843, 2130, 678], 203.2394429340719)
    sage: simultaneous_dioph_approx([e,pi], 203)         # 905 ms
    ((1843, 2130, 678), 203.239442934072)
    sage: 905/493. * 1000
    1835.69979716024

::

    sage: simultaneous_dioph_approx_pyx([e,pi], 204)     # 2.25 ms
    ([51892, 59973, 19090], 266.16775094991226)
    sage: simultaneous_dioph_approx([e,pi], 204)         # 25s (not tested)
    ((51892, 59973, 19090), 266.167750949912)
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

def simultaneous_dioph_approx_pyx(v, int Q, int start=1, int stop=-1):
    r"""
    Return a simultaneous diophantine approximation of vector ``v`` at distance
    ``1/Q``.

    It searches for all possibilities of denominators in the interval [1, Q^d]
    between start and stop.

    INPUT:

    - ``v`` -- list of real numbers to approximate
    - ``Q`` -- real number, Q>1
    - ``start`` -- integer (default: ``1``), starting value to check
    - ``stop`` -- integer (default: ``-1``), stoping value to check

    OUTPUT:

    - A tuple (u, R) where u=(p_1, ..., p_n, q) is a vector, q is an
      integer, R is a real number such that coordinates of vector p/q are
      at most 1/R from the coordinates of v.

    EXAMPLES::

        sage: from slabbe.diophantine_approx_pyx import simultaneous_dioph_approx_pyx
        sage: simultaneous_dioph_approx_pyx([e,pi], 2)
        (3, 3, 1), 3.54964677830384)
        sage: simultaneous_dioph_approx_pyx([e,pi], 4)
        (19, 22, 7), 35.7490143326079)
        sage: simultaneous_dioph_approx_pyx([e,pi], 35)
        (19, 22, 7), 35.7490143326079)
        sage: simultaneous_dioph_approx_pyx([e,pi], 36)
        (1843, 2130, 678), 203.239442934072)
        sage: simultaneous_dioph_approx_pyx([e,pi], 203)
        (1843, 2130, 678), 203.239442934072)

    We can start the next computation at step 678::

        sage: simultaneous_dioph_approx_pyx([e,pi], 204, start=678)
        (51892, 59973, 19090), 266.167750949912)

    TESTS::

        sage: simultaneous_dioph_approx_pyx([1,e,pi], 1)
        Traceback (most recent call last):
        ...
        ValueError: argument Q(=1) must be > 1

    If the interval is too short, then noting is found::

        sage: simultaneous_dioph_approx_pyx([e,pi], 204, start=500, stop=600)
        (None, None)
    """
    if not Q > 1:
        raise ValueError("argument Q(={}) must be > 1".format(Q))
    cdef int d = len(v)
    cdef double Qinv = 1. / Q
    cdef double un_moins_Qinv = 1 - Qinv
    cdef double a, error
    cdef list p, q_v, frac_q_v
    cdef double* vdouble = <double*>check_allocarray(d, sizeof(double))
    cdef int q
    for i in range(d):
        vdouble[i] = v[i]
    if stop == -1:
        stop = Q**d
    else:
        stop = min(stop, Q**d)
    for q in range(start, stop):
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
            sig_free(vdouble)
            return tuple(p), 1./error
    sig_free(vdouble)
    return None, None

