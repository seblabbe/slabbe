# -*- coding: utf-8 -*-
r"""

REFERENCES:

    Valérie Berthé. Autour du système de numération d’Ostrowski. Bull.
    Belg. Math. Soc.  Simon Stevin, 8(2):209–239, 2001. Journées Montoises
    d’Informatique Théorique (Marne-la-Vallée, 2000).

AUTHOR:

    - Sébastien Labbé, May 24, 2017
"""
#*****************************************************************************
#       Copyright (C) 2017 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.rings.continued_fraction import continued_fraction
from sage.rings.infinity import Infinity

def ostrowski(n, alpha):
    r"""
    INPUT:

    - ``n`` -- integer >= 0
    - ``alpha`` -- irrational real number > 0

    EXAMPLES::

        sage: ostrowski(5, golden_ratio)
        ([0, 0, 0, 0, 1], [1, 1, 2, 3, 5])
        sage: ostrowski(10, golden_ratio)
        ([0, 0, 1, 0, 0, 1], [1, 1, 2, 3, 5, 8])
        sage: ostrowski(123456, e)
        ([0, 0, 1, 0, 2, 0, 0, 5, 0, 0, 5, 0, 1, 6],
         [1, 1, 3, 4, 7, 32, 39, 71, 465, 536, 1001, 8544, 9545, 18089])
        sage: ostrowski(123456, pi)
        ([4, 11, 0, 211, 0, 0, 0, 1], 
         [1, 7, 106, 113, 33102, 33215, 66317, 99532])

    TESTS::

        sage: ostrowski(10, 4/5)
        Traceback (most recent call last):
        ...
        ValueError: alpha (=4/5) must be irrational

    ::

        sage: for i in range(10): print i, ostrowski(i, golden_ratio)
        0 ([], [])
        1 ([0, 1], [1, 1])
        2 ([0, 0, 1], [1, 1, 2])
        3 ([0, 0, 0, 1], [1, 1, 2, 3])
        4 ([0, 1, 0, 1], [1, 1, 2, 3])
        5 ([0, 0, 0, 0, 1], [1, 1, 2, 3, 5])
        6 ([0, 1, 0, 0, 1], [1, 1, 2, 3, 5])
        7 ([0, 0, 1, 0, 1], [1, 1, 2, 3, 5])
        8 ([0, 0, 0, 0, 0, 1], [1, 1, 2, 3, 5, 8])
        9 ([0, 1, 0, 0, 0, 1], [1, 1, 2, 3, 5, 8])
    """
    Q = [1]
    alpha_cf = continued_fraction(alpha)
    if alpha_cf.length() < Infinity:
        raise ValueError('alpha (={}) must be irrational'.format(alpha))
    i = 1
    while Q[-1] <= n:
        Q.append(alpha_cf.denominator(i))
        i += 1
    Q.pop()
    remainder = n
    D = []
    for q in reversed(Q):
        d = remainder // q
        D.insert(0, d)
        remainder = remainder % q
    return D,Q
        


    


