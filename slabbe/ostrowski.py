# -*- coding: utf-8 -*-
r"""
Ostrowski numeration

See [Ber2001]_.

REFERENCES:

.. [Ber2001] Valérie Berthé. Autour du système de numération d’Ostrowski.
   Bull.  Belg. Math. Soc.  Simon Stevin, 8(2):209–239, 2001. Journées
   Montoises d’Informatique Théorique (Marne-la-Vallée, 2000).

.. [Bou2015] Bourla, Avraham. « Irrational Base Counting ».
    arXiv:1511.02179 [math], 6 novembre 2015. http://arxiv.org/abs/1511.02179.

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
from __future__ import absolute_import, print_function

from sage.rings.continued_fraction import continued_fraction
from sage.rings.infinity import Infinity

def ostrowski_integer(n, alpha):
    r"""
    INPUT:

    - ``n`` -- integer >= 0
    - ``alpha`` -- irrational real number > 0

    EXAMPLES::

        sage: from slabbe.ostrowski import ostrowski_integer
        sage: ostrowski_integer(5, golden_ratio)
        ([0, 0, 0, 0, 1], [1, 1, 2, 3, 5])
        sage: ostrowski_integer(10, golden_ratio)
        ([0, 0, 1, 0, 0, 1], [1, 1, 2, 3, 5, 8])
        sage: ostrowski_integer(123456, e)
        ([0, 0, 1, 0, 2, 0, 0, 5, 0, 0, 5, 0, 1, 6],
         [1, 1, 3, 4, 7, 32, 39, 71, 465, 536, 1001, 8544, 9545, 18089])
        sage: ostrowski_integer(123456, pi)
        ([4, 11, 0, 211, 0, 0, 0, 1], 
         [1, 7, 106, 113, 33102, 33215, 66317, 99532])

    TESTS::

        sage: ostrowski_integer(10, 4/5)
        Traceback (most recent call last):
        ...
        ValueError: alpha (=4/5) must be irrational

    ::

        sage: for i in range(10): print(i, ostrowski_integer(i, golden_ratio))
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

    Digits of numbers from 1 to 24 in base sqrt(2)-1 where (q_k)_0^3=(1,2,5,12)
    appearing in [Bou2015]_:: 

        sage: rows = [[i]+ostrowski_integer(i, sqrt(2)-1)[0]+[0,0,0,0] for i in range(25)]
        sage: table(rows=rows,header_row='N c1 c2 c3 c4'.split())
          N    c1   c2   c3   c4
        +----+----+----+----+----+
          0    0    0    0    0
          1    1    0    0    0
          2    0    1    0    0
          3    1    1    0    0
          4    0    2    0    0
          5    0    0    1    0
          6    1    0    1    0
          7    0    1    1    0
          8    1    1    1    0
          9    0    2    1    0
          10   0    0    2    0
          11   1    0    2    0
          12   0    0    0    1
          13   1    0    0    1
          14   0    1    0    1
          15   1    1    0    1
          16   0    2    0    1
          17   0    0    1    1
          18   1    0    1    1
          19   0    1    1    1
          20   1    1    1    1
          21   0    2    1    1
          22   0    0    2    1
          23   1    0    2    1
          24   0    0    0    2

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
        

def ostrowski_real(beta, alpha, stop=10, verbose=False):
    r"""
    this is broken code

    EXAMPLES::

        sage: from slabbe.ostrowski import ostrowski_real
        sage: ostrowski_real(golden_ratio^-2, golden_ratio-1, stop=5)   # not tested
        Traceback (most recent call last):
        ...
        AssertionError: 0 <= b_2(=3) <= a_2(=1) is false

    ::

        sage: ostrowski_real(golden_ratio^-3, golden_ratio-1, stop=5)
        ([0, 0, 1, 0, 0],
        [golden_ratio - 1,
        golden_ratio - 2,
        2*golden_ratio - 3,
        3*golden_ratio - 5,
        5*golden_ratio - 8])
    """
    from sage.functions.other import floor

    alpha_cf = continued_fraction(alpha)
    if alpha_cf.length() < Infinity:
        raise ValueError('alpha (={}) must be irrational'.format(alpha))

    remainder = beta
    Digits = []
    Theta = []
    for i in range(stop):
        theta_i = alpha_cf.q(i)*alpha - alpha_cf.p(i)
        b_i = floor(remainder / theta_i * (-1)**i)
        assert 0 <= b_i <= alpha_cf.quotient(i), ("0 <= b_{i}(={}) <= a_{i}(={})"
                           " is false".format(b_i, alpha_cf.quotient(i),i=i))
        remainder -= b_i * theta_i
        if verbose:
            print(theta_i, b_i, remainder)

        Digits.append(b_i)
        Theta.append(theta_i)

    return Digits, Theta


def cf_positive_representation(beta, alpha):
    pass
    


