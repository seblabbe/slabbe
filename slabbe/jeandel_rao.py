# -*- coding: utf-8 -*-
r"""
Random generation of Jeandel-Rao tilings.

EXAMPLES::

    sage: from slabbe import random_jeandel_rao_tiling_rectangle
    sage: tiling = random_jeandel_rao_tiling_rectangle(4, 4)
    sage: tiling
    A wang tiling of a 4 x 4 rectangle
    sage: tiling.table()   # random
    [[1, 10, 4, 5], [1, 3, 3, 7], [0, 9, 10, 4], [0, 9, 3, 3]]

REFERENCES:

.. [Lab2019] S. Labbé. A Markov partition for Jeandel-Rao aperiodic
   Wang tilings. March 2019. https://arxiv.org/abs/1903.06137
"""
#*****************************************************************************
#       Copyright (C) 2017-2019 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
import math
from random import random

def plane_to_torus(m,n):
    r"""
    EXAMPLES::

        sage: from slabbe.jeandel_rao import plane_to_torus
        sage: plane_to_torus(0, 0)                       # abs tol 1e-10
        (0.000000000000000, 0.0)
        sage: plane_to_torus(0.324, .324)                # abs tol 1e-10
        (0.324000000000000, 0.324)
        sage: plane_to_torus(12.324, 12.324)             # abs tol 1e-10
        (0.615796067500630, 3.08793202250021)
        sage: plane_to_torus(100, 100)                   # abs tol 1e-10
        (1.3343685400050447, 3.021286236252207)

    """
    phi = .5 + .5 * math.sqrt(5)
    m = float(m)
    n = float(n)
    h = n // float(phi + 3)
    return (m - h) % float(phi), n % float(phi+3)

def torus_to_code(x,y):
    r"""
    Return in which atom of the partition associated to Jeandel-Rao tilings
    the point (x,y) falls in according to [Lab2019]_.

    EXAMPLES::

        sage: from slabbe.jeandel_rao import torus_to_code
        sage: torus_to_code(0,0)
        0
        sage: torus_to_code(0.23,3.5)
        5
        sage: torus_to_code(1.23,2.243)
        3

    ::

        sage: from slabbe.jeandel_rao import plane_to_torus, random_torus_point
        sage: torus_to_code(*plane_to_torus(14.4141, 89.14))
        9
        sage: torus_to_code(*random_torus_point())                  # random
        3

    """
    phi = .5 + .5 * math.sqrt(5)
    assert x >= 0, "x(={}) must be nonnegative".format(x)
    assert y >= 0, "y(={}) must be nonnegative".format(y)
    assert x < phi, "x(={}) must be less then phi".format(x)
    assert y < phi+3, "y(={}) must be less then phi+3".format(y)
    # test
    if y < 1:
        if x < (1./phi):
            if y <= phi*x:
                return 0
            else:
                return 1
        elif x < 1:
            if y <= (phi**2)*x -phi:
                return 0
            else:
                return 1
        else:
            if y <= phi*x - phi:
                return 0
            else:
                return 1
    if y <= phi*x + 1:
        if x <= 1./phi:
            return 9
        else:
            if y <= (phi**2)*x + 1 - phi:
                if x >= (1./phi) and x < 1:
                    return 9
                elif x >= 1:
                    if y <= phi*x + 1 - phi:
                        return 9
                    else:
                        return 3
            else:
                return 8
    if y > (phi*x + 1) and y <= phi*x + 2:
        if x <= 1:
            if y <= (phi**2)*x + 1:
                return 10
            else:
                return 7
        else:
            return 7
    if y > phi*x +2:
        if x <= 1:
            if y <= (phi**2)*x + 2:
                return 4
            if y <= (phi*x + phi +1) and y > (phi**2)*x + 2:
                return 2
            else:
                if y > (phi**2)*x + phi + 2:
                    return 6
                if y > (phi*x + 3) and x > (1./(phi**2)):
                    return 6
                else:
                    return 5
        else:
            return 6

def random_torus_point():
    r"""
    Return a random point in the rectangle `[0,\phi[\times[0,\phi+3[`.

    EXAMPLES::

        sage: from slabbe.jeandel_rao import random_torus_point
        sage: random_torus_point()                  # random
        (0.947478386174632, 2.62013791669977)
        sage: random_torus_point()                  # random
        (0.568010404619112, 0.933319012345482)
        sage: random_torus_point()                  # random
        (1.06782191679796, 4.58930423801758)

    ::

        sage: from slabbe.jeandel_rao import torus_to_code
        sage: torus_to_code(*random_torus_point())                  # random
        3
        sage: torus_to_code(*random_torus_point())                  # random
        7

    """
    phi = .5 + .5 * math.sqrt(5)
    return random() * phi, random() * (phi+3)

def random_jeandel_rao_tiling_rectangle(width, height, start=None):
    r"""
    Returns a jeandel rao tiling of a rectangle associated to a given
    (random) starting position on the torus.

    INPUT:

    - ``width`` -- integer
    - ``height`` -- integer
    - ``start`` -- pair of real numbers (default:``None``), if ``None``
      a random start point is chosen

    EXAMPLES::

        sage: from slabbe.jeandel_rao import random_jeandel_rao_tiling_rectangle
        sage: tiling = random_jeandel_rao_tiling_rectangle(4,4)
        sage: tiling
        A wang tiling of a 4 x 4 rectangle
        sage: tiling.table()   # random
        [[1, 10, 4, 5], [1, 3, 3, 7], [0, 9, 10, 4], [0, 9, 3, 3]]

    """
    if start is None:
        x0,y0 = random_torus_point()
    else: 
        x0,y0 = start
        x0,y0 = plane_to_torus(x0, y0)
    tiling = [[torus_to_code(*plane_to_torus(x0+a,y0+b)) 
                  for b in range(height)]
                  for a in range(width)]
    tiles = [(2,4,2,1), (2,2,2,0), (1,1,3,1), (1,2,3,2), (3,1,3,3), (0,1,3,1), 
             (0,0,0,1), (3,1,0,2), (0,2,1,2), (1,2,1,4), (3,3,1,2)]
    tiles = [map(str,t) for t in tiles]
    from collections import defaultdict
    color = defaultdict(lambda : 'white')
    color.update({0:'white', 1:'red', 2:'cyan', 3:'green', 4:'lightgray'})
    color.update({str(k):v for k,v in color.items()})
    from slabbe.wang_tiles import WangTiling
    return WangTiling(tiling, tiles=tiles, color=color)

