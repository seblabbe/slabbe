# -*- coding: utf-8 -*-
r"""
Matrix cocyles

EXAMPLES::

    ...

TODO:

    - Discrete geometry code should use this file

"""
#*****************************************************************************
#       Copyright (C) 2014-2015 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.matrix.constructor import matrix

def projection_matrix(dim_from=3, dim_to=2):
    r"""
    Return a projection matrix from R^d to R^l.

    INPUT:

    - ``dim_from` -- integer (default: ``3``)
    - ``dim_to` -- integer (default: ``2``)

    OUTPUT:

        matrix

    EXAMPLES::

        sage: from slabbe.matrix_cocycle import projection_matrix
        sage: projection_matrix(3,2)
        [-0.866025403784439  0.866025403784439  0.000000000000000]
        [-0.500000000000000 -0.500000000000000   1.00000000000000]
        sage: projection_matrix(2,3)
        [-0.577350269189626 -0.333333333333333]
        [ 0.577350269189626 -0.333333333333333]
        [ 0.000000000000000  0.666666666666667]
    """
    from math import sqrt
    from sage.rings.real_mpfr import RR
    sqrt3 = sqrt(3)
    if dim_from == 3 and dim_to == 2:
        return matrix(2,[-sqrt3,sqrt3,0,-1,-1,2],ring=RR)/2
    elif dim_from == 2 and dim_to == 3:
        return matrix(3,[-sqrt3,-1,sqrt3,-1,0,2],ring=RR)/3
    else:
        s = "for input dim_from={} and dim_to={}"
        raise NotImplementedError(s.format(dim_from, dim_to))

