# -*- coding: utf-8 -*-
r"""
Matrix functions

EXAMPLES::

    ...

TODO:

    - Discrete geometry code should use projection_matrix from here

"""
#*****************************************************************************
#       Copyright (C) 2014-2016 Sébastien Labbé <slabqc@gmail.com>
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

        sage: from slabbe.matrices import projection_matrix
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

def perron_right_eigenvector(M):
    r"""
    EXAMPLES::

        sage: from slabbe.matrices import perron_right_eigenvector
        sage: m = matrix(2,[-11,14,-26,29])
        sage: perron_right_eigenvector(m)
        (15.0000000000000, (0.35, 0.6499999999999999))
    """
    import numpy
    from sage.rings.real_mpfr import RR
    from sage.rings.all import CC
    from sage.modules.free_module_element import vector
    eig, vec = numpy.linalg.eig(M)
    index = abs(eig).argmax()
    rightv = vec.transpose()[index]
    if eig[index].imag == 0:
        eig_sage = RR(eig[index].real)
        vec_sage = vector(a.real for a in rightv)
    else:
        eig_sage = CC(eig[index])
        vec_sage = vector(CC, rightv)
    return eig_sage, vec_sage/sum(vec_sage)

def is_positive(M):
    r"""
    EXAMPLES::

        sage: from slabbe.matrices import is_positive
        sage: m = matrix(4, range(16))
        sage: is_positive(m)
        False
        sage: m = matrix(4, range(1,17))
        sage: is_positive(m)
        True
    """
    return all(a > 0 for a in M.list())

def is_nonnegative(M):
    r"""
    EXAMPLES::

        sage: from slabbe.matrices import is_nonnegative
        sage: m = matrix(4, range(16))
        sage: is_nonnegative(m)
        False
        sage: m = matrix(4, range(1,17))
        sage: is_nonnegative(m)
        True
    """
    return all(a >= 0 for a in M.list())

def is_primitive(M):
    r"""
    EXAMPLES::

        sage: from slabbe.matrices import is_primitive
        sage: m = matrix(2, [0,1,1,1])
        sage: is_primitive(m)
        True
        sage: m = matrix(2, [1,1,0,1])
        sage: is_primitive(m)
        False
    """
    if not is_nonnegative(M):
        return False
    power = M
    order = 1
    dim = M.nrows()
    max_order = (dim-1)**2 + 1
    while True:
        l = power.list()
        if len(l) == 0:
            return False
        try:
            l.index(0)
        except ValueError:
            return True
        if order > max_order:
            return False
        power *= power
        order += order
