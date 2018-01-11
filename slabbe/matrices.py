# -*- coding: utf-8 -*-
r"""
Matrix functions

EXAMPLES::

    ...

TODO:

    - Discrete geometry code should use projection_matrix from here

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
from sage.matrix.constructor import matrix
import heapq

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
    elif dim_from == 4 and dim_to == 2:
        return matrix(2,[-sqrt3,sqrt3,0,1,-1,-1,2,0],ring=RR)/2
    elif dim_from == 4 and dim_to == 3:
        sqrt2 = sqrt(2)
        return matrix([(1,-1,0,0), 
                       (0,0,1,-1),
                       (-1/sqrt2,-1/sqrt2,1/sqrt2,1/sqrt2)],
                       ring=RR)
    else:
        s = "for input dim_from={} and dim_to={}"
        raise NotImplementedError(s.format(dim_from, dim_to))

M3to2 = projection_matrix(3, 2)
M2to3 = projection_matrix(2, 3)
M4to2 = projection_matrix(4, 2)
M4to3 = projection_matrix(4, 3)

def perron_right_eigenvector(M):
    r"""
    EXAMPLES::

        sage: from slabbe.matrices import perron_right_eigenvector
        sage: m = matrix(2,[-11,14,-26,29])
        sage: perron_right_eigenvector(m)    # abs tol 0.0000001
        (15.0, (0.35, 0.65))
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
        sage: m = matrix(4, range(-8,8))
        sage: is_nonnegative(m)
        False
        sage: m = matrix(4, range(16))
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
def is_pisot(M):
    r"""
    EXAMPLES::

        sage: from slabbe.matrices import is_pisot
        sage: M = matrix(2,[1,1,0,1])
        sage: is_pisot(M)
        False

    ::

        sage: M = matrix(2,[0,1,1,1])
        sage: is_pisot(M)
        True
    """
    eig_module = map(abs, M.eigenvalues())
    largest, second_largest = heapq.nlargest(2, eig_module)
    return second_largest < 1

def conjugate_matrix_Z(M):
    r"""
    Return the conjugate matrix Z as defined in [1].

    EXAMPLES::

        sage: from slabbe.matrices import conjugate_matrix_Z
        sage: M = matrix(2, [11,29,14,-1])
        sage: conjugate_matrix_Z(M)       # abs tol 1e-8
        [11.674409930010482  27.69820597163912]
        [14.349386111618157  -1.67440993001048]
        sage: conjugate_matrix_Z(M)^2     # abs tol 1e-8
        [533.7440993001048 276.9820597163913]
        [143.4938611161816 400.2559006998952]

    ::

        sage: M = matrix(2, [-11,14,-26,29])
        sage: conjugate_matrix_Z(M)     # abs tol 1e-8
        [ 7.200000000000004  4.199999999999998]
        [ 7.799999999999995 10.800000000000002]
        sage: conjugate_matrix_Z(M) * 5     # abs tol 1e-8
        [ 36.00000000000002 20.999999999999993]
        [ 38.99999999999998 54.000000000000014]

    ::

        sage: M = matrix(2, [-11,26,-14,29]) / 15
        sage: conjugate_matrix_Z(M)     # abs tol 1e-8
        [ 0.5999999999999999  0.3999999999999999]
        [0.39999999999999986  0.5999999999999999]

    REFERENCES:

    [1] Labbé, Jean-Philippe, et Sébastien Labbé. « A Perron theorem for matrices
    with negative entries and applications to Coxeter groups ». arXiv:1511.04975
    [math], 16 novembre 2015. http://arxiv.org/abs/1511.04975.
    """
    from sage.matrix.special import identity_matrix
    eig, v = perron_right_eigenvector(M)
    v /= sum(v)
    nrows = M.nrows()
    UNtop = matrix([1]*nrows)
    v = matrix.column(v)
    id = identity_matrix(nrows)
    Z = eig*v*UNtop + (id - v*UNtop)*M
    return Z

def recurrence_matrix(coeffs):
    r"""
    Return the recurrence matrix of a relation, for example:

    INPUT: 

    - ``coeffs`` -- list of integers, for example if 
      R(n) = R(n-1) + 2 R(n-2) + 3R(n-3) + 4R(n-4) + 5R(n-5)
      then coeff must be [1,2,3,4,5]

    EXAMPLES::

        sage: from slabbe.matrices import recurrence_matrix
        sage: recurrence_matrix([1,2,3,4,5])
        [1 2 3 4 5]
        [1 0 0 0 0]
        [0 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 1 0]
    """
    poly = [-a for a in reversed(coeffs)]
    poly.append(1)
    return matrix.companion(poly, format='top')

def spectrum(M):
    r"""
    EXAMPLES::

        sage: from slabbe.matrices import spectrum, recurrence_matrix
        sage: M = recurrence_matrix([1,2,3,4,5])
        sage: spectrum(M)
        2.576021761956651?
    """
    return max(map(abs, M.eigenvalues()))

def map_coefficients_to_variable_index(M, x):
    r"""
    INPUT:

    - ``M`` -- matrix
    - ``x`` -- string, variable

    EXAMPLES::

        sage: from slabbe.matrices import map_coefficients_to_variable_index
        sage: M = matrix(2, range(4))
        sage: map_coefficients_to_variable_index(M, 's')
        [s_0 s_1]
        [s_2 s_3]
        sage: latex(_)
        \left(\begin{array}{rr}
        s_{0} & s_{1} \\
        s_{2} & s_{3}
        \end{array}\right)
    """
    from sage.calculus.var import var
    L = [var('{}_{}'.format(x,a)) for a in M.list()]
    return matrix(M.nrows(), L)


