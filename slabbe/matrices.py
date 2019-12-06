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

def column_norm_ratio(M, p=1):
    r"""
    Return the maximum of the ratio of the norm of two columns.

    INPUT:

    - ``p`` - default: 2 - ``p`` can be a real number greater than 1,
      infinity (``oo`` or ``Infinity``), or a symbolic expression.

      - `p=1`: the taxicab (Manhattan) norm
      - `p=2`: the usual Euclidean norm (the default)
      - `p=\infty`: the maximum entry (in absolute value)

    EXAMPLES::

        sage: from slabbe.matrices import column_norm_ratio
        sage: M = matrix(3, range(9))
        sage: column_norm_ratio(M)
        5/3

    """
    L = [column.norm(p) for column in M.columns()]
    m = min(L)
    if m == 0:
        return Infinity
    else:
        return max(L) / m

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

def Minkowski_embedding_without_sqrt2(self, B=None, prec=None):
    r"""
    This method is a modification of the ``Minkowski_embedding`` method of
    NumberField in sage (without sqrt2).

    EXAMPLES::

        sage: from slabbe.matrices import Minkowski_embedding_without_sqrt2
        sage: F.<alpha> = NumberField(x^3+2)
        sage: F.minkowski_embedding()
        [ 1.00000000000000 -1.25992104989487  1.58740105196820]
        [ 1.41421356237309 0.890898718140339 -1.12246204830937]
        [0.000000000000000  1.54308184421705  1.94416129723967]
        sage: Minkowski_embedding_without_sqrt2(F)
        [  1.00000000000000  -1.25992104989487   1.58740105196820]
        [  1.00000000000000  0.629960524947437 -0.793700525984099]
        [ 0.000000000000000   1.09112363597172   1.37472963699860]
        sage: Minkowski_embedding_without_sqrt2(F, [1, alpha+2, alpha^2-alpha])
        [ 1.00000000000000 0.740078950105127  2.84732210186307]
        [ 1.00000000000000  2.62996052494744 -1.42366105093154]
        [0.000000000000000  1.09112363597172 0.283606001026881]
        sage: Minkowski_embedding_without_sqrt2(F) * (alpha + 2).vector().column()
        [0.740078950105127]
        [ 2.62996052494744]
        [ 1.09112363597172]

    Tribo::

        sage: F.<beta> = NumberField(x^3-x^2-x-1)
        sage: F.minkowski_embedding()
        [  1.00000000000000   1.83928675521416   3.38297576790624]
        [  1.41421356237309 -0.593465355971987 -0.270804762516626]
        [ 0.000000000000000  0.857424571985895 -0.719625086862932]
        sage: Minkowski_embedding_without_sqrt2(F)
        [  1.00000000000000   1.83928675521416   3.38297576790624]
        [  1.00000000000000 -0.419643377607080 -0.191487883953119]
        [ 0.000000000000000  0.606290729207199 -0.508851778832738]

    Comprendre le problème de norme::

        sage: norme = lambda v:abs(v[0]) * (v[1]^2 + v[2]^2)
        sage: F.<beta> = NumberField(x^3-x^2-x-1)
        sage: M = Minkowski_embedding_without_sqrt2(F)
        sage: norme(M*vector((1,0,0)))
        1.00000000000000
        sage: norme(M*vector((1,0,-1)))
        4.00000000000000

    """
    r,s = self.signature()
    places = self.places(prec=prec)

    if B is None:
        B = [self.gen(0)**i for i in range(self.degree())]

    rows = []
    for i in range(r):
        rows.append([places[i](b) for b in B])
    for i in range(s):
        row_real = []
        row_imag = []
        for b in B:
            z = places[r+i](b)
            row_real.append(z.real())
            row_imag.append(z.imag())
        rows.append(row_real)
        rows.append(row_imag)

    from sage.matrix.constructor import matrix
    return matrix(rows)

def Minkowski_projection_pair(self, B=None, prec=None):
    r"""
    Return the projections to the expanding and contracting spaces.

    OUTPUT:

    - tuple (A, B) of matrices

    EXAMPLES::

        sage: from slabbe.matrices import Minkowski_projection_pair
        sage: F.<alpha> = NumberField(x^3+2)
        sage: Minkowski_projection_pair(F)
        (
        [  1.00000000000000  -1.25992104989487   1.58740105196820]
        [  1.00000000000000  0.629960524947437 -0.793700525984099]
        [ 0.000000000000000   1.09112363597172   1.37472963699860], []
        )
        sage: Minkowski_projection_pair(F, [1, alpha+2, alpha^2-alpha])
        (
        [ 1.00000000000000 0.740078950105127  2.84732210186307]
        [ 1.00000000000000  2.62996052494744 -1.42366105093154]
        [0.000000000000000  1.09112363597172 0.283606001026881], []
        )

    Tribo::

        sage: F.<beta> = NumberField(x^3-x^2-x-1)
        sage: Minkowski_projection_pair(F)
        (
        [1.000000000000000000000000000000 1.839286755214161132551852564671
        3.382975767906237494122708536521],
        [  1.00000000000000 -0.419643377607080 -0.191487883953119]
        [ 0.000000000000000  0.606290729207199 -0.508851778832738]
        )

    """
    r,s = self.signature()
    places = self.places(prec=prec)
    beta = self.gen()

    if B is None:
        B = [beta**i for i in range(self.degree())]

    rows_expanding = []
    rows_contracting = []

    for i in range(r):
        place = places[i]
        row = [place(b) for b in B]
        norm = place(beta).abs()
        if norm < 1:
            rows_contracting.append(row)
        elif norm > 1:
            rows_expanding.append(row)
        else:
            raise NotImplementedError

    for i in range(s):
        place = places[r+i]
        row_real = []
        row_imag = []
        for b in B:
            z = place(b)
            row_real.append(z.real())
            row_imag.append(z.imag())
        norm = place(beta).abs()
        if norm < 1:
            rows_contracting.append(row_real)
            rows_contracting.append(row_imag)
        elif norm > 1:
            rows_expanding.append(row_real)
            rows_expanding.append(row_imag)
        else:
            raise NotImplementedError

    from sage.matrix.constructor import matrix
    return (matrix(len(rows_expanding), self.degree(), rows_expanding),
            matrix(len(rows_contracting), self.degree(), rows_contracting))

def rauzy_projection(M, beta=None, prec=53):
    r"""
    Returns a projection matrix of the canonical basis using the Minkowski
    embedding associated to the left eigenvector of the given eigenvalue.

    INPUT:

    - ``beta`` - a real element of ``QQbar`` of degree >= 2 (default:
      ``None``).  The eigenvalue used for the projection.  It must be an
      eigenvalue of ``M``.  The one used by default is the maximal
      eigenvalue of ``M`` (usually a Pisot number), but matrices of order
      larger than 3 letters other interesting choices are sometimes
      possible.

    - ``prec`` - integer (default: ``53``), the number of bits used in the
      floating point representations of the coordinates.

    OUTPUT:

        matrix

    EXAMPLES:

    Fibonacci::

        sage: from slabbe.matrices import rauzy_projection
        sage: m = matrix(2,(1,1,1,0))
        sage: m
        [1 1]
        [1 0]
        sage: rauzy_projection(m)
        [ 1.000000000000000000000000000000 -1.618033988749894848204586834366]
        [ 1.000000000000000000000000000000 0.6180339887498948482045868343656]

    Tribonacci::

        sage: m = matrix(3, [1,1,1, 1,0,0, 0,1,0])
        sage: rauzy_projection(m)
        [  1.00000000000000  0.839286755214161  0.543689012692076]
        [  1.00000000000000  -1.41964337760708 -0.771844506346038]
        [ 0.000000000000000  0.606290729207199  -1.11514250803994]
        sage: matrix(2,(0,1,0, 0,0,-1))*rauzy_projection(m)
        [  1.00000000000000  -1.41964337760708 -0.771844506346038]
        [ 0.000000000000000 -0.606290729207199   1.11514250803994]

    which corresponds to the Rauzy fractal projection coded by Timo::

        sage: s = WordMorphism('1->12,2->13,3->1')
        sage: s.rauzy_fractal_projection()
        {'1': (1.00000000000000, 0.000000000000000),
         '2': (-1.41964337760708, -0.606290729207199),
         '3': (-0.771844506346038, 1.11514250803994)}

    TESTS::

        sage: t = WordMorphism('1->12,2->3,3->45,4->5,5->6,6->7,7->8,8->1')
        sage: m = matrix(t)
        sage: rauzy_projection(m).T
        [  1.00000000000000   1.00000000000000  0.000000000000000]
        [ 0.324717957244746  -1.66235897862237  0.562279512062301]
        [ 0.430159709001947  0.784920145499027  -1.30714127868205]
        [ 0.245122333753307   1.87743883312335  0.744861766619744]
        [ 0.324717957244746  -1.66235897862237  0.562279512062301]
        [ 0.430159709001947  0.784920145499027  -1.30714127868205]
        [ 0.569840290998053  0.215079854500973   1.30714127868205]
        [ 0.754877666246693 -0.877438833123346 -0.744861766619744]
        sage: t.rauzy_fractal_projection()
        {'1': (1.00000000000000, 0.000000000000000),
         '2': (-1.66235897862237, -0.562279512062301),
         '3': (0.784920145499027, 1.30714127868205),
         '4': (1.87743883312335, -0.744861766619744),
         '5': (-1.66235897862237, -0.562279512062301),
         '6': (0.784920145499027, 1.30714127868205),
         '7': (0.215079854500973, -1.30714127868205),
         '8': (-0.877438833123346, 0.744861766619744)}

    ::

        sage: E = t.incidence_matrix().eigenvalues()
        sage: x = [x for x in E if -0.8 < x < -0.7][0]
        sage: x
        -0.7548776662466928?
        sage: rauzy_projection(m, beta=x).T
        [  1.00000000000000   1.00000000000000  0.000000000000000]
        [ -1.75487766624669 -0.122561166876654  0.744861766619744]
        [  1.32471795724475 -0.662358978622373  0.562279512062301]
        [ -4.07959562349144 -0.460202188254281  0.182582254557443]
        [  3.07959562349144 -0.539797811745719 -0.182582254557443]
        [ -2.32471795724475 -0.337641021377627 -0.562279512062301]
        [  1.75487766624669  0.122561166876654 -0.744861766619744]
        [ -1.32471795724475  0.662358978622373 -0.562279512062301]
        sage: t.rauzy_fractal_projection(eig=x)
        {'1': (1.00000000000000, 0.000000000000000),
         '2': (-0.122561166876654, -0.744861766619744),
         '3': (-0.662358978622373, -0.562279512062301),
         '4': (-0.460202188254281, -0.182582254557443),
         '5': (-0.539797811745719, 0.182582254557443),
         '6': (-0.337641021377627, 0.562279512062301),
         '7': (0.122561166876654, 0.744861766619744),
         '8': (0.662358978622373, 0.562279512062301)}

    AUTHORS:

     - Timo Jolivet (2012-06-16) -- for substitutions in Sage
     - Sébastien Labbé (2018-03-08) -- for matrices, using Minkowski
       embedding
    """
    # Eigenvalue
    if beta is None:
        beta = max(M.eigenvalues(), key=abs)

    # One possibility is to do:
    # K,elt,hom = beta.as_number_field_element()

    # Left eigenvector vb in the number field Q(beta)
    from sage.rings.number_field.number_field import NumberField
    K = NumberField(beta.minpoly(), 'b')
    vb = (M-K.gen()).kernel().basis()[0]

    return Minkowski_embedding_without_sqrt2(K, vb)
    #return K.Minkowski_embedding(vb)

