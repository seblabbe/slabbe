# -*- coding: utf-8 -*-
r"""
Multidimensional Continued Fraction Algorithms (Cython code)

See also the Python code which provides more methods.

EXAMPLES::

    sage: from slabbe.mult_cont_frac_pyx import Brun
    sage: algo = Brun()

Orbit in the cone (with dual coordinates)::

    sage: algo.cone_orbit_list((10,23,15), 6)
    [(((10.0, 8.0, 15.0), (1.0, 1.0, 2.0)), 132),
     (((10.0, 8.0, 5.0), (3.0, 1.0, 2.0)), 213),
     (((2.0, 8.0, 5.0), (3.0, 4.0, 2.0)), 321),
     (((2.0, 3.0, 5.0), (3.0, 4.0, 6.0)), 132),
     (((2.0, 3.0, 2.0), (3.0, 10.0, 6.0)), 123),
     (((2.0, 1.0, 2.0), (3.0, 10.0, 16.0)), 132)]

Orbit in the simplex::

    sage: algo.simplex_orbit_list((10,23,15), 3)
    [(0.30303030303030304,
      0.24242424242424246,
      0.45454545454545453,
      0.25,
      0.25,
      0.5,
      132),
     (0.43478260869565216,
      0.3478260869565218,
      0.21739130434782603,
      0.5,
      0.16666666666666666,
      0.3333333333333333,
      213),
     (0.13333333333333328,
      0.5333333333333334,
      0.3333333333333333,
      0.33333333333333337,
      0.4444444444444445,
      0.22222222222222224,
      321)]

BENCHMARKS:

With slabbe-0.2 or earlier, 68.6 ms on my machine.
With slabbe-0.3.b1, 62.2 ms on my machine.
With slabbe-0.3.b2, 28.6 ms on my machine.
With slabbe-0.3.b2, 13.3 ms on priminfo in Liège::

    sage: from slabbe.mult_cont_frac_pyx import Brun
    sage: %time Brun().lyapunov_exponents(n_iterations=10^6)  # not tested
    (0.3049429393152174, -0.1120652699014143, 1.367495867105725)

With slabbe-0.3.b1, 74ms on my machine.
With slabbe-0.3.b2, 35ms on my machine.
With slabbe-0.3.b2, 17ms on priminfo in Liège::

    sage: from slabbe.mult_cont_frac_pyx import ARP
    sage: %time ARP().lyapunov_exponents(n_iterations=10^6)  # not tested
    (0.443493194984839, -0.17269097306340797, 1.3893881011394358)

With slabbe-0.2 or earlier, 3.71s at liafa, 4.58s on my machine.
With slabbe-0.3.b1, 3.93s on my machine.
With slabbe-0.3.b2, 1.93s on my machine.
With slabbe-0.3.b2, 1.22s on priminfo in Liège::

    sage: %time Brun().lyapunov_exponents(n_iterations=67000000) # not tested
    (0.30456433843239084, -0.1121770192467067, 1.36831961293987303)

With slabbe-0.3.b1, 4.83 s on my machine:
With slabbe-0.3.b2, 2.33 s on my machine:
With slabbe-0.3.b2, 1.56 s on priminfo in Liège::

    sage: %time ARP().lyapunov_exponents(n_iterations=67*10^6)   # not tested
    (0.44296596371477626, -0.17222952278277034, 1.3888098339168744)

With slabbe-0.2 or earlier, 660 ms on my machine.
With slabbe-0.3.b1, 640 ms on my machine (maybe this test could be made much
faster without using list...).
With slabbe-0.3.b2, 215 ms on priminfo in Liège::

    sage: %time L = Brun().simplex_orbit_list(n_iterations=10^6)   # not tested

.. TODO::

    - Ajout les algo de reuteneaour, nogueira, autres?
    - Allow 2d, 1d, 4d, algorithms
    - utilise les vecteurs pour les plus grandes dimensions?

    - Replace method ``_natural_extension_dict`` by ``simplex_orbit_filtered_list``
    - Use ``simplex_orbit_filtered_list`` for ``invariant_measure`` ?

    - Essayer d'utiliser @staticmethod pour call pour que ARP puisse
      apeler Poincare. Without the cython Error: Cannot call a static
      method on an instance variable.
      https://groups.google.com/d/topic/sage-support/DRI_s31D8ks/discussion

    - Trouver meilleur nom pour PointFiber + x + a (cotangent, bundle, fiber?)

    - Create memory inside algorithms (or points?) to perform computations...

    - See the TODO in the file

    - Move more methods into the PairPoint class

    - cone_orbit_list is broken because it does not return distinct
      points... (21 Feb 2017)

    - Code qsort and sort and permutation (21 Feb 2017) so that we can use
      this code in multidimensional case

    - In order for TestSuite(Brun()).run() to work properly, Brun must be an
      instance of SageObject. Is this a bug? But then using SageObject does
      not compile to c. Need to check this...

Question:

    - Comment factoriser le code sans utiliser les yield?
    - Comment faire un appel de fonction rapide (pour factoriser le code)

AUTHORS:

 - Sébastien Labbé, Invariant measures, Lyapounov exponents and natural
   extensions for a dozen of algorithms, October 2013.
 - Sébastien Labbé, Cleaning the code, Fall 2015
 - Sébastien Labbé, Making use of PairPoint to prepare for higher dimension, Fall 2016

"""
#*****************************************************************************
#       Copyright (C) 2013-2018 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
#from __future__ import absolute_import # creates a failing test in mult_cont_frac.py ??

from libc.math cimport log

from cysignals.signals cimport sig_check   # ctrl-c interrupt block support
from cysignals.memory cimport check_allocarray, sig_free

cdef double SQRT3SUR2 = 0.866025403784439

# random
# https://groups.google.com/d/msg/cython-users/jc-3UK2Ffoc/GJNVx-CzdKsJ
from libc.stdlib cimport rand, RAND_MAX
cdef double RAND_SCALE = 1.0 / RAND_MAX
cdef inline double random_double_32bit():
    return rand() * RAND_SCALE
cdef inline double random_double_64bit():
    return rand() * RAND_SCALE + rand() * RAND_SCALE * RAND_SCALE 

# qsort
# http://stackoverflow.com/questions/35095914/sorting-in-cython 
# from libc.stdlib cimport qsort
# ... declaring "const void *" type seems problematic

# http://stackoverflow.com/questions/8353076/how-do-i-pass-a-pointer-to-a-c-fun$
cdef extern from "stdlib.h":
    ctypedef void const_void "const void"
    void qsort(void *base, int nmemb, int size,
                int(*compar)(const_void *, const_void *)) nogil

cdef int cmp_double(const_void * pa, const_void * pb):
    r"""
    cmp of doubles
    """
    cdef double a = (<double *>pa)[0]
    cdef double b = (<double *>pb)[0]
    if a < b:
        return -1
    elif a > b:
        return 1
    else:
        return 0

cdef double* _KEY
cdef int cmp_int_KEY(const_void * pa, const_void * pb):
    r"""
    cmp of integers according to their values in global array _KEY 
    """
    cdef int a = (<int *>pa)[0]
    cdef int b = (<int *>pb)[0]
    if _KEY[a] < _KEY[b]:
        return -1
    elif _KEY[a] > _KEY[b]:
        return 1
    else:
        return 0

cdef class PairPoint:
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_pyx import PairPoint
        sage: PairPoint(3, (.2,.3,.4))
        ((0.2, 0.3, 0.4), (..., ..., ...))
        sage: PairPoint(3, a=(.2,.3,.4))
        ((..., ..., ...), (0.2, 0.3, 0.4))
    """
    cdef double* x
    cdef double* a
    cdef int* perm
    cdef int dim
    cdef double _tmp
    def __cinit__(self, int dim, *args, **kwds):
        self.x = <double*>check_allocarray(dim, sizeof(double))
        self.a = <double*>check_allocarray(dim, sizeof(double))
        self.perm = <int*>check_allocarray(dim, sizeof(int))
        self.dim = dim

    def __init__(self, int dim, x=None, a=None):
        cdef int i
        if x is None:
            for i in range(self.dim):
                self.x[i] = random_double_64bit()
        else:
            for i in range(self.dim):
                self.x[i] = x[i]
        if a is None:
            for i in range(self.dim):
                self.a[i] = random_double_64bit()
        else:
            for i in range(self.dim):
                self.a[i] = a[i]

    def __dealloc__(self):
        sig_free(self.x)
        sig_free(self.a)


    def to_dict(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: PairPoint(3, [1,2,3], [4,5,6]).to_dict()
            {'a': [4.0, 5.0, 6.0], 'x': [1.0, 2.0, 3.0]}
        """
        cdef int i
        return dict(x=[self.x[i] for i in range(self.dim)],
                    a=[self.a[i] for i in range(self.dim)])
    def to_tuple(self):
        r"""
        EXAMPLES::

             sage: from slabbe.mult_cont_frac_pyx import PairPoint
             sage: PairPoint(3, [1,2,3], [4,5,6]).to_tuple()
             ((1.0, 2.0, 3.0), (4.0, 5.0, 6.0))
        """
        cdef int i
        return (tuple(self.x[i] for i in range(self.dim)),
                tuple(self.a[i] for i in range(self.dim)))

    def __repr__(self):
        r"""
        EXAMPLES::

             sage: from slabbe.mult_cont_frac_pyx import PairPoint
             sage: PairPoint(5, range(5), range(5))
             ((0.0, 1.0, 2.0, 3.0, 4.0), (0.0, 1.0, 2.0, 3.0, 4.0))
        """
        return repr(self.to_tuple())

    cdef void sort(self):
        r"""
        # How to sort x and a arrays according to x?
        """
        if self.dim != 3:
            raise NotImplementedError("sort is implemented only for dim=3")
        if self.x[0] > self.x[1]:
            self._tmp = self.x[0]
            self.x[0] = self.x[1]
            self.x[1] = self._tmp
            self._tmp = self.a[0]
            self.a[0] = self.a[1]
            self.a[1] = self._tmp
        if self.x[1] > self.x[2]:
            self._tmp = self.x[1]
            self.x[1] = self.x[2]
            self.x[2] = self._tmp
            self._tmp = self.a[1]
            self.a[1] = self.a[2]
            self.a[2] = self._tmp
        if self.x[0] > self.x[1]:
            self._tmp = self.x[0]
            self.x[0] = self.x[1]
            self.x[1] = self._tmp
            self._tmp = self.a[0]
            self.a[0] = self.a[1]
            self.a[1] = self._tmp

    cpdef void sort_x(self): # nogil:
        r"""
        Sort array x independently of array a.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(4, [.4, .2, .3, .1], [4,3,2,1])
            sage: P.sort_x()
            sage: P
            ((0.1, 0.2, 0.3, 0.4), (4.0, 3.0, 2.0, 1.0))
        """
        qsort(self.x, self.dim, sizeof(double), cmp_double)

    cpdef int permutation(self):
        r"""
        http://stackoverflow.com/questions/17554242/how-to-obtain-the-index-permutation-after-the-sorting

        OUTPUT:

            int (the permutation, works well if self.dim < 10)

            Permutation gets written to self.perm

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(4, [.4, .2, .3, .1], [4,3,2,1])
            sage: P.permutation()
            4231
        """
        cdef int i
        for i in range(self.dim):
            self.perm[i] = i
        global _KEY
        _KEY = self.x
        qsort(self.perm, self.dim, sizeof(int), cmp_int_KEY)
        cdef rep = 0
        for i in range(self.dim):
            rep *= 10
            rep += self.perm[i] + 1
        return rep

    cdef PairPoint copy(self):
        r"""
        """
        cdef int i
        return PairPoint(self.dim, 
                         x=[self.x[i] for i in range(self.dim)],
                         a=[self.a[i] for i in range(self.dim)])

    cdef void copy_inplace(self, PairPoint P):
        r"""
        Copy P into self
        """
        cdef int i
        assert self.dim == P.dim, "dimension inconsistencies"
        for i in range(self.dim):
            self.x[i] = P.x[i]
        for i in range(self.dim):
            self.a[i] = P.a[i]

    cdef double min_x(self):
        r"""
        Return the minimum entry of vector x

        OUTPUT:

            double
        """
        cdef int i
        self._tmp = self.x[0]
        for i in range(1, self.dim):
            if self._tmp > self.x[i]:
                self._tmp = self.x[i]
        return self._tmp

    cdef double max_x(self):
        r"""
        Return the maximum entry of vector x

        OUTPUT:

            double
        """
        cdef int i
        self._tmp = self.x[0]
        for i in range(1, self.dim):
            if self._tmp < self.x[i]:
                self._tmp = self.x[i]
        return self._tmp

    cdef double norm_x(self, int p=1):
        r"""
        Return the p-norm of vector x

        INPUT:

        - ``p`` -- integer, 0 or 1

        OUTPUT:

            double
        """
        cdef int i
        if p == 1: # 1-norm
            self._tmp = 0
            for i in range(self.dim):
                self._tmp += abs(self.x[i])
            return self._tmp
        elif p == 0: # sup norm
            return self.max_x()

    cdef double norm_a(self, int p=1):
        r"""
        Return the p-norm of vector a

        INPUT:

        - ``p`` -- integer, 0 or 1

        OUTPUT:

            double

        .. TODO:: 

            How to deal with the scalar product case?
        """
        cdef int i
        if p == 1: # 1-norm
            self._tmp = 0
            for i in range(self.dim):
                self._tmp += abs(self.a[i])
            return self._tmp
        elif p == 0: # sup norm
            self._tmp = 0
            for i in range(self.dim):
                self._tmp = max(self._tmp, abs(self.a[i]))
            return self._tmp
        elif p == -1: # hypersurface
            return self.dot_product()

    cdef void normalize_x(self, double value):
        r"""
        Normalize vector x by dividing each entry by some value.

        INPUT:

        - ``value`` -- positive real number
        """
        cdef int i
        for i in range(self.dim):
            self.x[i] /= value

    cdef void normalize_a(self, double value):
        r"""
        Normalize vector a by dividing each entry by some value.

        INPUT:

        - ``value`` -- positive real number
        """
        cdef int i
        for i in range(self.dim):
            self.a[i] /= value

    cdef double dot_product(self):
        r"""
        TODO: peut-on utiliser _tmp ici? ou aura-t-on des probleme de
        copie?
        """
        cdef int i
        self._tmp = 0
        for i in range(self.dim):
            self._tmp += self.x[i] * self.a[i]
        return self._tmp

    cdef double dot_product_xx(self):
        cdef int i
        cdef double s = 0
        for i in range(self.dim):
            s += self.x[i] * self.x[i]
        return s

    cdef void gramm_schmidt(self):
        r"""
        Removes x component for vector a.
        """
        cdef int i
        cdef double p,s
        p = self.dot_product()
        s = self.dot_product_xx()
        for i in range(self.dim):
            self.a[i] -= p*self.x[i]/s

    cpdef int number_small_entries(self, double ratio, int p=1):
        r"""
        Returns the number of indices i such that x[i]/||x|| < ratio.
        """
        cdef int i,c=0
        cdef double norm_ratio = self.norm_x(p) * ratio
        for i in range(self.dim):
            if self.x[i] < norm_ratio:
                c += 1
        return c

    cdef act_by_diagonal_matrix(self):
        raise NotImplementedError

    cdef (double, double) projection3to2(self, int p=1):
        assert self.dim == 3, ("Dimension of point is {} but projection"
                               " implemented only for 3".format(self.dim))
        if p == 1:
            return (-SQRT3SUR2 * self.x[0] + SQRT3SUR2 * self.x[1],
                    -.5 * self.x[0] -.5 * self.x[1] + self.x[2])
        elif p == 0:
            return self.x[0], self.x[1]

    cdef int subcone(self):
        r"""
        .. NOTE::

            Calling this method seems slower than copying the same code inside
            _Poincare method.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
        """
        if self.x[0] <= self.x[1] <= self.x[2]:
            return 123
        elif self.x[0] <= self.x[2] <= self.x[1]:
            return 132
        elif self.x[1] <= self.x[2] <= self.x[0]:
            return 231
        elif self.x[1] <= self.x[0] <= self.x[2]:
            return 213
        elif self.x[2] <= self.x[0] <= self.x[1]:
            return 312
        elif self.x[2] <= self.x[1] <= self.x[0]:
            return 321
        else:
            return -1

    cdef void op_elem(self, int i, int j):
        r"""
        Elementary operation: `x_i = x_i - x_j` and `a_j = a_j + a_i`

        INPUT:

        - ``i`` -- integer, betwen 1 and dimension - 1
        - ``j`` -- integer, betwen 1 and dimension - 1

        .. NOTE::

            Calling this method seems slower than copying the same code inside
            _Poincare method.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
        """
        self.x[i] -= self.x[j]
        self.a[j] += self.a[i]

    cdef int _Poincare(self):
        r"""
        EXAMPLES::

        """
        cdef branch = self.permutation() # result written in self.perm
        cdef int i, s,t
        for i in range(self.dim-1, 0, -1):
            t = self.perm[i]
            s = self.perm[i-1]
            self.x[t] -= self.x[s]
            self.a[s] += self.a[t]
        return branch

    cdef int _Poincare3(self):
        r"""
        EXAMPLES::

        """
        cdef int i,j,k
        if self.x[0] <= self.x[1] <= self.x[2]:
            i = 0; j = 1; k = 2;
        elif self.x[0] <= self.x[2] <= self.x[1]:
            i = 0; j = 2; k = 1;
        elif self.x[1] <= self.x[0] <= self.x[2]:
            i = 1; j = 0; k = 2;
        elif self.x[2] <= self.x[0] <= self.x[1]:
            i = 2; j = 0; k = 1;
        elif self.x[1] <= self.x[2] <= self.x[0]:
            i = 1; j = 2; k = 0;
        elif self.x[2] <= self.x[1] <= self.x[0]:
            i = 2; j = 1; k = 0;
        self.x[k] -= self.x[j]
        self.x[j] -= self.x[i]
        self.a[j] += self.a[k]
        self.a[i] += self.a[j]
        return 100*i + 10*j + k + 111

    cdef int _Poincare_slower(self):
        r"""
        .. NOTE::

            This method is slower than _Poincare.

        EXAMPLES::

        """
        cdef int part,i,j,k
        part = self.subcone()
        i = (part // 100 % 10) - 1
        j = (part // 10 % 10) - 1
        k = (part % 10) - 1
        self.op_elem(k, j)
        self.op_elem(j, i)
        return part

    cdef int _Sorted_Poincare(self):
        r"""
        EXAMPLES::

        """
        # Apply the algo
        self.x[2] -= self.x[1]
        self.x[1] -= self.x[0]
        self.a[1] += self.a[2]
        self.a[0] += self.a[1]
        self.sort()
        return 200

    cdef int _Sorted_ArnouxRauzy(self):
        r"""
        EXAMPLES::

        """
        #Arnoux-Rauzy
        self.x[2] -= self.x[0] + self.x[1]
        self.a[0] += self.a[2]
        self.a[1] += self.a[2]
        self.sort()
        return 100

    cdef int _Sorted_ArnouxRauzyMulti(self):
        r"""
        EXAMPLES::

        """
        #Arnoux-Rauzy Multi
        cdef int m
        m = <int>(self.x[2] / (self.x[0] + self.x[1]))
        self.x[2] -= m * (self.x[0] + self.x[1])
        self.a[1] += m * self.a[2];
        self.a[0] += m * self.a[2];
        self.sort()
        return 100

    cdef int _Sorted_FullySubtractive(self):
        r"""
        EXAMPLES::

        """
        # Apply the algo
        self.x[1] -= self.x[0]
        self.x[2] -= self.x[0]
        self.a[0] += self.a[1] + self.a[2]
        self.sort()
        return 100


cdef class MCFAlgorithm(object):
    cdef int dim
    def __cinit__(self, int dim=3):
        self.dim = dim
    #def __init__(self, int dim=3):
    #    self.dim = dim
    def __reduce__(self):
        r"""
        Default pickle support

        TESTS::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: algo = Brun()
            sage: algo.__reduce__()
            (<type 'slabbe.mult_cont_frac_pyx.Brun'>, ())
        """
        return self.__class__, tuple()

    ########################################
    # METHODS IMPLEMENTED IN HERITED CLASSES
    ########################################
    cdef int call(self, PairPoint P) except *:
        r"""
        This method must be implemented in the inherited classes.

        EXAMPLES::

        """
        raise NotImplementedError
    def substitutions(self):
        r"""
        This method must be implemented in the inherited classes.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: Brun().substitutions()
            {123: WordMorphism: 1->1, 2->23, 3->3,
             132: WordMorphism: 1->1, 2->2, 3->32,
             213: WordMorphism: 1->13, 2->2, 3->3,
             231: WordMorphism: 1->1, 2->2, 3->31,
             312: WordMorphism: 1->12, 2->2, 3->3,
             321: WordMorphism: 1->1, 2->21, 3->3}

        """
        raise NotImplementedError

    def dual_substitutions(self):
        r"""
        This method must be implemented in the inherited classes.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: Brun().dual_substitutions()
            {123: WordMorphism: 1->1, 2->2, 3->32,
             132: WordMorphism: 1->1, 2->23, 3->3,
             213: WordMorphism: 1->1, 2->2, 3->31,
             231: WordMorphism: 1->13, 2->2, 3->3,
             312: WordMorphism: 1->1, 2->21, 3->3,
             321: WordMorphism: 1->12, 2->2, 3->3}
        """
        raise NotImplementedError

    def branches(self, int n_iterations=1000):
        r"""
        Returns the branches labels of the algorithm.

        This method is an heuristic and should be implemented in the
        inherited classes.

        EXAMPLES::

            sage: import slabbe.mult_cont_frac_pyx as mcf
            sage: mcf.Brun().branches()
            {123, 132, 213, 231, 312, 321}
            sage: mcf.ARP().branches()
            {1, 2, 3, 123, 132, 213, 231, 312, 321}
        """
        cdef unsigned int i         # loop counter
        cdef PairPoint P
        cdef int branch
        S = set()
        # Loop
        for i in range(n_iterations):

            # Check for Keyboard interupt
            sig_check()

            # random initial values
            P = PairPoint(self.dim)

            # Apply Algo
            branch = self.call(P)

            S.add(branch)
        return S
    ######################
    # TEST METHODS 
    ######################
    def _test_definition(self, int n_iterations=10000):
        r"""
        INPUT:

        - ``n_iterations`` -- integer

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: Brun()._test_definition(10000)
        """
        cdef double s,t             # temporary variables
        cdef unsigned int i         # loop counter
        cdef PairPoint P,R
        cdef int branch

        # Loop
        for i in range(n_iterations):

            sig_check() # Check for Keyboard interupt

            # random initial values
            P = PairPoint(self.dim)
            s = P.dot_product()

            R = P.copy()

            # Apply Algo
            try:
                branch = self.call(R)
            except ValueError:
                continue
            t = R.dot_product()

            if not abs(s - t) < 0.0000001:
                m = 'This algo does not preserve the scalar product\n'
                m += '{} != {}\n'.format(s,t)
                m += 'The problem is on branch {}\n'.format(branch)
                m += 'on the {}-th iteration\n'.format(i)
                m += 'INPUT: {}\n'.format(P)
                m += 'OUTPUT: {}\n'.format(R)
                raise Exception(m)

        return

    def _test_dual_substitution_definition(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: t = Brun()._test_dual_substitution_definition()
        """
        A = self.substitutions()
        B = self.dual_substitutions()
        assert set(A.keys()) == set(B.keys())
        for key in A:
            a = A[key].incidence_matrix()
            b = B[key].incidence_matrix()
            if not a == b.transpose():
                raise ValueError("Transpose of substitution {} do not "
                        "match with dual substitution for algo "
                        " {}".format(key, self.name()))
        return

    def _test_coherence(self, int n_iterations=1000):
        r"""
        Check coherence between substitutions and the algo.

        INPUT:

        - ``n_iterations`` -- integer

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: t = Brun()._test_coherence()
        """
        from sage.modules.free_module_element import vector
        cdef unsigned int i         # loop counter
        cdef PairPoint P,R
        cdef int branch

        A = dict((k,s.incidence_matrix()) for k,s in self.substitutions().iteritems())

        # Loop
        for i in range(n_iterations):

            # Check for Keyboard interupt
            sig_check()

            # random initial values
            P = PairPoint(self.dim)
            
            R = P.copy()  # TODO is the copy needed here?

            # Apply Algo
            try:
                branch = self.call(R)
            except ValueError:
                continue

            M = A[branch]

            # Check the algo
            s = M * vector((R.x[0], R.x[1], R.x[2]))
            t = vector((P.x[0], P.x[1], P.x[2]))
            if not abs(s - t) < 0.0000001:
                m = 'Incoherence between the definition of algo \n'
                m += 'and the associated substitutions.\n'
                m += 'M * v_{{n+1}} = {} !=\nv_n = {}\n'.format(s,t)
                m += 'The problem is on branch {}\n'.format(branch)
                m += 'on the {}-th iteration\n'.format(i)
                m += 'INPUT: {}\n'.format(P)
                m += 'OUTPUT: {}\n'.format(R)
                raise Exception(m)

            # Check the dual coordinates
            Mt = M.transpose()
            s = Mt * vector((P.a[0], P.a[1], P.a[2]))
            t = vector((R.a[0], R.a[1], R.a[2]))
            if not abs(s - t) < 0.0000001:
                m = 'Incoherence between the definition of algo (dual) \n'
                m += 'and the associated substitutions.\n'
                m += '{} != {}\n'.format(s,t)
                m += 'The problem is on branch {}\n'.format(branch)
                m += 'on the {}-th iteration\n'.format(i)
                m += 'INPUT: {}\n'.format(P)
                m += 'OUTPUT: {}\n'.format(R)
                raise Exception(m)
        return

    ######################
    # METHODS FOR THE USER:
    ######################
    def class_name(self):
        r"""
        The name of the class.

        .. NOTE::

            This might not be the same as the name of the algorithm.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Reverse, Brun, ARP
            sage: Reverse().class_name()
            'Reverse'
            sage: Brun().class_name()
            'Brun'
            sage: ARP().class_name()
            'ARP'
        """
        return self.__class__.__name__
    def name(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Reverse, Brun, ARP
            sage: Reverse().name()
            'Reverse'
            sage: Brun().name()
            'Brun'
            sage: ARP().name()
            "Arnoux-Rauzy-Poincar\\'e"
        """
        return self.class_name()

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Reverse, Brun
            sage: Reverse()
            Reverse 3-dimensional continued fraction algorithm
            sage: Brun()
            Brun 3-dimensional continued fraction algorithm
        """
        return "{} {}-dimensional continued fraction algorithm".format(self.name(), self.dim)

    def __call__(self, PairPoint P):
        r"""
        Wrapper for the cdef call method.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: Brun()(P)
            ((0.3, 0.6, 0.20000000000000007), (0.2, 0.6, 0.3))

        ::

            sage: P = PairPoint(3, [3,6,8], [.2,.3,.3])
            sage: Brun()(P)
            ((3.0, 6.0, 2.0), (0.2, 0.6, 0.3))
        """
        # TODO should __call__ return a copy or changes P?
        cdef PairPoint R = P.copy()
        self.call(R)
        return R

    def __richcmp__(self, other, op):
        r"""
        INPUT:

        - ``other`` -- 
        - ``op`` -- int, from 0 to 5

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun, ARP
            sage: Brun() == ARP()
            False
            sage: Brun() == Brun()
            True
        """
        # 0: <
        # 1: <=
        # 2: ==
        # 3: !=
        # 4: >
        # 5: >=
        if not isinstance(other, MCFAlgorithm):
            return NotImplemented
        if op == 2 or op == 3:
            return other.class_name() == self.class_name()
        else:
            return NotImplemented
    ######################
    # DYNAMICS METHODS
    ######################
    def coding_iterator(self, start):
        r"""
        INPUT:

        - ``start`` -- iterable of three real numbers

        OUTPUT:

            iterator

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import ARP
            sage: it = ARP().coding_iterator((1,e,pi))
            sage: [next(it) for _ in range(20)]
            [123, 2, 1, 123, 1, 231, 3, 3, 3, 3, 123, 1, 1, 1, 231, 2, 321, 2, 3, 312]

        ::

            sage: from slabbe.mult_cont_frac_pyx import Poincare
            sage: algo = Poincare(4)
            sage: it = algo.coding_iterator((1,e,pi,sqrt(2)))
            sage: [next(it) for _ in range(10)]
            [1423, 4312, 3241, 3412, 3142, 3214, 4312, 1342, 3412, 1342]
        """
        cdef PairPoint P = PairPoint(self.dim, start)
        cdef int branch

        while True:
            branch = self.call(P)
            yield branch
            P.normalize_x(P.norm_x())

    def simplex_orbit_iterator(self, start=None, int norm_xyz=0, int norm_uvw=1):
        r"""
        INPUT:

        - ``start`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``norm_xyz`` -- integer (default: ``0``), either ``0`` or ``1``, the
          norm used for the orbit of points `(x,y,z)` of the algo
        - ``norm_uvw`` -- integer (default: ``1``), either ``0`` or
          ``1`` or ``'hypersurfac'``, the norm used for the orbit of dual
          coordinates `(u,v,w)`.

        NOTE:

            This iterator is 10x slower because of the yield statement. So
            avoid using this when writing fast code. Just copy paste the
            loop or use simplex_orbit_list or simplex_orbit_filtered_list method.

        OUTPUT:

            iterator

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: it = Brun().simplex_orbit_iterator((.414578,.571324,.65513))
            sage: for _ in range(4): next(it)
            ((0.7256442929056017, 1.0, 0.14668734378391243), 
             (0.25, 0.5, 0.25), 
             123)
            ((1.0, 0.37808566783572695, 0.20214772612150184),
             (0.5, 0.3333333333333333, 0.16666666666666666),
             312)
            ((1.0, 0.6079385025908344, 0.32504111204194974),
             (0.3333333333333333, 0.5555555555555555, 0.1111111111111111),
             321)
            ((0.6449032192209051, 1.0, 0.534661171576946),
             (0.25, 0.6666666666666666, 0.08333333333333333),
             321)

        ::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: it = Brun().simplex_orbit_iterator((.414578,.571324,.65513), norm_xyz=1)
            sage: for _ in range(4): next(it)
            ((0.3875618393056797, 0.5340934161472103, 0.07834474454711005),
             (0.25, 0.5, 0.25),
             123)
             ((0.6328179140018012, 0.23925938363378257, 0.12792270236441622),
             (0.5, 0.3333333333333333, 0.16666666666666666),
             312)
            ((0.5173360300491189, 0.3145084914443481, 0.16815547850653312),
             (0.3333333333333333, 0.5555555555555555, 0.1111111111111111),
             321)
            ((0.2958862889959549, 0.45880727553726447, 0.24530643546678058),
             (0.25, 0.6666666666666666, 0.08333333333333333),
             321)

        """
        cdef double s           # temporary variables
        cdef PairPoint P
        cdef int branch, i
        if start is None:
            P = PairPoint(self.dim)
        else:
            P = PairPoint(self.dim, start, [1./self.dim for i in range(self.dim)])

        P.normalize_x(P.norm_x())

        # Loop
        while True:

            # Apply Algo
            branch = self.call(P)

            P.normalize_x(P.norm_x(p=norm_xyz))
            P.normalize_a(P.norm_a(p=norm_uvw))

            yield (P.x[0], P.x[1], P.x[2]), (P.a[0], P.a[1], P.a[2]), branch

    def simplex_orbit_list(self, start=None, int n_iterations=100, 
                                 int norm_xyz=1, int norm_uvw=1):
        r"""
        INPUT:

        - ``start`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``n_iterations`` - integer, number of iterations
        - ``norm_xyz`` -- integer (default: ``1``), either ``0`` or ``1``, the
          norm used for the orbit of points `(x,y,z)` of the algo
        - ``norm_uvw`` -- integer (default: ``1``), either ``0`` or
          ``1`` or ``'hypersurfac'``, the norm used for the orbit of dual
          coordinates `(u,v,w)`.

        OUTPUT:

            list

        .. NOTE::

            It could be 10 times faster because 10^6 iterations can be done in
            about 60ms on this machine. But for drawing images, it does not
            matter to be 10 times slower::

                sage: %time L = Brun().simplex_orbit_list(10^6)   # not tested
                CPU times: user 376 ms, sys: 267 ms, total: 643 ms
                Wall time: 660 ms

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: L = Brun().simplex_orbit_list(n_iterations=10^5)
            sage: L[-1]    # random
            (0.7307002153148079,
             1.0,
             0.31588474491578816,
             0.29055326655584235,
             0.4690741038784866,
             0.24037262956567113,
             321)

        """
        cdef double s           # temporary variables
        cdef PairPoint P
        cdef int i
        cdef int branch
        if start is None:
            P = PairPoint(self.dim)
        else:
            P = PairPoint(self.dim, start, [1./self.dim for i in range(self.dim)])

        P.normalize_x(P.norm_x(p=norm_xyz))

        L = []

        # Loop
        for i in range(n_iterations):

            # Check for Keyboard interupt
            sig_check()

            # Apply Algo
            branch = self.call(P)

            P.normalize_x(P.norm_x(p=norm_xyz))
            P.normalize_a(P.norm_a(p=norm_uvw))

            L.append( (P.x[0], P.x[1], P.x[2], P.a[0], P.a[1], P.a[2], branch))

        return L

    def simplex_orbit_filtered_list(self, start=None, int n_iterations=100,
            int norm_xyz=1, int norm_uvw=1,
            double xmin=-float('inf'), double xmax=float('inf'),
            double ymin=-float('inf'), double ymax=float('inf'),
            double umin=-float('inf'), double umax=float('inf'),
            double vmin=-float('inf'), double vmax=float('inf'),
            int ndivs=0):
        r"""
        Return a list of the orbit filtered to fit into a rectangle.

        INPUT:

        - ``start`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``n_iterations`` - integer, number of iterations
        - ``norm_xyz`` -- integer (default: ``1``), either ``0`` or ``1``, the
          norm used for the orbit of points `(x,y,z)` of the algo
        - ``norm_uvw`` -- integer (default: ``1``), either ``0`` or
          ``1`` or ``'hypersurfac'``, the norm used for the orbit of dual
          coordinates `(u,v,w)`.
        - ``xmin`` - double
        - ``ymin`` - double
        - ``umin`` - double
        - ``vmin`` - double
        - ``xmax`` - double
        - ``ymax`` - double
        - ``umax`` - double
        - ``vmax`` - double
        - ``ndvis`` - integer, number of divisions

        OUTPUT:

            list

        BENCHMARK::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: %time D = Brun().simplex_orbit_filtered_list(10^6) # not tested
            CPU times: user 366 ms, sys: 203 ms, total: 568 ms
            Wall time: 570 ms

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: start=(.414578,.571324,.65513)
            sage: D = Brun().simplex_orbit_filtered_list(start, 3)
            sage: D      # random
            [(0.3049590483124023,
              -0.36889249928767137,
              -0.21650635094610976,
              -0.125,
              312,
              312),
             (0.08651831333735083,
              -0.31784823591841554,
              -0.34641016151377557,
              -0.2,
              312,
              312),
             (-0.41045591033143647,
              -0.20171750067080554,
              -0.4330127018922195,
              -0.25000000000000006,
              312,
              231)]

        ::

            sage: Brun().simplex_orbit_filtered_list(n_iterations=3, norm_xyz=1,ndivs=1000)
            Traceback (most recent call last):
            ...
            ValueError: when ndivs is specified, you must provide a value
            for xmin, xmax, ymin, ymax, umin, umax, vmin and vmax

        ::

            sage: Brun().simplex_orbit_filtered_list(n_iterations=7,  # random
            ....:       norm_xyz=1, ndivs=100,
            ....:       xmin=-.866, xmax=.866, ymin=-.5, ymax=1.,
            ....:       umin=-.866, umax=.866, vmin=-.5, vmax=1.)
            [(30, 47, 50, 50, 132, 213),
             (15, 83, 33, 66, 213, 231),
             (18, 80, 38, 44, 231, 231),
             (22, 75, 41, 33, 231, 231),
             (30, 68, 43, 26, 231, 231),
             (44, 53, 44, 22, 231, 213),
             (41, 78, 24, 56, 213, 321)]

        """
        cdef double s,x,y,u,v           # temporary variables
        s = x = y = u = v = 0           # initialize to avoid a warning
        cdef PairPoint P
        cdef int branch
        cdef int previous_branch
        cdef int xa,ya,ua,va
        cdef int i
        cdef double xlen = xmax - xmin
        cdef double ylen = ymax - ymin
        cdef double ulen = umax - umin
        cdef double vlen = vmax - vmin
        if start is None:
            P = PairPoint(self.dim)
        else:
            P = PairPoint(self.dim, start, [1./self.dim for i in range(self.dim)])

        if ndivs and float('inf') in [-xmin, -ymin, -umin, -vmin, xmax, ymax, umax, vmax]:
            raise ValueError("when ndivs is specified, you must provide a"
                    " value for xmin, xmax, ymin, ymax, umin, umax, vmin"
                    " and vmax")

        P.normalize_x(P.norm_x(p=norm_xyz))

        L = []

        # Apply Algo once
        branch = self.call(P)

        # Loop
        for i in range(n_iterations):

            # Check for Keyboard interupt
            sig_check()

            P.normalize_x(P.norm_x(p=norm_xyz))
            P.normalize_a(P.norm_a(p=norm_uvw))

            # Projection
            if norm_xyz == 1:
                x = -SQRT3SUR2 * P.x[0] + SQRT3SUR2 * P.x[1]
                y = -.5 * P.x[0] -.5 * P.x[1] + P.x[2]
            elif norm_xyz == 0:
                x = P.x[0]
                y = P.x[1]

            if norm_uvw == 1:
                u = -SQRT3SUR2 * P.a[0] + SQRT3SUR2 * P.a[1]
                v = -.5 * P.a[0] -.5 * P.a[1] + P.a[2]
            elif norm_uvw == 0:
                u = P.a[0]
                v = P.a[1]

            # Apply Algo
            previous_branch = branch
            branch = self.call(P)

            # filter
            if not (xmin < x < xmax and ymin < y < ymax and
                    umin < u < umax and vmin < v < vmax):
                continue

            # ndivs
            if ndivs:
                xa = int( (x-xmin) / xlen * ndivs )
                ya = int( (ymax-y) / ylen * ndivs )
                ua = int( (u-umin) / ulen * ndivs )
                va = int( (vmax-v) / vlen * ndivs )
                L.append( (xa,ya,ua,va, previous_branch, branch))
            else:
                L.append( (x,y,u,v, previous_branch, branch))

        return L

    def cone_orbit_iterator(self, start=None):
        r"""
        INPUT:

        - ``start`` - initial vector (default: ``None``), if None, then
          initial point is random

        NOTE:

            This iterator is 10x slower because of the yield statement. So
            avoid using this when writing fast code. Just copy paste the
            loop or use simplex_orbit_list or simplex_orbit_filtered_list method.

        OUTPUT:

            iterator of tuples (PairPoint, integer)

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: it = Brun().cone_orbit_iterator((13,17,29))
            sage: for _ in range(10): next(it)
            (((13.0, 17.0, 12.0), (1.0, 2.0, 1.0)), 123)
            (((13.0, 4.0, 12.0), (3.0, 2.0, 1.0)), 312)
            (((1.0, 4.0, 12.0), (3.0, 2.0, 4.0)), 231)
            (((1.0, 4.0, 8.0), (3.0, 6.0, 4.0)), 123)
            (((1.0, 4.0, 4.0), (3.0, 10.0, 4.0)), 123)
            (((1.0, 4.0, 0.0), (3.0, 14.0, 4.0)), 123)
            (((1.0, 3.0, 0.0), (17.0, 14.0, 4.0)), 312)
            (((1.0, 2.0, 0.0), (31.0, 14.0, 4.0)), 312)
            (((1.0, 1.0, 0.0), (45.0, 14.0, 4.0)), 312)
            (((1.0, 0.0, 0.0), (59.0, 14.0, 4.0)), 312)
        """
        cdef int branch,i
        cdef PairPoint P
        if start is None:
            P = PairPoint(self.dim)
        else:
            P = PairPoint(self.dim, start, [1 for i in range(self.dim)])

        # Loop
        while True:
            sig_check() # Check for Keyboard interupt
            branch = self.call(P)
            yield P.copy(), branch

    def cone_orbit_list(self, start=None, int n_iterations=100):
        r"""
        INPUT:

        - ``start`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``n_iterations`` - integer, number of iterations

        OUTPUT:

            list of tuples (PairPoint, integer)

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: L = Brun().cone_orbit_list((10, 21, 37), 20)
            sage: L[-1]
            (((1.0, 0.0, 0.0), (68.0, 55.0, 658.0)), 231)

        .. TODO::

            Check for a fixed point stop then.

        """
        cdef PairPoint P
        cdef int branch, i
        if start is None:
            P = PairPoint(self.dim)
        else:
            P = PairPoint(self.dim, start, [1 for i in range(self.dim)])

        # Loop
        L = []
        for i in range(n_iterations):
            sig_check() # Check for Keyboard interupt
            branch = self.call(P)
            L.append( (P.copy(), branch) )
        return L

    def image(self, start, int n_iterations=1):
        r"""
        Return the image of a vector in R^3 after n iterations.

        INPUT:

        - ``start`` - initial vector
        - ``n_iterations`` - integer, number of iterations (default: 1)

        OUTPUT:

            tuple of three floats

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: Brun().image((10, 21, 37))
            (10.0, 21.0, 16.0)
            sage: Brun().image((10, 21, 37), 2)
            (10.0, 5.0, 16.0)
            sage: Brun().image((10, 21, 37), 3)
            (10.0, 5.0, 6.0)
            sage: Brun().image((10, 21, 37), 10)
            (1.0, 1.0, 0.0)
        """
        cdef PairPoint P = PairPoint(self.dim, start)
        cdef int i, branch

        for i in range(n_iterations):
            sig_check() # Check for Keyboard interupt
            branch = self.call(P)
        return (P.x[0], P.x[1], P.x[2])

    def _invariant_measure_dict(self, int n_iterations, int ndivs, v=None,
            int norm=1, verbose=False):
        r"""
        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``ndvis`` - integer (less than or equal to 256), number of divisions per dimension
        - ``v`` - initial vector (default: ``None``)
        - ``norm`` -- integer (default: ``1``), either ``0`` or ``1``, the p-norm
          used for the orbit of the algo
        - ``verbose`` -- bool (default: ``False``)

        OUTPUT:

            dict

        .. NOTE::

            This method should be erased and replaced by code in the spirit of
            simplex_orbit_list. Or otherwise read [1] and change cpdef int
            C[NDIVS][NDIVS][NDIVS]

            [1] https://groups.google.com/forum/?fromgroups=#!topic/sage-devel/NCBmj2KjwEM

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: D = Brun()._invariant_measure_dict(4, 10, verbose=True) # random
            0.0799404500357 0.199341464229 0.720718085735
            0 1
            0.0998433745026 0.248971884172 0.651184741325
            0 2
            0.132942259282 0.331508073966 0.535549666752
            1 3
            0.198868907918 0.495904379777 0.305226712305
            1 4

        ::

            sage: D = Brun()._invariant_measure_dict(100000, 5)
            sage: sorted(D.iteritems())
            [((0, 0), ...),
             ((0, 1), ...),
             ((0, 2), ...),
             ((0, 3), ...),
             ((0, 4), ...),
             ((1, 1), ...),
             ((1, 2), ...),
             ((1, 3), ...),
             ((2, 2), ...)]

        It is 1000 times faster using C counter instead of a python dict counter::

            sage: D = Brun()._invariant_measure_dict(1000000, 10) # 0.05s

        """
        # 146 works, 147 causes segmentation error!!!
        DEF NDIVS = 100
        assert ndivs <= NDIVS, "ndivs(=%s) must be less or equal to %s" % (ndivs, NDIVS)
        cdef double s
        cdef int i,j
        cdef int X,Y
        cdef int branch

        # initialization of the counter
        # change this to something else
        # see https://groups.google.com/forum/?fromgroups=#!topic/sage-devel/NCBmj2KjwEM
        cpdef int C[NDIVS][NDIVS]
        for j from 0 <= j <= ndivs:
            for i from 0 <= i <= ndivs:
                C[i][j] = 0

        cdef PairPoint P = PairPoint(self.dim)
        P.sort()
        P.normalize_x(P.norm_x(p=norm))

        for i in range(n_iterations):

            sig_check()            # Check for Keyboard interupt
            branch = self.call(P)  # Apply Algo

            P.normalize_x(P.norm_x(p=norm))

            # Increase by one the counter for that part
            X = int(P.x[0]*ndivs)
            Y = int(P.x[1]*ndivs)
            C[X][Y] += 1

            if verbose:
                print(P.x[0],P.x[1],P.x[2])
                print(X,Y)

        # Translate the counter into a python dict
        D = {}
        for j from 0 <= j <= ndivs:
            for i from 0 <= i <= ndivs:
                c = C[i][j]
                if c > 0:
                    D[(i,j)] = c
        return D

    def _natural_extension_dict(self, int n_iterations, int norm_xyz=0,
            int norm_uvw=1, verbose=False):
        r"""
        INPUT:

        - ``n_iterations`` -- integer
        - ``norm_xyz`` -- integer (default: ``0``), either ``0`` or ``1``, the
          norm used for the orbit of points `(x,y,z)` of the algo
        - ``norm_uvw`` -- integer (default: ``1``), either ``0`` or
          ``1`` or ``'hypersurfac'``, the norm used for the orbit of dual
          coordinates `(u,v,w)`.
        - ``verbose`` -- bool (default: ``False``)

        OUTPUT:

            dict, dict, dict, dict

        .. NOTE::

            This method should be erased and replaced by simplex_orbit_list

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import ARP
            sage: t = ARP()._natural_extension_dict(10000)
            sage: map(type, t)
            [<type 'collections.defaultdict'>,
             <type 'collections.defaultdict'>,
             <type 'collections.defaultdict'>,
             <type 'collections.defaultdict'>]

        ::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: t = Brun()._natural_extension_dict(3000, norm_xyz=1, norm_uvw=1)
            sage: map(len, t)
            [6, 6, 6, 6]

        """
        cdef double x,y,z           # vector (x,y,z)
        cdef double u,v,w           # vector (u,v,w)
        cdef double p,s,t           # temporary variables
        cdef unsigned int i         # loop counter
        cdef double x_new,y_new,z_new
        cdef double u_new,v_new,w_new
        cdef int branch

        cdef PairPoint P = PairPoint(self.dim)
        P.sort()
        P.normalize_x(P.norm_x(1))
        P.normalize_a(P.norm_a(1))
        cdef PairPoint R = P.copy()

        import collections
        domain_right = collections.defaultdict(list)
        image_right = collections.defaultdict(list)
        domain_left = collections.defaultdict(list)
        image_left = collections.defaultdict(list)

        # Loop
        for i in range(n_iterations):

            # Check for Keyboard interupt
            sig_check()

            # Apply Algo
            branch = self.call(R)

            if verbose:
                print("x=%f, y=%f, z=%f" % (R.x[0],R.x[1],R.x[2]))
                #print("u=%f, v=%f, w=%f" % (u,v,w))
                #s = x*u + y*v + z*w
                #print("scal prod <(x,y,z),(u,v,w)> = %f (after algo)" % s)

            R.normalize_x(R.norm_x(p=norm_xyz))
            R.normalize_a(R.norm_a(p=norm_uvw))

            # Projection
            if norm_xyz == 1:
                s = -SQRT3SUR2 * P.x[0] + SQRT3SUR2 * P.x[1]
                t = -.5 * P.x[0] -.5 * P.x[1] + P.x[2]
                domain_left[branch].append((s,t))
                s = -SQRT3SUR2 * R.x[0] + SQRT3SUR2 * R.x[1]
                t = -.5 * R.x[0] -.5 * R.x[1] + R.x[2]
                image_left[branch].append((s,t))
            elif norm_xyz == 0:
                domain_left[branch].append((P.x[0],P.x[1]))
                image_left[branch].append((R.x[0],R.x[1]))

            if norm_uvw == 1:
                s = -SQRT3SUR2 * P.a[0] + SQRT3SUR2 * P.a[1]
                t = -.5 * P.a[0] -.5 * P.a[1] + P.a[2]
                domain_right[branch].append((s,t))
                s = -SQRT3SUR2 * R.a[0] + SQRT3SUR2 * R.a[1]
                t = -.5 * R.a[0] -.5 * R.a[1] + R.a[2]
                image_right[branch].append((s,t))
            elif norm_uvw == 0:
                domain_right[branch].append((P.a[0],P.a[1]))
                image_right[branch].append((R.a[0],R.a[1]))

            P.copy_inplace(R)

        return domain_left, image_left, domain_right, image_right

    def lyapunov_exponents(self, start=None, int n_iterations=1000, verbose=False):
        r"""
        Return the lyapunov exponents (theta1, theta2, 1-theta2/theta1)

        See also the module ``slabbe.lyapunov`` for parallel computations.

        INPUT:

        - ``start`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``n_iterations`` -- integer
        - ``verbose`` -- bool (default: ``False``)

        OUTPUT:

            tuple of the first two liapounov exponents and the uniform
            approximation exponent:

            (theta1, theta2, 1-theta2/theta1)

        .. NOTE:: 
        
            the code of this method was translated from C to cython. The C
            version is from Vincent Delecroix.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: Brun().lyapunov_exponents(n_iterations=1000000)  # tol 0.01
            (0.3049429393152174, -0.1120652699014143, 1.367495867105725)

        ::

            sage: start = (0.2134134, 0.31618415, 0.414514985)
            sage: Brun().lyapunov_exponents(start=start, n_iterations=10^6)  # tol 0.01
            (0.3046809303742965, -0.1121152799778245, 1.3679760326322108)
        """
        cdef double theta1=0, theta2=0    # values of Lyapunov exponents
        cdef double theta1c=0, theta2c=0  # compensation (for Kahan summation algorithm)
        cdef double p,s,t           # temporary variables
        cdef unsigned int i         # loop counter
        cdef double critical_value=0.0001
        cdef int branch
        cdef PairPoint P

        # initial values
        P = PairPoint(self.dim, start)
        for i in range(self.dim):
            P.a[i] - .5
        P.sort()
        P.normalize_x(P.norm_x(p=1))
        P.gramm_schmidt()
        P.normalize_a(P.norm_a(p=1))

        if verbose:
            print("P = {}\nscal prod <x,a> = {}".format(P, P.dot_product()))

        # Loop
        for i in range(n_iterations):

            # Check for Keyboard interupt
            sig_check()

            # Apply Algo
            branch = self.call(P)

            if verbose:
                print("P = {}\nscal prod <x,a> = {} (after algo)".format(P, P.dot_product()))

            # Save some computations
            #if i % step == 0:
            if P.x[0] < critical_value:

                # Sum the first lyapunov exponent
                s = P.norm_x()
                p = -log(s) - theta1c
                t = theta1 + p
                theta1c = (t-theta1) - p   # mathematically 0 but not for a computer!!
                theta1 = t
                P.normalize_x(s)

                # Sum the second lyapunov exponent
                s = P.norm_a()
                p = log(s) - theta2c
                t = theta2 + p
                theta2c = (t-theta2) - p   # mathematically 0 but not for a computer!!
                theta2 = t

                # the following gramm shimdts seems to be useless, but it is not!!!
                P.gramm_schmidt()
                P.normalize_a(s)

        return theta1/n_iterations, theta2/n_iterations, 1-theta2/theta1

    def nsmall_entries_list(self, double ratio, start=None, int n_iterations=1000, int p=1):
        r"""
        INPUT:

        - ``ratio`` - real number, 0 < ratio < 1
        - ``start`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``n_iterations`` -- integer
        - ``p`` -- integer, p-norm

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Poincare
            sage: algo = Poincare(4)
            sage: algo.nsmall_entries_list(.1, (1,e,pi,sqrt(2)), n_iterations=20)
            [0, 1, 1, 1, 1, 0, 0, 1, 0, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3]

        ::

            sage: from slabbe.finite_word import run_length_encoding
            sage: L = algo.nsmall_entries_list(.01, (1,e,pi,sqrt(2)), n_iterations=1000)
            sage: run_length_encoding(L)
            [(0, 1), (1, 1), (0, 7), (1, 1), (0, 3), (1, 2), (2, 1), (3, 984)]

        """
        cdef unsigned int i         # loop counter
        cdef int branch
        cdef PairPoint P

        # initial values
        P = PairPoint(self.dim, start, [0 for i in range(self.dim)])
        P.normalize_x(P.norm_x(p=p))

        L = []
        for i in range(n_iterations):
            sig_check() # Check for Keyboard interupt

            branch = self.call(P) # Apply Algo
            P.normalize_x(P.norm_x(p=p))
            L.append(P.number_small_entries(ratio, p=p))
        return L

    def return_time_to_nsmall_entries(self, double ratio, int n, start=None, int p=1):
        r"""
        INPUT:

        - ``ratio`` - real number, 0 < ratio < 1
        - ``n`` - integer, number of small entries
        - ``start`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``p`` -- integer, p-norm

        OUTPUT:

            a tuple (integer, PairPoint)

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Poincare
            sage: algo = Poincare(4)
            sage: algo.return_time_to_nsmall_entries(.05, 0, (1,e,pi,sqrt(2)))
            (3,
             ((0.31830988618379064, 0.41509782135371204, 
               0.13474402056773493, 0.1318482718947624), 
              (0.0, 0.0, 0.0, 0.0)))

        ::

            sage: algo = Poincare(6)
            sage: start = (1,e,pi,sqrt(2),sqrt(3),sqrt(5))
            sage: algo.return_time_to_nsmall_entries(.05, 0, start)
            (5,
             ((0.3183098861837907, 0.153493436015088, 0.134744020567735,
             0.1318482718947624, 0.1011707373432389, 0.16043364799538504),
             (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)))

        """
        cdef unsigned int i=0         # loop counter
        cdef int branch
        cdef PairPoint P

        # initial values
        P = PairPoint(self.dim, start, [0 for i in range(self.dim)])
        P.normalize_x(P.norm_x(p=p))

        while True:
            sig_check() # Check for Keyboard interupt

            branch = self.call(P) # Apply Algo
            P.normalize_x(P.norm_x(p=p))
            if P.number_small_entries(ratio, p=p) == n:
                return i, P

            if P.min_x() == 0:
                raise ValueError("precision problem at "
                        "i={},\nP={}".format(i, P))

            i += 1

cdef class Brun(MCFAlgorithm):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_pyx import Brun
        sage: algo = Brun()
        sage: TestSuite(algo).run()
        sage: algo._test_dual_substitution_definition()
        sage: algo._test_coherence()
        sage: algo._test_definition()
    """
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: Brun()(P)
            ((0.3, 0.6, 0.20000000000000007), (0.2, 0.6, 0.3))
        """
        if P.x[0] <= P.x[1] <= P.x[2]:
            P.x[2] -= P.x[1]
            P.a[1] += P.a[2]
            return 123
        elif P.x[0] <= P.x[2] <= P.x[1]:
            P.x[1] -= P.x[2]
            P.a[2] += P.a[1]
            return 132
        elif P.x[1] <= P.x[2] <= P.x[0]:
            P.x[0] -= P.x[2]
            P.a[2] += P.a[0]
            return 231
        elif P.x[1] <= P.x[0] <= P.x[2]:
            P.x[2] -= P.x[0]
            P.a[0] += P.a[2]
            return 213
        elif P.x[2] <= P.x[0] <= P.x[1]:
            P.x[1] -= P.x[0]
            P.a[0] += P.a[1]
            return 312
        elif P.x[2] <= P.x[1] <= P.x[0]:
            P.x[0] -= P.x[1]
            P.a[1] += P.a[0]
            return 321
        else:
            raise ValueError('limit case: reach set of measure zero: {}'.format(P))

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: Brun().substitutions()
            {123: WordMorphism: 1->1, 2->23, 3->3,
             132: WordMorphism: 1->1, 2->2, 3->32,
             213: WordMorphism: 1->13, 2->2, 3->3,
             231: WordMorphism: 1->1, 2->2, 3->31,
             312: WordMorphism: 1->12, 2->2, 3->3,
             321: WordMorphism: 1->1, 2->21, 3->3}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {312: WordMorphism({1: [1, 2], 2: [2], 3: [3]}),
                321: WordMorphism({1: [1], 2: [2, 1], 3: [3]}),
                213: WordMorphism({1: [1, 3], 2: [2], 3: [3]}),
                231: WordMorphism({1: [1], 2: [2], 3: [3, 1]}),
                123: WordMorphism({1: [1], 2: [2, 3], 3: [3]}),
                132: WordMorphism({1: [1], 2: [2], 3: [3, 2]})}

    def dual_substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Brun
            sage: Brun().dual_substitutions()
            {123: WordMorphism: 1->1, 2->2, 3->32,
             132: WordMorphism: 1->1, 2->23, 3->3,
             213: WordMorphism: 1->1, 2->2, 3->31,
             231: WordMorphism: 1->13, 2->2, 3->3,
             312: WordMorphism: 1->1, 2->21, 3->3,
             321: WordMorphism: 1->12, 2->2, 3->3}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {321: WordMorphism({1: [1, 2], 2: [2], 3: [3]}),
                312: WordMorphism({1: [1], 2: [2, 1], 3: [3]}),
                231: WordMorphism({1: [1, 3], 2: [2], 3: [3]}),
                213: WordMorphism({1: [1], 2: [2], 3: [3, 1]}),
                132: WordMorphism({1: [1], 2: [2, 3], 3: [3]}),
                123: WordMorphism({1: [1], 2: [2], 3: [3, 2]})}

cdef class Reverse(MCFAlgorithm):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_pyx import Reverse
        sage: algo = Reverse()
        sage: TestSuite(algo).run()
        sage: algo._test_dual_substitution_definition()
        sage: algo._test_coherence()
        sage: algo._test_definition()
    """
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Reverse
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: Reverse()(P)
            ((0.55, 0.25, 0.04999999999999993), (0.6, 0.5, 0.5))
        """
        cdef PairPoint R = PairPoint(self.dim)
        if P.x[0] + P.x[1] < P.x[2]:
            P.x[2] -= P.x[0] + P.x[1]
            P.a[1] += P.a[2]
            P.a[0] += P.a[2]
            return 3
        elif P.x[0] + P.x[2] < P.x[1]:
            P.x[1] -= P.x[0] + P.x[2]
            P.a[2] += P.a[1]
            P.a[0] += P.a[1]
            return 2
        elif P.x[1] + P.x[2] < P.x[0]:
            P.x[0] -= P.x[1] + P.x[2]
            P.a[1] += P.a[0]
            P.a[2] += P.a[0]
            return 1
        else:
            # R.x[0] = 0.629960524947437 * (-P.x[0] + P.x[1] + P.x[2])
            # R.x[1] = 0.629960524947437 * ( P.x[0] - P.x[1] + P.x[2])
            # R.x[2] = 0.629960524947437 * ( P.x[0] + P.x[1] - P.x[2])
            # # 0.793700525984100 = 1/2*4^(1/3)
            # # 0.629960524947437 = 1/4*4^(2/3)
            # R.a[0] = 0.793700525984100 * (P.a[1] + P.a[2])
            # R.a[1] = 0.793700525984100 * (P.a[0] + P.a[2])
            # R.a[2] = 0.793700525984100 * (P.a[0] + P.a[1])
            R.x[0] = 0.5 * (-P.x[0] + P.x[1] + P.x[2])
            R.x[1] = 0.5 * ( P.x[0] - P.x[1] + P.x[2])
            R.x[2] = 0.5 * ( P.x[0] + P.x[1] - P.x[2])
            R.a[0] = P.a[1] + P.a[2]
            R.a[1] = P.a[0] + P.a[2]
            R.a[2] = P.a[0] + P.a[1]
            P.copy_inplace(R)
            return 4

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Reverse
            sage: Reverse().substitutions()
            {1: WordMorphism: 1->1, 2->21, 3->31,
             2: WordMorphism: 1->12, 2->2, 3->32,
             3: WordMorphism: 1->13, 2->23, 3->3,
             4: WordMorphism: 1->23, 2->31, 3->12}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {1:  WordMorphism({1: [1], 2: [2, 1], 3: [3, 1]}),
                2:  WordMorphism({1: [1, 2], 2: [2], 3: [3, 2]}),
                3:  WordMorphism({1: [1, 3], 2: [2, 3], 3: [3]}),
                4:  WordMorphism({1: [2,3], 2: [3,1], 3: [1,2]})}

    def dual_substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Reverse
            sage: Reverse().dual_substitutions()
            {1: WordMorphism: 1->123, 2->2, 3->3,
             2: WordMorphism: 1->1, 2->231, 3->3,
             3: WordMorphism: 1->1, 2->2, 3->312,
             4: WordMorphism: 1->23, 2->13, 3->12}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {1: WordMorphism({1: [1,2,3], 2: [2], 3: [3]}),
                2: WordMorphism({1: [1], 2: [2,3,1], 3: [3]}),
                3: WordMorphism({1: [1], 2: [2], 3: [3,1,2]}),
                4:  WordMorphism({1: [2,3], 2: [1,3], 3: [1,2]})}

cdef class ARP(MCFAlgorithm):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_pyx import ARP
        sage: algo = ARP()
        sage: TestSuite(algo).run()
        sage: algo._test_dual_substitution_definition()
        sage: algo._test_coherence()
        sage: algo._test_definition()
    """
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import ARP
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: ARP()(P)
            ((0.3, 0.3, 0.20000000000000007), (0.8, 0.6, 0.3))
        """
        if P.x[0] + P.x[1] < P.x[2]:
            P.x[2] -= P.x[0] + P.x[1]
            P.a[1] += P.a[2]
            P.a[0] += P.a[2]
            return 3
        elif P.x[0] + P.x[2] < P.x[1]:
            P.x[1] -= P.x[0] + P.x[2]
            P.a[2] += P.a[1]
            P.a[0] += P.a[1]
            return 2
        elif P.x[1] + P.x[2] < P.x[0]:
            P.x[0] -= P.x[1] + P.x[2]
            P.a[1] += P.a[0]
            P.a[2] += P.a[0]
            return 1
        else:
            return P._Poincare()
    def name(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import ARP
            sage: ARP().name()
            "Arnoux-Rauzy-Poincar\\'e"
        """
        return r"Arnoux-Rauzy-Poincar\'e"

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import ARP
            sage: ARP().substitutions()
            {1: WordMorphism: 1->1, 2->21, 3->31,
             2: WordMorphism: 1->12, 2->2, 3->32,
             3: WordMorphism: 1->13, 2->23, 3->3,
             123: WordMorphism: 1->123, 2->23, 3->3,
             132: WordMorphism: 1->132, 2->2, 3->32,
             213: WordMorphism: 1->13, 2->213, 3->3,
             231: WordMorphism: 1->1, 2->231, 3->31,
             312: WordMorphism: 1->12, 2->2, 3->312,
             321: WordMorphism: 1->1, 2->21, 3->321}

        """
        from sage.combinat.words.morphism import WordMorphism
        return {1:  WordMorphism({1: [1], 2: [2, 1], 3: [3, 1]}),
                2:  WordMorphism({1: [1, 2], 2: [2], 3: [3, 2]}),
                3:  WordMorphism({1: [1, 3], 2: [2, 3], 3: [3]}),
                312: WordMorphism({1: [1, 2], 2: [2], 3: [3, 1, 2]}),
                321: WordMorphism({1: [1], 2: [2, 1], 3: [3, 2, 1]}),
                213: WordMorphism({1: [1, 3], 2: [2, 1, 3], 3: [3]}),
                231: WordMorphism({1: [1], 2: [2, 3, 1], 3: [3, 1]}),
                123: WordMorphism({1: [1, 2, 3], 2: [2, 3], 3: [3]}),
                132: WordMorphism({1: [1, 3, 2], 2: [2], 3: [3, 2]})}

    def dual_substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import ARP
            sage: ARP().dual_substitutions()
            {1: WordMorphism: 1->123, 2->2, 3->3,
             2: WordMorphism: 1->1, 2->231, 3->3,
             3: WordMorphism: 1->1, 2->2, 3->312,
             123: WordMorphism: 1->1, 2->21, 3->321,
             132: WordMorphism: 1->1, 2->231, 3->31,
             213: WordMorphism: 1->12, 2->2, 3->312,
             231: WordMorphism: 1->132, 2->2, 3->32,
             312: WordMorphism: 1->13, 2->213, 3->3,
             321: WordMorphism: 1->123, 2->23, 3->3}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {1: WordMorphism({1: [1,2,3], 2: [2], 3: [3]}),
                2: WordMorphism({1: [1], 2: [2,3,1], 3: [3]}),
                3: WordMorphism({1: [1], 2: [2], 3: [3,1,2]}),
                213: WordMorphism({1: [1, 2], 2: [2], 3: [3, 1, 2]}),
                123: WordMorphism({1: [1], 2: [2, 1], 3: [3, 2, 1]}),
                312: WordMorphism({1: [1, 3], 2: [2, 1, 3], 3: [3]}),
                132: WordMorphism({1: [1], 2: [2, 3, 1], 3: [3, 1]}),
                321: WordMorphism({1: [1, 2, 3], 2: [2, 3], 3: [3]}),
                231: WordMorphism({1: [1, 3, 2], 2: [2], 3: [3, 2]})}

cdef class ArnouxRauzy(MCFAlgorithm):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_pyx import ArnouxRauzy
        sage: algo = ArnouxRauzy()
        sage: TestSuite(algo).run()
        sage: algo._test_dual_substitution_definition()
        sage: algo._test_coherence()
        sage: algo._test_definition()
    """
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import ArnouxRauzy
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: ArnouxRauzy()(P)
            Traceback (most recent call last):
            ...
            ValueError: Arnoux is not defined on 
            ((0.3, 0.6, 0.8), (0.2, 0.3, 0.3))

        :: 

            sage: P = PairPoint(3, [.3,.2,.8], [.2,.3,.3])
            sage: ArnouxRauzy()(P)
            ((0.3, 0.2, 0.30000000000000004), (0.5, 0.6, 0.3))
        """
        if P.x[0] + P.x[1] < P.x[2]:
            P.x[2] -= P.x[0] + P.x[1]
            P.a[1] += P.a[2]
            P.a[0] += P.a[2]
            return 3
        elif P.x[0] + P.x[2] < P.x[1]:
            P.x[1] -= P.x[0] + P.x[2]
            P.a[2] += P.a[1]
            P.a[0] += P.a[1]
            return 2
        elif P.x[1] + P.x[2] < P.x[0]:
            P.x[0] -= P.x[1] + P.x[2]
            P.a[1] += P.a[0]
            P.a[2] += P.a[0]
            return 1
        else:
            raise ValueError('Arnoux is not defined on {}'.format(P))

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import ArnouxRauzy
            sage: ArnouxRauzy().substitutions()
            {1: WordMorphism: 1->1, 2->21, 3->31,
             2: WordMorphism: 1->12, 2->2, 3->32,
             3: WordMorphism: 1->13, 2->23, 3->3}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {1:  WordMorphism({1: [1], 2: [2, 1], 3: [3, 1]}),
                2:  WordMorphism({1: [1, 2], 2: [2], 3: [3, 2]}),
                3:  WordMorphism({1: [1, 3], 2: [2, 3], 3: [3]})}

    def dual_substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import ArnouxRauzy
            sage: ArnouxRauzy().dual_substitutions()
            {1: WordMorphism: 1->123, 2->2, 3->3,
             2: WordMorphism: 1->1, 2->231, 3->3,
             3: WordMorphism: 1->1, 2->2, 3->312}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {1: WordMorphism({1: [1,2,3], 2: [2], 3: [3]}),
                2: WordMorphism({1: [1], 2: [2,3,1], 3: [3]}),
                3: WordMorphism({1: [1], 2: [2], 3: [3,1,2]})}
cdef class Poincare(MCFAlgorithm):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_pyx import Poincare
        sage: algo = Poincare()
        sage: TestSuite(algo).run()
        sage: algo._test_dual_substitution_definition()
        sage: algo._test_coherence()
        sage: algo._test_definition()
    """
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Poincare
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: Poincare()(P)
            ((0.3, 0.3, 0.20000000000000007), (0.8, 0.6, 0.3))
        """
        return P._Poincare()

    def name(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Poincare
            sage: Poincare().name()
            "Poincar\\'e"
        """
        return r"Poincar\'e"

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Poincare
            sage: Poincare().substitutions()
            {123: WordMorphism: 1->123, 2->23, 3->3,
             132: WordMorphism: 1->132, 2->2, 3->32,
             213: WordMorphism: 1->13, 2->213, 3->3,
             231: WordMorphism: 1->1, 2->231, 3->31,
             312: WordMorphism: 1->12, 2->2, 3->312,
             321: WordMorphism: 1->1, 2->21, 3->321}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {312: WordMorphism({1: [1, 2], 2: [2], 3: [3, 1, 2]}),
                321: WordMorphism({1: [1], 2: [2, 1], 3: [3, 2, 1]}),
                213: WordMorphism({1: [1, 3], 2: [2, 1, 3], 3: [3]}),
                231: WordMorphism({1: [1], 2: [2, 3, 1], 3: [3, 1]}),
                123: WordMorphism({1: [1, 2, 3], 2: [2, 3], 3: [3]}),
                132: WordMorphism({1: [1, 3, 2], 2: [2], 3: [3, 2]})}

    def dual_substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Poincare
            sage: Poincare().dual_substitutions()
            {123: WordMorphism: 1->1, 2->21, 3->321,
             132: WordMorphism: 1->1, 2->231, 3->31,
             213: WordMorphism: 1->12, 2->2, 3->312,
             231: WordMorphism: 1->132, 2->2, 3->32,
             312: WordMorphism: 1->13, 2->213, 3->3,
             321: WordMorphism: 1->123, 2->23, 3->3}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {213: WordMorphism({1: [1, 2], 2: [2], 3: [3, 1, 2]}),
                123: WordMorphism({1: [1], 2: [2, 1], 3: [3, 2, 1]}),
                312: WordMorphism({1: [1, 3], 2: [2, 1, 3], 3: [3]}),
                132: WordMorphism({1: [1], 2: [2, 3, 1], 3: [3, 1]}),
                321: WordMorphism({1: [1, 2, 3], 2: [2, 3], 3: [3]}),
                231: WordMorphism({1: [1, 3, 2], 2: [2], 3: [3, 2]})}

cdef class Selmer(MCFAlgorithm):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_pyx import Selmer
        sage: algo = Selmer()
        sage: TestSuite(algo).run()
        sage: algo._test_dual_substitution_definition()
        sage: algo._test_coherence()
        sage: algo._test_definition()
    """
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Selmer
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: Selmer()(P)
            ((0.3, 0.6, 0.5), (0.5, 0.3, 0.3))
        """
        if P.x[0] <= P.x[1] <= P.x[2]:
            P.x[2] -= P.x[0]
            P.a[0] += P.a[2]
            return 123
        elif P.x[0] <= P.x[2] <= P.x[1]:
            P.x[1] -= P.x[0]
            P.a[0] += P.a[1]
            return 132
        elif P.x[1] <= P.x[2] <= P.x[0]:
            P.x[0] -= P.x[1]
            P.a[1] += P.a[0]
            return 231
        elif P.x[1] <= P.x[0] <= P.x[2]:
            P.x[2] -= P.x[1]
            P.a[1] += P.a[2]
            return 213
        elif P.x[2] <= P.x[0] <= P.x[1]:
            P.x[1] -= P.x[2]
            P.a[2] += P.a[1]
            return 312
        elif P.x[2] <= P.x[1] <= P.x[0]:
            P.x[0] -= P.x[2]
            P.a[2] += P.a[0]
            return 321
        else:
            raise ValueError('limit case: reach set of measure zero: {}'.format(P))

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Selmer
            sage: Selmer().substitutions()
            {123: WordMorphism: 1->13, 2->2, 3->3,
             132: WordMorphism: 1->12, 2->2, 3->3,
             213: WordMorphism: 1->1, 2->23, 3->3,
             231: WordMorphism: 1->1, 2->21, 3->3,
             312: WordMorphism: 1->1, 2->2, 3->32,
             321: WordMorphism: 1->1, 2->2, 3->31}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {132: WordMorphism({1: [1, 2], 2: [2], 3: [3]}),
                231: WordMorphism({1: [1], 2: [2, 1], 3: [3]}),
                123: WordMorphism({1: [1, 3], 2: [2], 3: [3]}),
                321: WordMorphism({1: [1], 2: [2], 3: [3, 1]}),
                213: WordMorphism({1: [1], 2: [2, 3], 3: [3]}),
                312: WordMorphism({1: [1], 2: [2], 3: [3, 2]})}

    def dual_substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Selmer
            sage: Selmer().dual_substitutions()
            {123: WordMorphism: 1->1, 2->2, 3->31,
             132: WordMorphism: 1->1, 2->21, 3->3,
             213: WordMorphism: 1->1, 2->2, 3->32,
             231: WordMorphism: 1->12, 2->2, 3->3,
             312: WordMorphism: 1->1, 2->23, 3->3,
             321: WordMorphism: 1->13, 2->2, 3->3}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {231: WordMorphism({1: [1, 2], 2: [2], 3: [3]}),
                132: WordMorphism({1: [1], 2: [2, 1], 3: [3]}),
                321: WordMorphism({1: [1, 3], 2: [2], 3: [3]}),
                123: WordMorphism({1: [1], 2: [2], 3: [3, 1]}),
                312: WordMorphism({1: [1], 2: [2, 3], 3: [3]}),
                213: WordMorphism({1: [1], 2: [2], 3: [3, 2]})}

cdef class FullySubtractive(MCFAlgorithm):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_pyx import FullySubtractive
        sage: algo = FullySubtractive()
        sage: TestSuite(algo).run()
        sage: algo._test_dual_substitution_definition()
        sage: algo._test_coherence()
        sage: algo._test_definition()
    """
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import FullySubtractive
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: FullySubtractive()(P)
            ((0.3, 0.3, 0.5), (0.8, 0.3, 0.3))
        """
        if P.x[0] <= P.x[1] and P.x[0] <= P.x[2]:
            P.x[1] -= P.x[0]
            P.x[2] -= P.x[0]
            P.a[0] += P.a[1] + P.a[2]
            return 1
        elif P.x[1] <= P.x[0] and P.x[1] <= P.x[2]:
            P.x[0] -= P.x[1]
            P.x[2] -= P.x[1]
            P.a[1] += P.a[0] + P.a[2]
            return 2
        elif P.x[2] <= P.x[0] and P.x[2] <= P.x[1]:
            P.x[0] -= P.x[2]
            P.x[1] -= P.x[2]
            P.a[2] += P.a[0] + P.a[1]
            return 3
        else:
            raise ValueError('limit case: reach set of measure zero: {}'.format(P))

    def name(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import FullySubtractive
            sage: FullySubtractive().name()
            'Fully Subtractive'
        """
        return r"Fully Subtractive"

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import FullySubtractive
            sage: FullySubtractive().substitutions()
            {1: WordMorphism: 1->123, 2->2, 3->3,
             2: WordMorphism: 1->1, 2->231, 3->3,
             3: WordMorphism: 1->1, 2->2, 3->312}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {1: WordMorphism({1: [1,2,3], 2: [2], 3: [3]}),
                2: WordMorphism({1: [1], 2: [2,3,1], 3: [3]}),
                3: WordMorphism({1: [1], 2: [2], 3: [3,1,2]})}

    def dual_substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import FullySubtractive
            sage: FullySubtractive().dual_substitutions()
            {1: WordMorphism: 1->1, 2->21, 3->31,
             2: WordMorphism: 1->12, 2->2, 3->32,
             3: WordMorphism: 1->13, 2->23, 3->3}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {1:  WordMorphism({1: [1], 2: [2, 1], 3: [3, 1]}),
                2:  WordMorphism({1: [1, 2], 2: [2], 3: [3, 2]}),
                3:  WordMorphism({1: [1, 3], 2: [2, 3], 3: [3]})}

cdef class Cassaigne(MCFAlgorithm):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_pyx import Cassaigne
        sage: algo = Cassaigne()
        sage: TestSuite(algo).run()
        sage: algo._test_dual_substitution_definition()
        sage: algo._test_coherence()
        sage: algo._test_definition()
    """
    cdef int call(self, PairPoint P) except *:
        r"""
        This algorithm was provided by Julien Cassaigne during a meeting of
        the ANR DynA3S on October 12, 2015 held in Paris. It is inspired
        from a method to generate words of complexity 2n+1 on a three
        letter alphabet of arbitrary letter frequencies.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Cassaigne
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: Cassaigne()(P)
            ((0.6, 0.3, 0.5), (0.3, 0.5, 0.3))
        """
        cdef double tmp
        if P.x[0] >= P.x[2] :
            P.x[0] -= P.x[2]
            tmp = P.x[1]
            P.x[1] = P.x[2]
            P.x[2] = tmp
            tmp = P.a[1]
            P.a[1] = P.a[0] + P.a[2]
            P.a[2] = tmp
            return 1
        elif P.x[0] < P.x[2] :
            P.x[2] -= P.x[0]
            tmp = P.x[1]
            P.x[1] = P.x[0]
            P.x[0] = tmp
            tmp = P.a[1]
            P.a[1] = P.a[0] + P.a[2]
            P.a[0] = tmp
            return 2
        else:
            raise ValueError('limit case: reach set of measure zero: {}'.format(P))

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Cassaigne
            sage: Cassaigne().substitutions()
            {1: WordMorphism: 1->1, 2->13, 3->2, 
             2: WordMorphism: 1->2, 2->13, 3->3}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {1:  WordMorphism({1: [1], 2: [1, 3], 3: [2]}),
                2:  WordMorphism({1: [2], 2: [1, 3], 3: [3]})}

    def dual_substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Cassaigne
            sage: Cassaigne().dual_substitutions()
            {1: WordMorphism: 1->12, 2->3, 3->2, 
             2: WordMorphism: 1->2, 2->1, 3->23}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {1:  WordMorphism({1: [1,2], 2: [3], 3: [2]}),
                2:  WordMorphism({1: [2], 2: [1], 3: [2,3]})}

cdef class Sorted_Brun(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_Brun
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: Sorted_Brun()(P)
            ((0.20000000000000007, 0.3, 0.6), (0.3, 0.2, 0.6))

        ::

            sage: P = PairPoint(3, [.3,.45,.8], [.2,.3,.3])
            sage: D = {'x':.3,'y':.45,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: Sorted_Brun()(P)
            ((0.3, 0.35000000000000003, 0.45), (0.2, 0.3, 0.6))

         ::

            sage: P = PairPoint(3, [.3,.3,.8], [.2,.3,.3])
            sage: Sorted_Brun()(P)
            ((0.3, 0.3, 0.5), (0.2, 0.6, 0.3))

        """
        P.x[2] -= P.x[1]
        P.a[1] += P.a[2]
        P.sort()
        return 100

cdef class Sorted_BrunMulti(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_BrunMulti
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.3,.8], [.2,.3,.3])
            sage: Sorted_BrunMulti()(P)
            ((0.20000000000000007, 0.3, 0.3), (0.3, 0.2, 0.8999999999999999))
        """
        cdef int m = <int>(P.x[2] / P.x[1])
        P.x[2] -= m*P.x[1]
        P.a[1] += m*P.a[2]
        P.sort()
        return 100

cdef class Sorted_Selmer(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_Selmer
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: Sorted_Selmer()(P)
            ((0.3, 0.5, 0.6), (0.5, 0.3, 0.3))
        """
        P.x[2] -= P.x[0]
        P.a[0] += P.a[2]
        P.sort()
        return 100

cdef class Sorted_FullySubtractive(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_FullySubtractive
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.6,.8], [.2,.3,.3])
            sage: Sorted_FullySubtractive()(P)
            ((0.3, 0.3, 0.5), (0.8, 0.3, 0.3))
        """
        # Apply the algo
        P.x[1] -= P.x[0]
        P.x[2] -= P.x[0]
        P.a[0] += P.a[1] + P.a[2]
        P.sort()
        return 100

cdef class Sorted_ArnouxRauzy(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_ArnouxRauzy
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.4,.8], [.2,.3,.4])
            sage: Sorted_ArnouxRauzy()(P)  # tol
            ((0.1, 0.3, 0.4), (0.4, 0.6, 0.7))

        ::

            sage: P = PairPoint(3, [.3,.7,.8], [.2,.3,.4])
            sage: Sorted_ArnouxRauzy()(P)
            ((-0.19999999999999996, 0.3, 0.7), (0.4, 0.6000000000000001, 0.7))

        """
        #Arnoux-Rauzy
        P.x[2] -= P.x[0] + P.x[1]
        P.a[0] += P.a[2]
        P.a[1] += P.a[2]
        P.sort()
        return 100

cdef class Sorted_ArnouxRauzyMulti(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

        """
        #Arnoux-Rauzy Multi
        cdef int m
        m = <int>(P.x[2] / (P.x[0] + P.x[1]))
        P.x[2] -= m * (P.x[0] + P.x[1])
        P.a[1] += m * P.a[2];
        P.a[0] += m * P.a[2];
        P.sort()
        return 100

cdef class Sorted_ARP(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_ARP
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.4,.8], [.2,.3,.4])
            sage: Sorted_ARP()(P)        # tol
            ((0.1, 0.3, 0.4), (0.4, 0.6, 0.7))

        ::

            sage: P = PairPoint(3, [.3,.7,.8], [.2,.3,.4])
            sage: Sorted_ARP()(P)       # tol
            ((0.1, 0.3, 0.4), (0.4, 0.9, 0.7))
        """
        # Apply the algo
        if P.x[2] > P.x[0] + P.x[1]:
            return P._Sorted_ArnouxRauzy()
        else:
            return P._Sorted_Poincare()

cdef class Sorted_ARPMulti(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_ARPMulti
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.5,.8], [.2,.3,.3])
            sage: Sorted_ARPMulti()(P)
            ((0.2, 0.3, 0.30000000000000004), (0.6, 0.8, 0.3))
        """
        # Apply the algo
        if P.x[2] > P.x[0] + P.x[1]:
            return P._Sorted_ArnouxRauzyMulti()
        else:
            return P._Sorted_Poincare()

cdef class Sorted_Poincare(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_Poincare
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.7,.8], [.2,.3,.4])
            sage: Sorted_Poincare()(P)       # tol
            ((0.1, 0.3, 0.4), (0.4, 0.9, 0.7))

        """
        # Apply the algo
        P.x[2] -= P.x[1]
        P.x[1] -= P.x[0]
        P.a[1] += P.a[2]
        P.a[0] += P.a[1]
        P.sort()
        return 200

cdef class Sorted_ARrevert(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_ARrevert
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.3,.8], [.2,.3,.3])
            sage: Sorted_ARrevert()(P)
            ((0.2, 0.3, 0.3), (0.3, 0.5, 0.6))
        """
        cdef PairPoint R = PairPoint(self.dim)
        cdef double z_x_y = P.x[2] - P.x[0] - P.x[1]
        if z_x_y > P.x[1]:
            R.x[0] = P.x[0]
            R.x[1] = P.x[1]
            R.x[2] = z_x_y
            R.a[0] = P.a[0] + P.a[2]
            R.a[1] = P.a[1] + P.a[2]
            R.a[2] = P.a[2]
            P.copy_inplace(R)
            return 100
        elif z_x_y > P.x[0]:
            R.x[0] = P.x[0]
            R.x[1] = z_x_y
            R.x[2] = P.x[1]
            R.a[0] = P.a[0] + P.a[2]
            R.a[1] = P.a[2]
            R.a[2] = P.a[1] + P.a[2]
            P.copy_inplace(R)
            return 200
        elif z_x_y > 0:
            R.x[0] = z_x_y
            R.x[1] = P.x[0]
            R.x[2] = P.x[1]
            R.a[0] = P.a[2]
            R.a[1] = P.a[0] + P.a[2]
            R.a[2] = P.a[1] + P.a[2]
            P.copy_inplace(R)
            return 300
        else:
            # Revert
            R.x[0] = (P.x[0] + P.x[1] - P.x[2])/2
            R.x[1] = (P.x[0] - P.x[1] + P.x[2])/2
            R.x[2] = (-P.x[0] + P.x[1] + P.x[2])/2
            R.a[0] = P.a[0] + P.a[1]
            R.a[1] = P.a[0] + P.a[2]
            R.a[2] = P.a[1] + P.a[2]
            P.copy_inplace(R)
            return 400

cdef class Sorted_ARrevertMulti(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_ARrevertMulti
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.3,.8], [.2,.3,.3])
            sage: Sorted_ARrevertMulti()(P)
            ((0.20000000000000007, 0.3, 0.3), (0.3, 0.5, 0.6))
        """
        cdef PairPoint R = PairPoint(self.dim)
        cdef int m = <int>(P.x[2] / (P.x[0] + P.x[1]))
        cdef double z_mxy = P.x[2] - m * (P.x[0] + P.x[1])
        if m == 0:
            # Revert
            R.x[0] = (P.x[0] + P.x[1] - P.x[2])/2
            R.x[1] = (P.x[0] - P.x[1] + P.x[2])/2
            R.x[2] = (-P.x[0] + P.x[1] + P.x[2])/2
            R.a[0] = P.a[0] + P.a[1]
            R.a[1] = P.a[0] + P.a[2]
            R.a[2] = P.a[1] + P.a[2]
            P.copy_inplace(R)
            return 100
        elif z_mxy > P.x[1]:
            R.x[0] = P.x[0]
            R.x[1] = P.x[1]
            R.x[2] = z_mxy
            R.a[0] = P.a[0] + m*P.a[2]
            R.a[1] = P.a[1] + m*P.a[2]
            R.a[2] = P.a[2]
            P.copy_inplace(R)
            return 200
        elif z_mxy > P.x[0]:
            R.x[0] = P.x[0]
            R.x[1] = z_mxy
            R.x[2] = P.x[1]
            R.a[0] = P.a[0] + m*P.a[2]
            R.a[1] = P.a[2]
            R.a[2] = P.a[1] + m*P.a[2]
            P.copy_inplace(R)
            return 300
        else:
            R.x[0] = z_mxy
            R.x[1] = P.x[0]
            R.x[2] = P.x[1]
            R.a[0] = P.a[2]
            R.a[1] = P.a[0] + m*P.a[2]
            R.a[2] = P.a[1] + m*P.a[2]
            P.copy_inplace(R)
            return 400

cdef class Sorted_ARMonteil(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_ARMonteil
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.3,.8], [.2,.3,.3])
            sage: Sorted_ARMonteil()(P)         # tol
            ((0.2, 0.3, 0.3), (0.3, 0.5, 0.6))
        """
        cdef PairPoint R = PairPoint(self.dim)

        # Apply the algo
        if P.x[2] > P.x[0] + P.x[1]:
            return P._Sorted_ArnouxRauzy()
        else:
            # Monteil
            R.x[0] = P.x[0] + P.x[1] - P.x[2]
            R.x[1] = -P.x[0] + P.x[2]
            R.x[2] = -P.x[1] + P.x[2]
            R.a[0] = P.a[0] + P.a[1] + P.a[2]
            R.a[1] = P.a[1] + P.a[2]
            R.a[2] = P.a[0] + P.a[2]
            P.copy_inplace(R)
            P.sort()
            return 200

cdef class Sorted_Delaunay(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        Donné par Xavier Provençal (inspiré de en fait) le 3 février 2014.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import Sorted_Delaunay
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.3,.8], [.2,.3,.3])
            sage: Sorted_Delaunay()(P)
            ((0.0, 0.3, 0.8), (0.6, 0.5, 0.3))
        """
        cdef PairPoint R = PairPoint(self.dim)
        # Apply the algo
        if P.x[2] > P.x[0] + P.x[1]:
            # Genre de semi revert
            R.x[0] = P.x[0]
            R.x[1] = P.x[1] - P.x[0]
            R.x[2] = P.x[0] - P.x[1] + P.x[2]
            R.a[0] = P.a[0] + P.a[1]
            R.a[1] = P.a[1] + P.a[2]
            R.a[2] = P.a[2]
            P.copy_inplace(R)
            P.sort()
            return 200
        else:
            return P._Sorted_FullySubtractive()
cdef class JacobiPerron(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import JacobiPerron
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.3,.8], [.2,.3,.3])
            sage: JacobiPerron()(P)
            ((0.0, 0.20000000000000007, 0.3), (0.3, 0.3, 1.0999999999999999))
        """
        cdef int m,n                # temporary integer variables
        cdef double r,s,t           # temporary variables

        # Apply the algo
        m = int(P.x[2] / P.x[0])
        n = int(P.x[1] / P.x[0])
        t = P.x[2] - m*P.x[0]
        s = P.x[1] - n*P.x[0]
        r = P.x[0]

        P.x[2] = r
        P.x[1] = t
        P.x[0] = s

        t = P.a[2]
        s = P.a[1]
        r = m*P.a[2] + n*P.a[1] + P.a[0]

        P.a[2] = r
        P.a[1] = t
        P.a[0] = s

        return 100
cdef class JacobiPerronAdditif(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import JacobiPerronAdditif
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.3,.8], [.2,.3,.3])
            sage: JacobiPerronAdditif()(P)
            ((0.3, 0.3, 0.5), (0.5, 0.3, 0.3))
        """
        cdef PairPoint R = PairPoint(self.dim)
        # Apply the algo
        if P.x[0] < P.x[1]:
            R.x[0] = P.x[0]
            R.x[1] = P.x[1] - P.x[0]
            R.x[2] = P.x[2]
            R.a[0] = P.a[0] + P.a[1]
            R.a[1] = P.a[1]
            R.a[2] = P.a[2]
            P.copy_inplace(R)
            return 100
        elif P.x[0] < P.x[2]:
            R.x[0] = P.x[0]
            R.x[1] = P.x[1]
            R.x[2] = P.x[2] - P.x[0]
            R.a[0] = P.a[0] + P.a[2]
            R.a[1] = P.a[1]
            R.a[2] = P.a[2]
            P.copy_inplace(R)
            return 200
        elif P.x[0] > P.x[1] and P.x[0] > P.x[2]:
            R.x[0] = P.x[1]
            R.x[1] = P.x[2]
            R.x[2] = P.x[0]
            R.a[0] = P.a[1]
            R.a[1] = P.a[2]
            R.a[2] = P.a[0]
            P.copy_inplace(R)
            return 300
        else:
            raise ValueError("jacobi not defined for (x,y,z)=(%s,%s,%s)"%(P.x[0],P.x[1],P.x[2]))
cdef class JacobiPerronAdditifv2(MCFAlgorithm):
    cdef int call(self, PairPoint P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_pyx import JacobiPerronAdditifv2
            sage: from slabbe.mult_cont_frac_pyx import PairPoint
            sage: P = PairPoint(3, [.3,.3,.8], [.2,.3,.3])
            sage: JacobiPerronAdditifv2()(P)
            ((0.3, 0.3, 0.5), (0.5, 0.3, 0.3))
        """
        cdef PairPoint R = PairPoint(self.dim)
        # Apply the algo
        if P.x[0] < P.x[1]:
            R.x[0] = P.x[0]
            R.x[1] = P.x[1] - P.x[0]
            R.x[2] = P.x[2] - P.x[0]
            R.a[0] = P.a[0] + P.a[1] + P.a[2]
            R.a[1] = P.a[1]
            R.a[2] = P.a[2]
            P.copy_inplace(R)
            return 100
        elif P.x[0] < P.x[2]:
            R.x[0] = P.x[0]
            R.x[1] = P.x[1]
            R.x[2] = P.x[2] - P.x[0]
            R.a[0] = P.a[0] + P.a[2]
            R.a[1] = P.a[1]
            R.a[2] = P.a[2]
            P.copy_inplace(R)
            return 200
        elif P.x[0] > P.x[1] and P.x[0] > P.x[2]:
            R.x[0] = P.x[1]
            R.x[1] = P.x[2]
            R.x[2] = P.x[0]
            R.a[0] = P.a[1]
            R.a[1] = P.a[2]
            R.a[2] = P.a[0]
            P.copy_inplace(R)
            return 300
        else:
            raise ValueError("jacobi not defined for (x,y,z)=(%s,%s,%s)"%(P.x[0],P.x[1],P.x[2]))

cdef inline (double, double) projection3to2(double x, double y, double z):
    cdef double s = -SQRT3SUR2 * x + SQRT3SUR2 * y
    cdef double t = -.5 * x -.5 * y + z
    return s,t

cdef inline double max_array(double[:] tab, int dim):
    cdef double s = 0
    cdef int i
    for i in range(dim):
        s = max(s, abs(tab[i]))
    return s

cdef inline double sum_array(double[:] tab, int dim):
    cdef double s = 0
    cdef int i
    for i in range(dim):
        s += abs(tab[i])
    return s
