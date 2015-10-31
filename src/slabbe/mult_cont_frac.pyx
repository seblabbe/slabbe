# -*- coding: utf-8 -*-
r"""
Multidimensional Continued Fraction Algorithms

EXAMPLES::

    sage: from slabbe.mult_cont_frac import Brun
    sage: Brun().natural_extension_plot(3000, norm_left='1', axis_off=True, savefig=False)
    <matplotlib.figure.Figure object at ...>

Construction of an s-adic word::

    sage: from slabbe.mult_cont_frac import ARP
    sage: from itertools import repeat
    sage: D = ARP().substitutions()
    sage: it = ARP().coding_iterator((1,e,pi))
    sage: words.s_adic(it, repeat(1), D)
    word: 1232323123233231232332312323123232312323...

It can be done with one line::

    sage: ARP().s_adic_word((1,e,pi))
    word: 1232323123233231232332312323123232312323...

TODO:

    - Ajout les algo de reuteneaour, nogueira, autres?
    - make a other class for 2d and and 1d methods
    - faire une fonction max
    - utilise les vecteurs pour les plus grandes dimensions?
    - or avoid creating a new pairpoint R, the copy is already done by
      default

    - Read [1] and change cpdef int C[NDIVS][NDIVS][NDIVS]
    - Replace method ``natural_extention_dict`` by ``orbit_filtered_list``
    - Use ``orbit_filtered_list`` for ``invariant_measure`` ?

    - Essayer d'utiliser @staticmethod pour call pour que ARP puisse
      apeler Poincare. Without the cython Error: Cannot call a static
      method on an instance variable.
      https://groups.google.com/d/topic/sage-support/DRI_s31D8ks/discussion

    - invariant_measure_dict ne fonctionne pas bien avec norm 1
    - rapatrier w0_distance_tijdeman
    - rapatrier plot2d_distance
    - rapatrier save_latex_and_dat_tables

Question:

    - Comment factoriser le code sans utiliser les yield?
    - Comment faire un appel de fonction rapide (pour factoriser le code)

[1] https://groups.google.com/forum/?fromgroups=#!topic/sage-devel/NCBmj2KjwEM

Gain d'optimisation in Cython:

    - Declare x,y,z as double instead of list : factor of 4
    - Faire les calculs dans une boucle : factor of 2
    - Do not use yield : factor of 10

AUTHORS:

 - Vincent Delecroix, C code, Computation of Lyapounov exponents for Brun
   algorithm, June 2013.
 - Sebastien Labbe, Invariant measures, Lyapounov exponents and natural
   extensions for a dozen of algorithms, October 2013.

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
from sage.misc.prandom import random

#include "sage/ext/cdefs.pxi"
#include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"  # ctrl-c interrupt block support

cdef struct PairPoint3d:
    double x
    double y
    double z
    double u
    double v
    double w
    int branch

cdef double SQRT3SUR2 = 0.866025403784439

PGF_COLORS = ["red", "green", "blue", "cyan", "brown", "gray", "orange", "pink",
"yellow", "black", "white", "darkgray", "lightgray",
"lime", "olive", "magenta", "purple", "teal", "violet"]

cdef class MCFAlgorithm(object):
    ########################################
    # METHODS IMPLEMENTED IN HERITED CLASSES
    ########################################
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        This method must be implemented in the inherited classes.

        EXAMPLES::

        """
        raise NotImplementedError
    def substitutions(self):
        r"""
        This method must be implemented in the inherited classes.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
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

            sage: from slabbe.mult_cont_frac import Brun
            sage: Brun().dual_substitutions()
            {123: WordMorphism: 1->1, 2->23, 3->3,
             132: WordMorphism: 1->1, 2->2, 3->32,
             213: WordMorphism: 1->13, 2->2, 3->3,
             231: WordMorphism: 1->1, 2->2, 3->31,
             312: WordMorphism: 1->12, 2->2, 3->3,
             321: WordMorphism: 1->1, 2->21, 3->3}

        """
        raise NotImplementedError

    def branches(self, int n_iterations=1000):
        r"""
        Returns the branches labels of the algorithm.

        This method is an heuristic and should be implemented in the
        inherited classes.

        EXAMPLES::

            sage: import slabbe.mult_cont_frac as mcf
            sage: mcf.Brun().branches()
            {123, 132, 213, 231, 312, 321}
            sage: mcf.ARP().branches()
            {1, 2, 3, 123, 132, 213, 231, 312, 321}
        """
        cdef unsigned int i         # loop counter
        cdef PairPoint3d P
        S = set()
        # Loop
        for i from 0 <= i < n_iterations:

            # Check for Keyboard interupt
            sig_check()

            # random initial values
            P.x = random(); P.y = random(); P.z = random();
            P.u = random(); P.v = random(); P.w = random();

            # Apply Algo
            P = self.call(P)

            S.add(P.branch)
        return S
    ######################
    # TEST METHODS 
    ######################
    def _test_definition(self, int n_iterations=10000):
        r"""
        INPUT:

        - ``n_iterations`` -- integer

        OUTPUT:

        bool

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ARP, Brun
            sage: t = ARP()._test_definition(10000)
            sage: t = Brun()._test_definition(10000)

        """
        cdef double s,t             # temporary variables
        cdef unsigned int i         # loop counter

        cdef PairPoint3d P, R

        # Loop
        for i from 0 <= i < n_iterations:

            # Check for Keyboard interupt
            sig_check()

            # random initial values
            P.x = random(); P.y = random(); P.z = random();
            P.u = random(); P.v = random(); P.w = random();
            s = P.x*P.u + P.y*P.v + P.z*P.w

            # Apply Algo
            R = self.call(P)
            t = R.x*R.u + R.y*R.v + R.z*R.w

            if not abs(s - t) < 0.0000001:
                m = 'This algo does not preserve the scalar product\n'
                m += '{} != {}\n'.format(s,t)
                m += 'The problem is on branch {}\n'.format(R.branch)
                m += 'INPUT: ({}, {}, {}, {}, {}, {})\n'.format(P.x,P.y,P.z,P.u,P.v,P.w)
                m += 'OUTPUT: ({}, {}, {}, {}, {}, {})\n'.format(R.x,R.y,R.z,R.u,R.v,R.w)
                raise Exception(m)

        return True

    def _test_substitution_definition(self):
        r"""
        OUTPUT:

        bool

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ARP, Brun
            sage: t = ARP()._test_substitution_definition()
            sage: t = Brun()._test_substitution_definition()

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
        return True

    ######################
    # METHODS FOR THE USER:
    ######################
    def name(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Reverse, Brun
            sage: Reverse().name()
            'Reverse'
            sage: Brun().name()
            'Brun'
        """
        return self.__class__.__name__
    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Reverse, Brun
            sage: Reverse()
            Reverse 3-dimensional continued fraction algorithm
            sage: Brun()
            Brun 3-dimensional continued fraction algorithm
        """
        return "{} 3-dimensional continued fraction algorithm".format(self.name())

    def __call__(self, PairPoint3d P):
        r"""
        Wrapper for the cdef call method.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
            sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Brun()(D)
            sage: sorted(E.iteritems())
            [('branch', 123),
             ('u', 0.2),
             ('v', 0.6),
             ('w', 0.3),
             ('x', 0.3),
             ('y', 0.6),
             ('z', 0.20000000000000007)]

        """
        return self.call(P)
    def matrix_cocycle(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ARP
            sage: ARP().matrix_cocycle()
            Cocycle with 9 gens over Regular language over ['A1', 'A2', 'A3',
            'P12', 'P13', 'P21', 'P23', 'P31', 'P32']
            defined by: Automaton with 7 states

        ::

            sage: from slabbe.mult_cont_frac import Sorted_Brun
            sage: Sorted_Brun().matrix_cocycle()
            Cocycle with 3 gens over Language of finite words over alphabet
            ['B1', 'B2', 'B3']

        ::

            sage: from slabbe.mult_cont_frac import Brun
            sage: Brun().matrix_cocycle()
            Cocycle with 6 gens over Regular language over 
            ['B12', 'B13', 'B21', 'B23', 'B31', 'B32']
            defined by: Automaton with 6 states
        """
        from matrix_cocycle import cocycles
        try:
            f = getattr(cocycles, self.name())
        except AttributeError:
            msg = "Matrix cocyle not implemented for {}"
            msg = msg.format(self.name())
            raise NotImplementedError(msg)
        return f()

    def coding_iterator(self, start):
        r"""
        INPUT:

        - ``start`` -- iterable of three real numbers

        OUTPUT:

            iterator

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ARP
            sage: it = ARP().coding_iterator((1,e,pi))
            sage: [next(it) for _ in range(20)]
            [123, 2, 1, 123, 1, 231, 3, 3, 3, 3, 123, 1, 1, 1, 231, 2, 321, 2, 3, 312]

        """
        cdef double s             # temporary variables
        cdef PairPoint3d P

        # initial values
        P.x, P.y, P.z = start
        P.u = random(); P.v = random(); P.w = random();

        # Normalize (x,y,z)
        s = P.x + P.y + P.z
        P.x /= s; P.y /= s; P.z /= s

        # Loop
        while True:
            P = self.call(P)
            yield P.branch

            # Normalize (x,y,z)
            s = P.x + P.y + P.z
            P.x /= s; P.y /= s; P.z /= s

    def orbit_iterator(self, start=None, norm_left='sup', norm_right='1'):
        r"""
        INPUT:

        - ``start`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``norm_left`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit of the algo
        - ``norm_right`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'`` or ``'hypersurfac'``, the norm used for the orbit in the natural extension

        NOTE:

            This iterator is 10x slower because of the yield statement. So
            avoid using this when writing fast code. Just copy paste the
            loop or use orbit_list or orbit_filtered_list method.

        OUTPUT:

            iterator

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
            sage: it = Brun().orbit_iterator((.414578,.571324,.65513))
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

            sage: from slabbe.mult_cont_frac import Brun
            sage: it = Brun().orbit_iterator((.414578,.571324,.65513), norm_left='1')
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
        cdef PairPoint3d P
        if start is None:
            P.x = random(); P.y = random(); P.z = random()
        else:
            P.x = start[0]; P.y = start[1]; P.z = start[2]
        P.u = 1./3
        P.v = 1./3
        P.w = 1./3

        # Normalize (x,y,z)
        s = P.x + P.y + P.z
        P.x /= s; P.y /= s; P.z /= s

        # Loop
        while True:

            # Apply Algo
            P = self.call(P)

            # Normalize (xnew,ynew,znew)
            if norm_left == '1':
                s = P.x + P.y + P.z # norm 1
            elif norm_left == 'sup':
                s = max(P.x, P.y, P.z) # norm sup
            else:
                raise ValueError("Unknown value for norm_left(=%s)" %norm_left)
            P.x /= s; P.y /= s; P.z /= s

            # Normalize (unew,vnew,wnew)
            if norm_right == '1':
                s = P.u + P.v + P.w    # norm 1
            elif norm_right == 'sup':
                s = max(P.u, P.v, P.w) # norm sup
            elif norm_right == 'hypersurface':
                s = P.x*P.u + P.y*P.v + P.z*P.w # hypersurface
            else:
                raise ValueError("Unknown value for norm_right(=%s)" %norm_right)
            P.u /= s; P.v /= s; P.w /= s

            yield (P.x, P.y, P.z), (P.u, P.v, P.w), P.branch

    def s_adic_word(self, v=None, first_letter=1):
        r"""
        INPUT:

        - ``v`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``first_letter`` - letter (default: ``1``)

        OUTPUT:

            word

        .. NOTE::

            The code is currently based on ``coding_iterator`` method. It
            could be made faster if based on ``orbit_list`` instead which
            is not using ``yield`` statements.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun, ARP
            sage: Brun().s_adic_word((.414578,.571324,.65513))
            word: 1232312312323123123312323123123312323123...
            sage: Brun().s_adic_word((1,e,pi))
            word: 1232323123233231232332312323123232312323...
            sage: ARP().s_adic_word((1,e,pi))
            word: 1232323123233231232332312323123232312323...

        Cassaigne substitutions do not allow an s-adic construction::

            sage: Cassaigne().s_adic_word((1,e,pi))
            Traceback (most recent call last):
            ...
            ValueError: impossible choice of letter for 6-th iteration in
            the s-adic construction

        TESTS::

            sage: v = ARP().s_adic_word((1,e,pi))
            sage: w = Brun().s_adic_word((1,e,pi))
            sage: v.longest_common_prefix(w, 'finite').length()
            212
        """
        from sage.combinat.words.word_generators import words
        if v is None:
            v = (random(), random(), random())
        it = self.coding_iterator(v)
        letters_it = self._s_adic_word_letter_iterator(v, first_letter)
        D = self.substitutions()
        return words.s_adic(it, letters_it, D)

    def _s_adic_word_letter_iterator(self, v, first_letter=1):
        r"""
        INPUT:

        - ``v`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``first_letter`` - letter (default: ``1``)

        OUTPUT:

            iterator of letters

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun, ARP
            sage: it = Brun()._s_adic_word_letter_iterator((.414578,.571324,.65513), 1)
            sage: [next(it) for _ in range(20)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        ::

            sage: it = ARP()._s_adic_word_letter_iterator((.414578,.571324,.65513), 2)
            sage: [next(it) for _ in range(20)]
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

        ::

            sage: it = Reverse()._s_adic_word_letter_iterator((.414578,.571324,.65513), 2)
            sage: [next(it) for _ in range(20)]
            [1, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

        ::

            sage: it = Cassaigne()._s_adic_word_letter_iterator((.414578,.571324,.65513), 2)
            sage: [next(it) for _ in range(20)]
            Traceback (most recent call last):
            ...
            ValueError: impossible choice of letter for 6-th iteration in
            the s-adic construction
        """
        from collections import defaultdict
        from sage.misc.prandom import choice
        from should_be_in_sage import image_first_letter_to_letter
        D = self.substitutions()
        # Construct the letter to letter function
        F = {key:image_first_letter_to_letter(s) for key,s in D.iteritems()}
        # Build the letters iterator
        cur_letter = first_letter
        for i,key in enumerate(self.coding_iterator(v)):
            possible = F[key][cur_letter]
            if not possible:
                raise ValueError("impossible choice of letter for {}-th "
                        "iteration in the s-adic construction".format(i))
            cur_letter = choice(possible)
            yield cur_letter

    def orbit_list(self, int n_iterations, start=None, norm_left='sup', norm_right='1'):
        r"""
        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``start`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``norm_left`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit of the algo
        - ``norm_right`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'`` or ``'hypersurfac'``, the norm used for the orbit in the natural extension

        OUTPUT:

            list

        BENCHMARK:

        It could be 10 times faster because 10^6 iterations can be done in
        about 60ms on this machine. But for drawing images, it does not
        matter to be 10 times slower::

            sage: %time L = Brun().orbit_list(10^6)   # not tested
            CPU times: user 376 ms, sys: 267 ms, total: 643 ms
            Wall time: 660 ms

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
            sage: L = Brun().orbit_list(10^5)
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
        cdef PairPoint3d P
        cdef int i
        if start is None:
            P.x = random(); P.y = random(); P.z = random()
        else:
            P.x = start[0]; P.y = start[1]; P.z = start[2]
        P.u = 1./3
        P.v = 1./3
        P.w = 1./3

        # Normalize (x,y,z)
        s = P.x + P.y + P.z
        P.x /= s; P.y /= s; P.z /= s

        L = []

        # Loop
        for i from 0 <= i < n_iterations:

            # Check for Keyboard interupt
            sig_check()

            # Apply Algo
            P = self.call(P)

            # Normalize (xnew,ynew,znew)
            if norm_left == '1':
                s = P.x + P.y + P.z # norm 1
            elif norm_left == 'sup':
                s = max(P.x, P.y, P.z) # norm sup
            else:
                raise ValueError("Unknown value for norm_left(=%s)" %norm_left)
            P.x /= s; P.y /= s; P.z /= s

            # Normalize (unew,vnew,wnew)
            if norm_right == '1':
                s = P.u + P.v + P.w    # norm 1
            elif norm_right == 'sup':
                s = max(P.u, P.v, P.w) # norm sup
            elif norm_right == 'hypersurface':
                s = P.x*P.u + P.y*P.v + P.z*P.w # hypersurface
            else:
                raise ValueError("Unknown value for norm_right(=%s)" %norm_right)
            P.u /= s; P.v /= s; P.w /= s

            L.append( (P.x, P.y, P.z, P.u, P.v, P.w, P.branch))

        return L

    def e_one_star_patch(self, v, n):
        r"""
        Return the n-th iterated patch of normal vector v.

        INPUT:

        - ``v`` -- vector, the normal vector
        - ``n`` -- integer

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Meester
            sage: Meester().e_one_star_patch((1,e,pi), 4)
            Patch of 21 faces
        """
        from sage.combinat.e_one_star import E1Star, Patch, Face
        from sage.misc.misc_c import prod
        if v is None:
            v = (random(), random(), random())
        it = self.coding_iterator(v)
        keys = [next(it) for _ in range(n)]
        D = self.dual_substitutions()
        L = prod(D[key] for key in reversed(keys))
        dual_sub = E1Star(L)
        cube = Patch([Face((1,0,0),1), Face((0,1,0),2), Face((0,0,1),3)])
        return dual_sub(cube)

    def orbit_filtered_list(self, int n_iterations, start=None,
            norm_left='1', norm_right='1',
            double xmin=-float('inf'), double xmax=float('inf'),
            double ymin=-float('inf'), double ymax=float('inf'),
            double umin=-float('inf'), double umax=float('inf'),
            double vmin=-float('inf'), double vmax=float('inf'),
            int ndivs=0):
        r"""
        Return a list of the orbit filtered to fit into a rectangle.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``start`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``norm_left`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit of the algo
        - ``norm_right`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'`` or ``'hypersurfac'``, the norm used for the orbit in the natural extension
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

        BENCHMARK:

            sage: from slabbe.mult_cont_frac import Brun
            sage: %time D = Brun().orbit_filtered_list(10^6) # not tested
            CPU times: user 366 ms, sys: 203 ms, total: 568 ms
            Wall time: 570 ms

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
            sage: D = Brun().orbit_filtered_list(3, start=(.414578,.571324,.65513))
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

            sage: Brun().orbit_filtered_list(3, norm_left='1',ndivs=1000)
            Traceback (most recent call last):
            ...
            ValueError: when ndivs is specified, you must provide a value
            for xmin, xmax, ymin, ymax, umin, umax, vmin and vmax

        ::

            sage: Brun().orbit_filtered_list(7, norm_left='1', ndivs=100, # random
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
        cdef PairPoint3d P
        cdef int previous_branch
        cdef int xa,ya,ua,va
        cdef int i
        cdef double xlen = xmax - xmin
        cdef double ylen = ymax - ymin
        cdef double ulen = umax - umin
        cdef double vlen = vmax - vmin
        if start is None:
            P.x = random(); P.y = random(); P.z = random()
        else:
            P.x = start[0]; P.y = start[1]; P.z = start[2]
        P.u = 1./3
        P.v = 1./3
        P.w = 1./3

        if ndivs and float('inf') in [-xmin, -ymin, -umin, -vmin, xmax, ymax, umax, vmax]:
            raise ValueError("when ndivs is specified, you must provide a"
                    " value for xmin, xmax, ymin, ymax, umin, umax, vmin"
                    " and vmax")

        # Normalize (x,y,z)
        s = P.x + P.y + P.z
        P.x /= s; P.y /= s; P.z /= s

        L = []

        # Apply Algo once
        P = self.call(P)

        # Loop
        for i from 0 <= i < n_iterations:

            # Check for Keyboard interupt
            sig_check()

            # Normalize (xnew,ynew,znew)
            if norm_left == '1':
                s = P.x + P.y + P.z # norm 1
            elif norm_left == 'sup':
                s = max(P.x, P.y, P.z) # norm sup
            else:
                raise ValueError("Unknown value for norm_left(=%s)" %norm_left)
            P.x /= s; P.y /= s; P.z /= s

            # Normalize (unew,vnew,wnew)
            if norm_right == '1':
                s = P.u + P.v + P.w    # norm 1
            elif norm_right == 'sup':
                s = max(P.u, P.v, P.w) # norm sup
            elif norm_right == 'hypersurface':
                s = P.x*P.u + P.y*P.v + P.z*P.w # hypersurface
            else:
                raise ValueError("Unknown value for norm_right(=%s)" %norm_right)
            P.u /= s; P.v /= s; P.w /= s

            # Projection
            if norm_left == '1':
                x = -SQRT3SUR2 * P.x + SQRT3SUR2 * P.y
                y = -.5 * P.x -.5 * P.y + P.z
            elif norm_left == 'sup':
                x = P.x
                y = P.y

            if norm_right == '1':
                u = -SQRT3SUR2 * P.u + SQRT3SUR2 * P.v
                v = -.5 * P.u -.5 * P.v + P.w
            elif norm_right == 'sup':
                u = P.u
                v = P.v

            # Apply Algo
            previous_branch = P.branch
            P = self.call(P)

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
                L.append( (xa,ya,ua,va, previous_branch, P.branch))
            else:
                L.append( (x,y,u,v, previous_branch, P.branch))

        return L

    def invariant_measure_dict(self, int n_iterations, int ndivs, v=None,
            str norm='1', verbose=False):
        r"""
        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``ndvis`` - integer (less than or equal to 256), number of divisions per dimension
        - ``v`` - initial vector (default: ``None``)
        - ``norm`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit of the algo
        - ``verbose`` -- bool (default: ``False``)

        OUTPUT:

        dict

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
            sage: D = Brun().invariant_measure_dict(4, 10, verbose=True) # random
            0.0799404500357 0.199341464229 0.720718085735
            0 1
            0.0998433745026 0.248971884172 0.651184741325
            0 2
            0.132942259282 0.331508073966 0.535549666752
            1 3
            0.198868907918 0.495904379777 0.305226712305
            1 4

        ::

            sage: D = Brun().invariant_measure_dict(100000, 5)
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

            sage: D = Brun().invariant_measure_dict(1000000, 10) # 0.05s

        """
        # 146 works, 147 causes segmentation error!!!
        DEF NDIVS = 100
        assert ndivs <= NDIVS, "ndivs(=%s) must be less or equal to %s" % (ndivs, NDIVS)
        cdef double s
        cdef int i,j
        cdef int X,Y

        # initialization of the counter
        # change this to something else
        # see https://groups.google.com/forum/?fromgroups=#!topic/sage-devel/NCBmj2KjwEM
        cpdef int C[NDIVS][NDIVS]
        for j from 0 <= j <= ndivs:
            for i from 0 <= i <= j:
                C[i][j] = 0

        cdef PairPoint3d P
        P.x = random()
        P.y = random()
        P.z = random()
        P.u = .3
        P.v = .3
        P.w = .3
        P.branch = 999

        # Order (x,y,z)
        if P.y > P.z: P.z,P.y = P.y,P.z
        if P.x > P.z: P.x,P.y,P.z = P.y,P.z,P.x
        elif P.x > P.y: P.x,P.y = P.y,P.x

        # Normalize (x,y,z)
        s = P.x + P.y + P.z
        P.x /= s; P.y /= s; P.z /= s

        for i from 0 <= i < n_iterations:

            # Check for Keyboard interupt
            sig_check()

            # Apply Algo
            P = self.call(P)

            # Normalize (xnew,ynew,znew)
            if norm== '1':
                s = P.x + P.y + P.z # norm 1
            elif norm== 'sup':
                s = max(P.x, P.y, P.z) # norm sup
            else:
                raise ValueError("Unknown value for norm(=%s)" %norm)
            P.x /= s; P.y /= s; P.z /= s

            # Increase by one the counter for that part
            X = int(P.x*ndivs)
            Y = int(P.y*ndivs)
            C[X][Y] += 1

            if verbose:
                print P.x,P.y,P.z
                print X,Y

        # Translate the counter into a python dict
        D = {}
        for j from 0 <= j <= ndivs:
            for i from 0 <= i <= j:
                c = C[i][j]
                if c > 0:
                    D[(i,j)] = c
        return D

    def natural_extention_dict(self, int n_iterations, norm_left='sup',
            norm_right='1', verbose=False):
        r"""
        INPUT:

        - ``n_iterations`` -- integer
        - ``norm_left`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit of the algo
        - ``norm_right`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'`` or ``'hypersurfac'``, the norm used for the orbit in the natural extension
        - ``verbose`` -- bool (default: ``False``)

        OUTPUT:

        dict, dict, dict, dict

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ARP
            sage: t = ARP().natural_extention_dict(10000)
            sage: map(type, t)
            [<type 'collections.defaultdict'>,
             <type 'collections.defaultdict'>,
             <type 'collections.defaultdict'>,
             <type 'collections.defaultdict'>]

        """
        cdef double x,y,z           # vector (x,y,z)
        cdef double u,v,w           # vector (u,v,w)
        cdef double p,s,t           # temporary variables
        cdef unsigned int i         # loop counter
        cdef double x_new,y_new,z_new
        cdef double u_new,v_new,w_new

        # random initial values
        x = random(); y = random(); z = random();
        u = 1./3; v = 1./3; w = 1./3;

        # Order (x,y,z)
        if y > z: z,y = y,z
        if x > z: x,y,z = y,z,x
        elif x > y: x,y = y,x

        # Normalize (x,y,z)
        s = x + y + z
        x /= s; y /= s; z /= s

        cdef PairPoint3d P,R
        P.x = x
        P.y = y
        P.z = z
        P.u = u
        P.v = v
        P.w = w

        import collections
        domain_right = collections.defaultdict(list)
        image_right = collections.defaultdict(list)
        domain_left = collections.defaultdict(list)
        image_left = collections.defaultdict(list)

        # Loop
        for i from 0 <= i < n_iterations:

            # Check for Keyboard interupt
            sig_check()

            # Apply Algo
            R = self.call(P)

            if verbose:
                print("x=%f, y=%f, z=%f" % (R.x,R.y,R.z))
                #print("u=%f, v=%f, w=%f" % (u,v,w))
                #s = x*u + y*v + z*w
                #print("scal prod <(x,y,z),(u,v,w)> = %f (after algo)" % s)

            # Normalize (xnew,ynew,znew)
            if norm_left == '1':
                s = R.x + R.y + R.z # norm 1
            elif norm_left == 'sup':
                s = max(R.x, R.y, R.z) # norm sup
            else:
                raise ValueError("Unknown value for norm_left(=%s)" %norm_left)
            R.x /= s; R.y /= s; R.z /= s

            # Normalize (unew,vnew,wnew)
            if norm_right == '1':
                s = R.u + R.v + R.w    # norm 1
            elif norm_right == 'sup':
                s = max(R.u, R.v, R.w) # norm sup
            elif norm_right == 'hypersurface':
                s = R.x*R.u + R.y*R.v + R.z*R.w # hypersurface
            else:
                raise ValueError("Unknown value for norm_right(=%s)" %norm_right)
            R.u /= s; R.v /= s; R.w /= s

            # Projection
            if norm_left == '1':
                s = -SQRT3SUR2 * P.x + SQRT3SUR2 * P.y
                t = -.5 * P.x -.5 * P.y + P.z
                domain_left[R.branch].append((s,t))
                s = -SQRT3SUR2 * R.x + SQRT3SUR2 * R.y
                t = -.5 * R.x -.5 * R.y + R.z
                image_left[R.branch].append((s,t))
            elif norm_left == 'sup':
                domain_left[R.branch].append((P.x,P.y))
                image_left[R.branch].append((R.x,R.y))

            if norm_right == '1':
                s = -SQRT3SUR2 * P.u + SQRT3SUR2 * P.v
                t = -.5 * P.u -.5 * P.v + P.w
                domain_right[R.branch].append((s,t))
                s = -SQRT3SUR2 * R.u + SQRT3SUR2 * R.v
                t = -.5 * R.u -.5 * R.v + R.w
                image_right[R.branch].append((s,t))
            elif norm_right == 'sup':
                domain_right[R.branch].append((P.u,P.v))
                image_right[R.branch].append((R.u,R.v))

            P = R

        return domain_left, image_left, domain_right, image_right

    def lyapounov_exponents(self, int n_iterations=1000, verbose=False):
        r"""
        INPUT:

        - ``n_iterations`` -- integer
        - ``verbose`` -- bool (default: ``False``)

        OUTPUT:

        liapounov exponents

        theta1, theta2, theta2/theta1

        .. NOTE:: 
        
            the code of this method was translated from C to cython. The C
            version is from Vincent Delecroix.

        EXAMPLES:

        Some benchmarks (on my machine)::

            sage: from slabbe.mult_cont_frac import Brun
            sage: Brun().lyapounov_exponents(1000000)  # 68.6 ms # tolerance 0.003
            (0.3049429393152174, -0.1120652699014143, 1.367495867105725)

        Cython code on liafa is as fast as C on my machine::

            sage: Brun().lyapounov_exponents(67000000)    # 3.71s # tolerance 0.001
            (0.30452120021265766, -0.11212586210856369, 1.36820379674801734)

        Cython code on my machine is almost as fast as C on my machine::

            sage: Brun().lyapounov_exponents(67000000) # 4.58 s # tolerance 0.001
            (0.30456433843239084, -0.1121770192467067, 1.36831961293987303)

        """
        from math import log
        cdef double theta1=0, theta2=0    # values of Lyapunov exponents
        cdef double theta1c=0, theta2c=0  # compensation (for Kahan summation algorithm)
        cdef double x,y,z           # vector (x,y,z)
        cdef double u,v,w           # vector (u,v,w)
        cdef double p,s,t           # temporary variables
        cdef unsigned int i         # loop counter
        cdef double critical_value=0.0001

        # random initial values
        x = random(); y = random(); z = random();
        u = random() - .5; v = random() - .5; w = random() - .5;

        # Order (x,y,z)
        if y > z: z,y = y,z
        if x > z: x,y,z = y,z,x
        elif x > y: x,y = y,x

        # Normalize (x,y,z)
        s = x + y + z
        x /= s; y /= s; z /= s

        # Gram Shmidtt on (u,v,w)
        p = x*u + y*v + z*w
        s = x*x + y*y + z*z
        u -= p*x/s; v -= p*y/s; w -= p*z/s

        # Normalize (u,v,w)
        s = abs(u) + abs(v) + abs(w);
        u /= s; v /= s; w /= s

        if verbose:
            print("x=%f, y=%f, z=%f" % (x,y,z))
            print("u=%f, v=%f, w=%f" % (u,v,w))
            s = x*u + y*v + z*w
            print("scal prod <(x,y,z),(u,v,w)> = %f" % s)

        cdef PairPoint3d P
        P.x = x
        P.y = y
        P.z = z
        P.u = u
        P.v = v
        P.w = w

        # Loop
        for i from 0 <= i < n_iterations:

            # Check for Keyboard interupt
            sig_check()

            # Apply Algo
            P = self.call(P)

            if verbose:
                print("x=%f, y=%f, z=%f" % (P.x,P.y,P.z))
                print("u=%f, v=%f, w=%f" % (P.u,P.v,P.w))
                s = P.x*P.u + P.y*P.v + P.z*P.w
                print("scal prod <(x,y,z),(u,v,w)> = %f (after algo)" % s)

            # Save some computations
            #if i % step == 0:
            if P.x < critical_value:

                # Sum the first lyapounov exponent
                s = P.x + P.y + P.z
                p = -log(s) - theta1c
                t = theta1 + p
                theta1c = (t-theta1) - p   # mathematically 0 but not for a computer!!
                theta1 = t
                P.x /= s; P.y /= s; P.z /= s;

                # Sum the second lyapounov exponent
                s = abs(P.u) + abs(P.v) + abs(P.w)
                p = log(s) - theta2c
                t = theta2 + p
                theta2c = (t-theta2) - p   # mathematically 0 but not for a computer!!
                theta2 = t

                # the following gramm shimdts seems to be useless, but it is not!!!
                p = P.x*P.u + P.y*P.v + P.z*P.w
                s = P.x*P.x + P.y*P.y + P.z*P.z
                P.u -= p*P.x/s; P.v -= p*P.y/s; P.w -= p*P.z/s
                s = abs(P.u) + abs(P.v) + abs(P.w)
                P.u /= s; P.v /= s; P.w /= s

        return theta1/n_iterations, theta2/n_iterations, 1-theta2/theta1

    ######################
    # DRAWINGS METHODS (python):
    ######################
    def invariant_measure_plot(self, n_iterations, ndivs, norm='sup',
            savefig=True):
        r"""
        Return a matplotlib graph of the invariant measure.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``ndvis`` - integer, number of divisions per dimension
        - ``norm`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``savefig`` - boolean (default: ``True``)

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Reverse, Brun
            sage: Reverse().invariant_measure_plot(1000000, 80, savefig=False)
            <matplotlib.figure.Figure object at ...>
            sage: Brun().invariant_measure_plot(1000000, 40, norm='sup',savefig=False)
            <matplotlib.figure.Figure object at ...>

        """
        D = self.invariant_measure_dict(n_iterations, ndivs, norm=norm)
        mx,my = map(max, zip(*D.keys()))
        len_D = float(len(D))
        the_mean = n_iterations / len_D

        X = [[i for i in range(mx+1)] for j in range(my+1)]
        Y = [[j for i in range(mx+1)] for j in range(my+1)]
        Z = [[D.get((i,j),0)/the_mean for i in range(mx+1)] for j in range(my+1)]

        title = "Algo=%s, nbiter=%s, ndivs=%s" % (self.name(), n_iterations, ndivs)
        filename = 'mesure_%s_iter%s_div%s.png' % (self.name(), n_iterations, ndivs)

        fig = self._plot_wideframe(X,Y,Z,title)
        if savefig:
            fig.savefig(filename)
            print "Creation of the file %s" % filename
        return fig

    def invariant_measure_inverse_plot(self, n_iterations, ndivs,
            norm='sup', savefig=True):
        r"""
        Return a matplotlib graph of the inverse of the invariant measure.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``ndvis`` - integer, number of divisions per dimension
        - ``norm`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``savefig`` - boolean (default: ``True``)

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Reverse
            sage: Reverse().invariant_measure_inverse_plot(1000000, 80, savefig=False)
            <matplotlib.figure.Figure object at ...>

        """
        D = self.invariant_measure_dict(n_iterations, ndivs, norm=norm)
        mx,my = map(max, zip(*D.keys()))
        len_D = float(len(D))
        the_mean = n_iterations / len_D

        E = dict((key,the_mean/value) for key,value in D.iteritems())

        X = [[i for i in range(mx+1)] for j in range(my+1)]
        Y = [[j for i in range(mx+1)] for j in range(my+1)]
        Z = [[E.get((i,j),0) for i in range(mx+1)] for j in range(my+1)]

        title = "Algo=%s, nbiter=%s, ndivs=%s" % (self.name(), n_iterations, ndivs)
        filename = 'mesure_inversed_%s_iter%s_div%s.png' % (self.name(), n_iterations, ndivs)

        fig = self._plot_wideframe(X,Y,Z,title)
        if savefig:
            fig.savefig(filename)
            print "Creation of the file %s" % filename
        return fig

    def _plot_wideframe(self, X,Y,Z, title):
        r"""
        EXAMPLES::
        """
        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        ax.set_title(title)
        return fig

    def natural_extension_plot(self, n_iterations, norm_left='1',
            norm_right='1', axis_off=False, savefig=True):
        r"""
        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``norm_left`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``norm_right`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``axis_off`` - boolean (default: ``False``), 
        - ``savefig`` - boolean (default: ``True``)

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_ARP
            sage: Sorted_ARP().natural_extension_plot(1000, savefig=False)
            <matplotlib.figure.Figure object at ...>
        """
        import matplotlib
        import matplotlib.pyplot as plt
        t = self.natural_extention_dict(n_iterations, norm_left=norm_left,
                norm_right=norm_right)
        domain_left, image_left, domain_right, image_right = t
        c = dict(zip(domain_left.keys(), ['b','r','g','c','m','y','k']))

        def create_sub_plot(axx, D, title, norm):
            for key, value in D.iteritems():
                value = D[key]
                X,Y = zip(*value)
                A = axx.plot(X, Y, 'o', markersize=2, color=c[key], label=key)
            axx.legend(markerscale=3,fontsize="xx-small")
            axx.set_title(title)

        fig, (ax,bx,cx,dx) = plt.subplots(1,4,sharex=True,sharey=True)
        fig.set_figheight(2)
        fig.set_figwidth(12)
        if axis_off:
            ax.set_axis_off()
            bx.set_axis_off()
            cx.set_axis_off()
            dx.set_axis_off()

        create_sub_plot(ax, domain_left,  "Algo IN",    norm=norm_left)
        create_sub_plot(bx, domain_right,    "NatExt IN",  norm=norm_right)
        create_sub_plot(cx, image_left, "Algo OUT",   norm=norm_left)
        create_sub_plot(dx, image_right,   "NatExt OUT", norm=norm_right)

        title = "Algo=%s, nbiter=%s" % (self.name(), n_iterations)
        fig.suptitle(title)
        if savefig:
            filename = 'nat_ext_%s_iter%s.png' % (self.name(), n_iterations)
            fig.savefig(filename)
            #plt.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)
            print "Creation of the file %s" % filename
        return fig

    def natural_extension_tikz(self, n_iterations, norm_left='1',
            norm_right='1', marksize=0.2, legend_marksize=2, 
            group_size="4 by 1"):
        r"""

        INPUT:

        - ``n_iterations`` -- integer, number of iterations
        - ``norm_left`` -- string (default: ``'1'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``norm_right`` -- string (default: ``'1'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``marksize`` -- tikz marksize (default:``0.2``)
        - ``legend_marksize`` -- tikz legend marksize (default:``2``)
        - ``group_size`` -- string (default:``"4 by 1"``)

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
            sage: s = Brun().natural_extension_tikz(1000)
            sage: view(s, tightpage=True)   # not tested
        """

        t = self.natural_extention_dict(n_iterations, norm_left=norm_left,
                norm_right=norm_right)
        domain_left, image_left, domain_right, image_right = t
        sqrt3 = 1.73205080756888
        r = 1.08
        color_dict = dict(zip(domain_left.keys(), PGF_COLORS))

        lines = []
        lines.append(r"\begin{tikzpicture}[scale=.7]")
        lines.append(r"\begin{groupplot}")
        lines.append(r"[group style={group size=%s}," % group_size)
        lines.append(r"height=7cm,width=8cm,")
        lines.append(r"xmin=-1.1,xmax=1.1,ymin=-.6,ymax=1.20,")
        lines.append(r"hide axis]")
        for data in [domain_left, domain_right, image_left, image_right]:
            lines.append(r"\nextgroupplot")
            lines.append(r"\draw[dashed] ")
            lines.append(r"(axis cs:%s, %s)" % (-r*sqrt3/2,r*-.5))
            lines.append(r" node[left] {$\mathbf{e}_1$} -- ")
            lines.append("(axis cs:%s, %s)" % (r*sqrt3/2,r*-.5))
            lines.append(r" node[right] {$\mathbf{e}_2$} -- ")
            lines.append("(axis cs:%s, %s)" % (0, r))
            lines.append(r" node[above] {$\mathbf{e}_3$} -- cycle;")
            for key,value in data.iteritems():
                lines.append(r"\addplot+[")
                lines.append(r"legend image post style={mark size=%s}," % legend_marksize)
                lines.append(r"only marks,mark=*,")
                lines.append(r"mark size=%s," % marksize)
                lines.append(r"mark options={color=%s}] " % color_dict[key])
                lines.append(r"coordinates {%s};" % '\n'.join(map(str, value)))
                lines.append(r"\addlegendentry{%s}" % key)
        lines.append(r"\end{groupplot}")
        if group_size == "4 by 1":
            lines.append(r"\draw[draw=none] (group c1r1.center) -- ")
            lines.append(r"node {$\times$}  (group c2r1.center);")
            lines.append(r"\draw[draw=none] (group c2r1.center) -- ")
            lines.append(r"node {$\to$}     (group c3r1.center);")
            lines.append(r"\draw[draw=none] (group c3r1.center) -- ")
            lines.append(r"node {$\times$}  (group c4r1.center);")
        lines.append(r"\end{tikzpicture}")
        from sage.misc.latex import LatexExpr
        return LatexExpr('\n'.join(lines))

    def natural_extension_part_tikz(self, n_iterations, part=3, 
                                    norm_left='1', norm_right='1',
                                    marksize='1pt', limit_nb_points=None,
                                    verbose=False):
        r"""
        Return a pgfplots or some part of an orbit in the natural
        extension.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``part`` - integer, taking value 0, 1, 2 or 3
        - ``norm_left`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``norm_right`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``marksize`` -- string (default: ``'1pt'``), pgfplots mark size value
        - ``limit_nb_points`` -- None or integer (default: ``None``), limit
          number of points per patch
        - ``verbose`` -- string (default: ``False``)

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ARP
            sage: s = ARP().natural_extension_part_tikz(1000, part=3)
            sage: view(s, tightpage=True)    # not tested

        also ::

            sage: s = ARP().natural_extension_part_tikz(1000, part=3, marksize='.2pt', limit_nb_points=1200, verbose=True)
            Taking ... points for key 1
            Taking ... points for key 2
            Taking ... points for key 3
            Taking ... points for key 132
            Taking ... points for key 321
            Taking ... points for key 231
            Taking ... points for key 213
            Taking ... points for key 312
            Taking ... points for key 123
        """
        t = self.natural_extention_dict(n_iterations, norm_left=norm_left,
                norm_right=norm_right)
        domain_left, image_left, domain_right, image_right = t
        sqrt3 = 1.73205080756888
        r = 1.08
        color_dict = dict(zip(domain_left.keys(), PGF_COLORS))
        data = t[part]

        s = ''
        s += "\\begin{tikzpicture}[scale=.7]\n"
        s += ("\\begin{axis}[height=7cm,width=8cm,\n"
               "xmin=-1.1,xmax=1.1,ymin=-.6,ymax=1.20,\n"
               "hide axis]\n")
        s += ("\\draw[dashed] \n"
              "(axis cs:%s, %s)" % (-r*sqrt3/2,r*-.5) +
              " node[left] {$\\mathbf{e}_1$} -- \n"
              "(axis cs:%s, %s)" % (r*sqrt3/2,r*-.5)  +
              " node[right] {$\\mathbf{e}_2$} -- \n"
              "(axis cs:%s, %s)" % (0, r)             +
              " node[above] {$\\mathbf{e}_3$} -- cycle;\n")
        for key,value in data.iteritems():
            if limit_nb_points and len(value) > limit_nb_points:
                if verbose:
                    print "Taking only {} points instead of {} for key {}".format(
                            limit_nb_points, len(value), key)
                value = value[:limit_nb_points]
            elif verbose:
                print "Taking {} points for key {}".format(len(value), key)

            s += "\\addplot+[only marks,mark=*,mark options={color=%s}," % color_dict[key]
            s += "mark size=%s]\n" % marksize
            s += "coordinates {\n"
            s += '\n'.join(map(str, value))
            s += "};\n" 
            s += "\\addlegendentry{%s}\n " % key
        s += "\\end{axis}\n"
        s += "\\end{tikzpicture}\n"
        from sage.misc.latex import LatexExpr
        return LatexExpr(s)

    def natural_extension_part_png(self, n_iterations, draw,
                                    norm_left='1', norm_right='1',
                                    xrange=(-.866, .866),
                                    yrange=(-.5, 1.),
                                    urange=(-.866, .866),
                                    vrange=(-.5, 1.),
                                    color_dict=None,
                                    branch_order=None,
                                    ndivs=1024,
                                    verbose=False):
        r"""
        Return a png or some part of an orbit in the natural extension.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``draw`` -- string (default: ``'image_right'``), possible values
          are:

          - ``'domain_left'`` - use x and y ranges
          - ``'domain_right'`` - use u and v ranges
          - ``'image_left'`` - use x and y ranges
          - ``'image_right'`` - use u and v ranges

        - ``norm_left`` -- string (default: ``'1'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``norm_right`` -- string (default: ``'1'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``xrange`` -- tuple (default: ``(-.866, .866)``), interval of
          values for x
        - ``yrange`` -- tuple (default: ``(-.5, 1.)``), interval of
          values for y
        - ``urange`` -- tuple (default: ``(-.866, .866)``), interval of
          values for u
        - ``vrange`` -- tuple (default: ``(-.5, 1.)``), interval of
          values for v
        - ``color_dict`` -- dict (default: ``None``), dict from branches
          int to color (as RGB tuples)
        - ``branch_order`` -- list (default: ``None``), list of branches int
        - ``ndivs`` -- int (default: ``1024``), number of pixels
        - ``verbose`` -- string (default: ``False``)

        BENCHMARK:

        A minute (1min 13s) for a picture with 10^7 points::

            sage: from slabbe.mult_cont_frac import ARP
            sage: c = {}
            sage: c[1] = c[2] = c[3] = [0,0,0]
            sage: c[12] = c[13] = c[23] = c[21] = c[31] = c[32] = [255,0,0]
            sage: b = [1,2,3,12,13,21,23,31,32]
            sage: P = ARP().natural_extension_part_png(10^7, draw='image_right',  # not tested
            ....:   branch_order=b, color_dict=c, urange=(-.6,.6), vrange=(-.6,.6))  # not tested

        Half a minute (27s) for a picture zoomed in the orbit of 10^8
        points::

            sage: P = ARP().natural_extension_part_png(10^8, draw='image_right', # not tested
            ....:   branch_order=b, color_dict=c, urange=(.2,.3), vrange=(.2,.3))   # not tested

        EXAMPLES::

            sage: c = {}
            sage: c[1] = c[2] = c[3] = [0,0,0]
            sage: c[123] = c[132] = c[231] = c[213] = c[312] = c[321] = [255,0,0]
            sage: b = [1,2,3,123,132,213,231,312,321]
            sage: opt = dict(urange=(-.6,.6), vrange=(-.6,.6), color_dict=c, branch_order=b)
            sage: P = ARP().natural_extension_part_png(10^5, draw='domain_left', **opt)
            sage: P = ARP().natural_extension_part_png(10^5, draw='domain_right', **opt)
            sage: P = ARP().natural_extension_part_png(10^5, draw='image_left', **opt)
            sage: P = ARP().natural_extension_part_png(10^5, draw='image_right', **opt)
            sage: P.show() # not tested

        """
        L = self.orbit_filtered_list(n_iterations,
            norm_left=norm_left, norm_right=norm_right,
            xmin=xrange[0], xmax=xrange[1],
            ymin=yrange[0], ymax=yrange[1],
            umin=urange[0], umax=urange[1],
            vmin=vrange[0], vmax=vrange[1],
            ndivs=ndivs)
        if branch_order is None:
            branch_order = []
            raise NotImplementedError
        if color_dict is None:
            from random import randint
            color_dict = {}
            for key in branch_order:
                color_dict[key] = [randint(0,255),randint(0,255),randint(0,255)]

        #http://stackoverflow.com/questions/434583/what-is-the-fastest-way-to-draw-an-image-from-discrete-pixel-values-in-python
        import numpy as np
        import scipy.misc as smp

        # Create a 1024x1024x3 array of 8 bit unsigned integers
        data = np.zeros( (ndivs,ndivs,3), dtype=np.uint8 )
        data += 255   # white as default color

        if draw.startswith('domain'):
            L.sort(key=lambda a:branch_order.index(a[5]))
        elif draw.startswith('image'):
            L.sort(key=lambda a:branch_order.index(a[4]))
        else:
            raise ValueError("Unkown value for draw(={})".format(draw))

        if draw == 'domain_left':
            for x,y,u,v,prev_br,next_br in L:
                data[y,x] = color_dict[next_br]
        elif draw == 'domain_right':
            for x,y,u,v,prev_br,next_br in L:
                data[v,u] = color_dict[next_br]
        elif draw == 'image_left':
            for x,y,u,v,prev_br,next_br in L:
                data[y,x] = color_dict[prev_br]
        elif draw == 'image_right':
            for x,y,u,v,prev_br,next_br in L:
                data[v,u] = color_dict[prev_br]
        else:
            raise ValueError("Unkown value for draw(={})".format(draw))

        img = smp.toimage( data )       # Create a PIL image
        #img.show()                      # View in default viewer
        return img

    def measure_evaluation(self, n_iterations, draw,
                                norm_left='1', norm_right='1',
                                xrange=(-.866, .866),
                                yrange=(-.5, 1.),
                                urange=(-.866, .866),
                                vrange=(-.5, 1.),
                                ndivs=1024,
                                verbose=False):
        r"""

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``draw`` -- string (default: ``'image_right'``), possible values
          are:

          - ``'domain_left'`` - use x and y ranges
          - ``'domain_right'`` - use u and v ranges
          - ``'image_left'`` - use x and y ranges
          - ``'image_right'`` - use u and v ranges

        - ``norm_left`` -- string (default: ``'1'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``norm_right`` -- string (default: ``'1'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``xrange`` -- tuple (default: ``(-.866, .866)``), interval of
          values for x
        - ``yrange`` -- tuple (default: ``(-.5, 1.)``), interval of
          values for y
        - ``urange`` -- tuple (default: ``(-.866, .866)``), interval of
          values for u
        - ``vrange`` -- tuple (default: ``(-.5, 1.)``), interval of
          values for v
        - ``ndivs`` -- int (default: ``1024``), number of pixels
        - ``verbose`` -- string (default: ``False``)

        BENCHMARK::

            sage: from slabbe.mult_cont_frac import ARP
            sage: opt = dict(urange=(-.15,.25), vrange=(-.05,.05))
            sage: ARP().measure_evaluation(10^8, draw='right', ndivs=100, **opt) # optional long
            0.435...
            sage: ARP().measure_evaluation(10^8, draw='right', ndivs=1000, **opt) # optional long
            0.357...
            sage: ARP().measure_evaluation(10^8, draw='right', ndivs=2000, **opt) # optional long
            0.293...
            sage: ARP().measure_evaluation(10^8, draw='right', ndivs=4000, **opt) # optional long
            0.177...

        """
        L = self.orbit_filtered_list(n_iterations,
            norm_left=norm_left, norm_right=norm_right,
            xmin=xrange[0], xmax=xrange[1],
            ymin=yrange[0], ymax=yrange[1],
            umin=urange[0], umax=urange[1],
            vmin=vrange[0], vmax=vrange[1],
            ndivs=ndivs)

        if draw.endswith('left'):
            S = set((p[0],p[1]) for p in L)
        elif draw.endswith('right'):
            S = set((p[2],p[3]) for p in L)
        else:
            raise ValueError("Unkown value for draw(={})".format(draw))

        if verbose:
            print "nombre diterations dans la fenetre : ", len(L)
            print "{} pixels touchés parmi limage {}^2 ".format(len(S), ndivs)

        print float(len(S) / ndivs**2)

    def lyapounov_exponents_sample(self, ntimes, n_iterations=1000):
        r"""
        Return two lists of values for theta1 and theta2

        INPUT:

        - ``ntimes`` -- integer, number of orbits
        - ``n_iterations`` -- integer, length of each orbit

        OUTPUT:

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
            sage: T1, T2, U = Brun().lyapounov_exponents_sample(10, 100000)
            sage: T1[0]
            0.303680940345907
            sage: T2[0]
            -0.1119022164698561 
            sage: U[0]
            1.3684861366090153

        """
        Theta1 = []
        Theta2 = []
        Uniform = []
        for _ in range(ntimes):
            l1,l2,approx = self.lyapounov_exponents(n_iterations)
            Theta1.append(l1)
            Theta2.append(l2)
            Uniform.append(approx)
        return Theta1, Theta2, Uniform

    def lyapounov_exponents_table(self, ntimes, n_iterations=1000):
        r"""
        Return a table of values of 1 - theta2/theta1.

        INPUT:

        - ``ntimes`` -- integer, number of orbits
        - ``n_iterations`` -- integer, length of each orbit

        OUTPUT:

        liapounov exponents

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_Brun
            sage: Sorted_Brun().lyapounov_exponents_table(10, 100000) # abs tol 0.005
            n          | 10
            iterations | 100000
            min        | 1.364250082762
            mean       | 1.367620100939
            max        | 1.369748622574
            std        | 0.00155230985651

        """
        T1, T2, U = self.lyapounov_exponents_sample(ntimes, n_iterations)
        import numpy as np
        a = np.array(U)
        header = ['n', 'iterations', 'min','mean','max','std']
        cols = (ntimes, n_iterations, a.min(), a.mean(), a.max(), a.std()),
        from sage.misc.table import table
        t = table(columns=cols,header_column=header)
        return t

cdef class Brun(MCFAlgorithm):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Brun
        sage: algo = Brun()
        sage: TestSuite(algo).run()
    """
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
            sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Brun()(D)
            sage: sorted(E.iteritems())
            [('branch', 123),
             ('u', 0.2),
             ('v', 0.6),
             ('w', 0.3),
             ('x', 0.3),
             ('y', 0.6),
             ('z', 0.20000000000000007)]
        """
        if P.x < P.y < P.z:
            P.z -= P.y
            P.v += P.w
            P.branch = 123
        elif P.x < P.z < P.y:
            P.y -= P.z
            P.w += P.v
            P.branch = 132
        elif P.y < P.z < P.x:
            P.x -= P.z
            P.w += P.u
            P.branch = 231
        elif P.y < P.x < P.z:
            P.z -= P.x
            P.u += P.w
            P.branch = 213
        elif P.z < P.x < P.y:
            P.y -= P.x
            P.u += P.v
            P.branch = 312
        elif P.z < P.y < P.x:
            P.x -= P.y
            P.v += P.u
            P.branch = 321
        else:
            raise ValueError('limit case: reach set of measure zero')
        return P

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
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

            sage: from slabbe.mult_cont_frac import Brun
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
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Reverse
            sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Reverse()(D)
            sage: sorted(E.iteritems())
            [('branch', 4),
             ('u', 0.47622031559045996),
             ('v', 0.39685026299205),
             ('w', 0.39685026299205),
             ('x', 0.6929565774421808),
             ('y', 0.3149802624737185),
             ('z', 0.06299605249474362)]
        """
        cdef PairPoint3d R
        if P.x + P.y < P.z:
            P.z -= P.x + P.y
            P.v += P.w
            P.u += P.w
            P.branch = 3
            return P
        elif P.x + P.z < P.y:
            P.y -= P.x + P.z
            P.w += P.v
            P.u += P.v
            P.branch = 2
            return P
        elif P.y + P.z < P.x:
            P.x -= P.y + P.z
            P.v += P.u
            P.w += P.u
            P.branch = 1
            return P
        else:
            R.x = 0.629960524947437 * (-P.x + P.y + P.z)
            R.y = 0.629960524947437 * ( P.x - P.y + P.z)
            R.z = 0.629960524947437 * ( P.x + P.y - P.z)
            # 0.793700525984100 = 1/2*4^(1/3)
            # 0.629960524947437 = 1/4*4^(2/3)
            R.u = 0.793700525984100 * (P.v + P.w)
            R.v = 0.793700525984100 * (P.u + P.w)
            R.w = 0.793700525984100 * (P.u + P.v)
            R.branch = 4
            return R

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Reverse
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

            sage: from slabbe.mult_cont_frac import ArnouxRauzy
            sage: ArnouxRauzy().dual_substitutions()
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
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ARP
            sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = ARP()(D)
            sage: sorted(E.iteritems())
            [('branch', 123),
             ('u', 0.8),
             ('v', 0.6),
             ('w', 0.3),
             ('x', 0.3),
             ('y', 0.3),
             ('z', 0.20000000000000007)]
        """
        if P.x + P.y < P.z:
            P.z -= P.x + P.y
            P.v += P.w
            P.u += P.w
            P.branch = 3
            return P
        elif P.x + P.z < P.y:
            P.y -= P.x + P.z
            P.w += P.v
            P.u += P.v
            P.branch = 2
            return P
        elif P.y + P.z < P.x:
            P.x -= P.y + P.z
            P.v += P.u
            P.w += P.u
            P.branch = 1
            return P
        else:
            return _Poincare(P)

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ARP
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

            sage: from slabbe.mult_cont_frac import ARP
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
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ArnouxRauzy
            sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = ArnouxRauzy()(D)
            Error: arnoux rauzy not defined on input:
            {'branch': 999, 'u': 0.2, 'w': 0.3, 'v': 0.3, 'y': 0.6, 'x': 0.3, 'z': 0.8}

        :: 

            sage: D = {'x':.3,'y':.2,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = ArnouxRauzy()(D)
            sage: sorted(E.iteritems())
            [('branch', 3),
             ('u', 0.5),
             ('v', 0.6),
             ('w', 0.3),
             ('x', 0.3),
             ('y', 0.2),
             ('z', 0.30000000000000004)]
        """
        if P.x + P.y < P.z:
            P.z -= P.x + P.y
            P.v += P.w
            P.u += P.w
            P.branch = 3
            return P
        elif P.x + P.z < P.y:
            P.y -= P.x + P.z
            P.w += P.v
            P.u += P.v
            P.branch = 2
            return P
        elif P.y + P.z < P.x:
            P.x -= P.y + P.z
            P.v += P.u
            P.w += P.u
            P.branch = 1
            return P
        else:
            print ("Error: arnoux rauzy not defined on input:\n%s"%P)
    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ArnouxRauzy
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

            sage: from slabbe.mult_cont_frac import ArnouxRauzy
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
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Poincare
            sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Poincare()(D)
            sage: sorted(E.iteritems())
            [('branch', 123),
             ('u', 0.8),
             ('v', 0.6),
             ('w', 0.3),
             ('x', 0.3),
             ('y', 0.3),
             ('z', 0.20000000000000007)]
        """
        return _Poincare(P)

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Poincare
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

            sage: from slabbe.mult_cont_frac import Poincare
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
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Selmer
            sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Selmer()(D)
            sage: sorted(E.iteritems())
            [('branch', 123),
             ('u', 0.5),
             ('v', 0.3),
             ('w', 0.3),
             ('x', 0.3),
             ('y', 0.6),
             ('z', 0.5)]
        """
        if P.x < P.y < P.z:
            P.z -= P.x
            P.u += P.w
            P.branch = 123
        elif P.x < P.z < P.y:
            P.y -= P.x
            P.u += P.v
            P.branch = 132
        elif P.y < P.z < P.x:
            P.x -= P.y
            P.v += P.u
            P.branch = 231
        elif P.y < P.x < P.z:
            P.z -= P.y
            P.v += P.w
            P.branch = 213
        elif P.z < P.x < P.y:
            P.y -= P.z
            P.w += P.v
            P.branch = 312
        elif P.z < P.y < P.x:
            P.x -= P.z
            P.w += P.u
            P.branch = 321
        else:
            raise ValueError('limit case: reach set of measure zero')
        return P

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Selmer
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

            sage: from slabbe.mult_cont_frac import Selmer
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

cdef class Meester(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Meester
            sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Meester()(D)
            sage: sorted(E.iteritems())
            [('branch', 1),
             ('u', 0.8),
             ('v', 0.3),
             ('w', 0.3),
             ('x', 0.3),
             ('y', 0.3),
             ('z', 0.5)]
        """
        if P.x <= P.y and P.x <= P.z:
            P.y -= P.x
            P.z -= P.x
            P.u += P.v + P.w
            P.branch = 1
        elif P.y <= P.x and P.x <= P.z:
            P.x -= P.y
            P.z -= P.y
            P.v += P.u + P.w
            P.branch = 2
        elif P.z <= P.x and P.z <= P.y:
            P.x -= P.z
            P.y -= P.z
            P.w += P.u + P.v
            P.branch = 3
        else:
            raise ValueError('limit case: reach set of measure zero')
        return P

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Meester
            sage: Meester().substitutions()
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

            sage: from slabbe.mult_cont_frac import Meester
            sage: Meester().dual_substitutions()
            {1: WordMorphism: 1->1, 2->21, 3->31,
             2: WordMorphism: 1->12, 2->2, 3->32,
             3: WordMorphism: 1->13, 2->23, 3->3}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {1:  WordMorphism({1: [1], 2: [2, 1], 3: [3, 1]}),
                2:  WordMorphism({1: [1, 2], 2: [2], 3: [3, 2]}),
                3:  WordMorphism({1: [1, 3], 2: [2, 3], 3: [3]})}

cdef class Cassaigne(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        This algorithm was provided by Julien Cassaigne during a meeting of
        the ANR DynA3S on October 12, 2015 held in Paris. It is inspired
        from a method to generate words of complexity 2n+1 on a three
        letter alphabet of arbitrary letter frequencies.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Cassaigne
            sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Cassaigne()(D)
            sage: sorted(E.iteritems())
            [('branch', 2),
             ('u', 0.3),
             ('v', 0.5),
             ('w', 0.3),
             ('x', 0.6),
             ('y', 0.3),
             ('z', 0.5)]
        """
        cdef double tmp
        if P.x > P.z :
            P.x -= P.z
            tmp = P.y
            P.y = P.z
            P.z = tmp
            tmp = P.v
            P.v = P.u + P.w
            P.w = tmp
            P.branch = 1
        elif P.x < P.z :
            P.z -= P.x
            tmp = P.y
            P.y = P.x
            P.x = tmp
            tmp = P.v
            P.v = P.u + P.w
            P.u = tmp
            P.branch = 2
        else:
            raise ValueError('limit case: reach set of measure zero')
        return P

    def substitutions(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Cassaigne
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

            sage: from slabbe.mult_cont_frac import Cassaigne
            sage: Cassaigne().dual_substitutions()
            {1: WordMorphism: 1->12, 2->3, 3->2, 
             2: WordMorphism: 1->2, 2->1, 3->23}
        """
        from sage.combinat.words.morphism import WordMorphism
        return {1:  WordMorphism({1: [1,2], 2: [3], 3: [2]}),
                2:  WordMorphism({1: [2], 2: [1], 3: [2,3]})}

cdef class Sorted_Brun(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_Brun
            sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Sorted_Brun()(D)
            sage: sorted(E.iteritems())
            [('branch', 106),
             ('u', 0.3),
             ('v', 0.2),
             ('w', 0.6),
             ('x', 0.20000000000000007),
             ('y', 0.3),
             ('z', 0.6)]
            sage: D = {'x':.3,'y':.45,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Sorted_Brun()(D)
            sage: sorted(E.iteritems())
            [('branch', 102),
             ('u', 0.2),
             ('v', 0.3),
             ('w', 0.6),
             ('x', 0.3),
             ('y', 0.35000000000000003),
             ('z', 0.45)]
            sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Sorted_Brun()(D)
            sage: sorted(E.iteritems())
            [('branch', 100),
             ('u', 0.2),
             ('v', 0.6),
             ('w', 0.3),
             ('x', 0.3),
             ('y', 0.3),
             ('z', 0.5)]

        """
        P.z -= P.y
        P.v += P.w
        P.branch = 100
        return Sort(P)

cdef class Sorted_BrunMulti(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_BrunMulti
            sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Sorted_BrunMulti()(D)
            sage: sorted(E.iteritems())
            [('branch', 106),
             ('u', 0.3),
             ('v', 0.2),
             ('w', 0.8999999999999999),
             ('x', 0.20000000000000007),
             ('y', 0.3),
             ('z', 0.3)]
        """
        cdef int m = <int>(P.z / P.y)
        P.z -= m*P.y
        P.v += m*P.w
        P.branch = 100
        return Sort(P)

cdef class Sorted_Selmer(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_Selmer
            sage: D = {'x':.2,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Sorted_Selmer()(D)
            sage: sorted(E.iteritems())
            [('branch', 100),
             ('u', 0.5),
             ('v', 0.3),
             ('w', 0.3),
             ('x', 0.2),
             ('y', 0.3),
             ('z', 0.6000000000000001)]
        """
        P.z -= P.x
        P.u += P.w
        P.branch = 100
        return Sort(P)

cdef class Sorted_Meester(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_Meester
            sage: D = {'x':.5,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Sorted_Meester()(D)
            sage: sorted(E.iteritems())
            [('branch', 103),
             ('u', 0.3),
             ('v', 0.3),
             ('w', 0.8),
             ('x', 0.09999999999999998),
             ('y', 0.30000000000000004),
             ('z', 0.5)]
        """
        # Apply the algo
        P.y -= P.x
        P.z -= P.x
        P.u += P.v + P.w
        P.branch = 100
        return Sort(P)

cdef class Sorted_ArnouxRauzy(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_ArnouxRauzy
            sage: D = dict(x=.3,y=.4,z=.8,u=.2,v=.3,w=.4,branch=999)
            sage: E = Sorted_ArnouxRauzy()(D)
            sage: sorted(E.iteritems())
            [('branch', 106),
             ('u', 0.4),
             ('v', 0.6000000000000001),
             ('w', 0.7),
             ('x', 0.10000000000000009),
             ('y', 0.3),
             ('z', 0.4)]

        ::

            sage: D = dict(x=.3,y=.7,z=.8,u=.2,v=.3,w=.4,branch=999)
            sage: E = Sorted_ArnouxRauzy()(D)
            sage: sorted(E.iteritems())
            [('branch', 106),
             ('u', 0.4),
             ('v', 0.6000000000000001),
             ('w', 0.7),
             ('x', -0.19999999999999996),
             ('y', 0.3),
             ('z', 0.7)]

        """
        #Arnoux-Rauzy
        cdef PairPoint3d R
        R.x = P.x
        R.y = P.y
        R.z = P.z - (P.x + P.y)
        R.u = P.u + P.w
        R.v = P.v + P.w
        R.w = P.w
        R.branch = 100
        return Sort(R)

cdef class Sorted_ArnouxRauzyMulti(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

        """
        #Arnoux-Rauzy Multi
        cdef int m
        m = <int>(P.z / (P.x + P.y))
        P.z -= m * (P.x + P.y)
        P.v += m * P.w;
        P.u += m * P.w;
        P.branch = 100
        return Sort(P)

cdef class Sorted_ARP(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_ARP
            sage: D = dict(x=.3,y=.4,z=.8,u=.2,v=.3,w=.4,branch=999)
            sage: E = Sorted_ARP()(D)
            sage: sorted(E.iteritems())
            [('branch', 106),
             ('u', 0.4),
             ('v', 0.6000000000000001),
             ('w', 0.7),
             ('x', 0.10000000000000009),
             ('y', 0.3),
             ('z', 0.4)]

        ::

            sage: D = dict(x=.3,y=.7,z=.8,u=.2,v=.3,w=.4,branch=999)
            sage: E = Sorted_ARP()(D)
            sage: sorted(E.iteritems())
            [('branch', 206),
             ('u', 0.4),
             ('v', 0.9),
             ('w', 0.7),
             ('x', 0.10000000000000009),
             ('y', 0.3),
             ('z', 0.39999999999999997)]
        """
        # Apply the algo
        if P.z > P.x + P.y:
            return _Sorted_ArnouxRauzy(P)
        else:
            return _Sorted_Poincare(P)

cdef class Sorted_ARPMulti(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_ARPMulti
            sage: D = {'x':.3,'y':.5,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Sorted_ARPMulti()(D)
            sage: sorted(E.iteritems())
            [('branch', 201),
             ('u', 0.6),
             ('v', 0.8),
             ('w', 0.3),
             ('x', 0.2),
             ('y', 0.3),
             ('z', 0.30000000000000004)]
        """
        # Apply the algo
        if P.z > P.x + P.y:
            return _Sorted_ArnouxRauzyMulti(P)
        else:
            return _Sorted_Poincare(P)

cdef class Sorted_Poincare(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_Poincare
            sage: D = dict(x=.3,y=.7,z=.8,u=.2,v=.3,w=.4,branch=999)
            sage: E = Sorted_Poincare()(D)
            sage: sorted(E.iteritems())
            [('branch', 206),
             ('u', 0.4),
             ('v', 0.9),
             ('w', 0.7),
             ('x', 0.10000000000000009),
             ('y', 0.3),
             ('z', 0.39999999999999997)]

        ::

            sage: E = Sorted_Poincare()(D)
            sage: sorted(E.iteritems())
            [('branch', 206),
             ('u', 0.4),
             ('v', 0.9),
             ('w', 0.7),
             ('x', 0.10000000000000009),
             ('y', 0.3),
             ('z', 0.39999999999999997)]

        """
        # Apply the algo
        cdef PairPoint3d R
        R.x = P.x
        R.y = P.y - P.x
        R.z = P.z - P.y
        R.u = P.u + P.v + P.w
        R.v = P.v + P.w
        R.w = P.w
        R.branch = 200
        return Sort(R)

cdef class Sorted_ARrevert(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_ARrevert
            sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Sorted_ARrevert()(D)
            sage: sorted(E.iteritems())
            [('branch', 300),
             ('u', 0.3),
             ('v', 0.5),
             ('w', 0.6),
             ('x', 0.2),
             ('y', 0.3),
             ('z', 0.3)]
        """
        cdef PairPoint3d R
        cdef double z_x_y = P.z - P.x - P.y
        if z_x_y > P.y:
            R.x = P.x
            R.y = P.y
            R.z = z_x_y
            R.u = P.u + P.w
            R.v = P.v + P.w
            R.w = P.w
            R.branch = 100
        elif z_x_y > P.x:
            R.x = P.x
            R.y = z_x_y
            R.z = P.y
            R.u = P.u + P.w
            R.v = P.w
            R.w = P.v + P.w
            R.branch = 200
        elif z_x_y > 0:
            R.x = z_x_y
            R.y = P.x
            R.z = P.y
            R.u = P.w
            R.v = P.u + P.w
            R.w = P.v + P.w
            R.branch = 300
        else:
            # Revert
            R.x = (P.x + P.y - P.z)/2
            R.y = (P.x - P.y + P.z)/2
            R.z = (-P.x + P.y + P.z)/2
            R.u = P.u + P.v
            R.v = P.u + P.w
            R.w = P.v + P.w
            R.branch = 400
        return R

cdef class Sorted_ARrevertMulti(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_ARrevertMulti
            sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Sorted_ARrevertMulti()(D)
            sage: sorted(E.iteritems())
            [('branch', 400),
             ('u', 0.3),
             ('v', 0.5),
             ('w', 0.6),
             ('x', 0.20000000000000007),
             ('y', 0.3),
             ('z', 0.3)]
        """
        cdef PairPoint3d R
        cdef int m = <int>(P.z / (P.x + P.y))
        cdef double z_mxy = P.z - m * (P.x + P.y)
        if m == 0:
            # Revert
            R.x = (P.x + P.y - P.z)/2
            R.y = (P.x - P.y + P.z)/2
            R.z = (-P.x + P.y + P.z)/2
            R.u = P.u + P.v
            R.v = P.u + P.w
            R.w = P.v + P.w
            R.branch = 100
        elif z_mxy > P.y:
            R.x = P.x
            R.y = P.y
            R.z = z_mxy
            R.u = P.u + m*P.w
            R.v = P.v + m*P.w
            R.w = P.w
            R.branch = 200
        elif z_mxy > P.x:
            R.x = P.x
            R.y = z_mxy
            R.z = P.y
            R.u = P.u + m*P.w
            R.v = P.w
            R.w = P.v + m*P.w
            R.branch = 300
        else:
            R.x = z_mxy
            R.y = P.x
            R.z = P.y
            R.u = P.w
            R.v = P.u + m*P.w
            R.w = P.v + m*P.w
            R.branch = 400
        return R

cdef class Sorted_ARMonteil(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_ARMonteil
            sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Sorted_ARMonteil()(D)
            sage: sorted(E.iteritems())
            [('branch', 106),
             ('u', 0.3),
             ('v', 0.5),
             ('w', 0.6),
             ('x', 0.2),
             ('y', 0.3),
             ('z', 0.3)]
        """
        cdef PairPoint3d R

        # Apply the algo
        if P.z > P.x + P.y:
            # Arnoux-Rauzy
            R.z = P.z - P.y - P.x
            R.x = P.x
            R.y = P.y
            R.u = P.u + P.w
            R.v = P.v + P.w
            R.w = P.w
            R.branch = 100
        else:
            # Monteil
            R.x = P.x + P.y - P.z
            R.y = -P.x + P.z
            R.z = -P.y + P.z
            R.u = P.u + P.v + P.w
            R.v = P.v + P.w
            R.w = P.u + P.w
            R.branch = 200

        return Sort(R)
cdef class Sorted_Delaunay(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        Donné par Xavier Provençal (inspiré de en fait) le 3 février 2014.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Sorted_Delaunay
            sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = Sorted_Delaunay()(D)
            sage: sorted(E.iteritems())
            [('branch', 201),
             ('u', 0.6),
             ('v', 0.5),
             ('w', 0.3),
             ('x', 0.0),
             ('y', 0.3),
             ('z', 0.8)]
        """
        cdef PairPoint3d R
        # Apply the algo
        if P.z > P.x + P.y:
            # Genre de semi revert
            R.x = P.x
            R.y = P.y - P.x
            R.z = P.x - P.y + P.z
            R.u = P.u + P.v
            R.v = P.v + P.w
            R.w = P.w
            R.branch = 200
            return Sort(R)
        else:
            return _Sorted_Meester(P)
cdef class JacobiPerron(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import JacobiPerron
            sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = JacobiPerron()(D)
            sage: sorted(E.iteritems())
            [('branch', 100),
             ('u', 0.3),
             ('v', 0.3),
             ('w', 1.0999999999999999),
             ('x', 0.0),
             ('y', 0.20000000000000007),
             ('z', 0.3)]
        """
        cdef PairPoint3d R
        cdef int m,n                # temporary integer variables
        cdef double r,s,t           # temporary variables

        R.branch = 100

        # Apply the algo
        m = int(P.z / P.x)
        n = int(P.y / P.x)
        t = P.z - m*P.x
        s = P.y - n*P.x
        r = P.x

        R.z = r
        R.y = t
        R.x = s

        t = P.w
        s = P.v
        r = m*P.w + n*P.v + P.u

        R.w = r
        R.v = t
        R.u = s

        return R
cdef class JacobiPerronAdditif(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import JacobiPerronAdditif
            sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = JacobiPerronAdditif()(D)
            sage: sorted(E.iteritems())
            [('branch', 200),
             ('u', 0.5),
             ('v', 0.3),
             ('w', 0.3),
             ('x', 0.3),
             ('y', 0.3),
             ('z', 0.5)]
        """
        cdef PairPoint3d R
        # Apply the algo
        if P.x < P.y:
            R.branch = 100
            R.x = P.x
            R.y = P.y - P.x
            R.z = P.z
            R.u = P.u + P.v
            R.v = P.v
            R.w = P.w
        elif P.x < P.z:
            R.branch = 200
            R.x = P.x
            R.y = P.y
            R.z = P.z - P.x
            R.u = P.u + P.w
            R.v = P.v
            R.w = P.w
        elif P.x > P.y and P.x > P.z:
            R.branch = 300
            R.x = P.y
            R.y = P.z
            R.z = P.x
            R.u = P.v
            R.v = P.w
            R.w = P.u
        else:
            raise ValueError("jacobi not defined for (x,y,z)=(%s,%s,%s)"%(P.x,P.y,P.z))
        return R
cdef class JacobiPerronAdditifv2(MCFAlgorithm):
    cdef PairPoint3d call(self, PairPoint3d P) except *:
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import JacobiPerronAdditifv2
            sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
            sage: E = JacobiPerronAdditifv2()(D)
            sage: sorted(E.iteritems())
            [('branch', 200),
             ('u', 0.5),
             ('v', 0.3),
             ('w', 0.3),
             ('x', 0.3),
             ('y', 0.3),
             ('z', 0.5)]
        """
        cdef PairPoint3d R
        # Apply the algo
        if P.x < P.y:
            R.branch = 100
            R.x = P.x
            R.y = P.y - P.x
            R.z = P.z - P.x
            R.u = P.u + P.v + P.w
            R.v = P.v
            R.w = P.w
        elif P.x < P.z:
            R.branch = 200
            R.x = P.x
            R.y = P.y
            R.z = P.z - P.x
            R.u = P.u + P.w
            R.v = P.v
            R.w = P.w
        elif P.x > P.y and P.x > P.z:
            R.branch = 300
            R.x = P.y
            R.y = P.z
            R.z = P.x
            R.u = P.v
            R.v = P.w
            R.w = P.u
        else:
            raise ValueError("jacobi not defined for (x,y,z)=(%s,%s,%s)"%(P.x,P.y,P.z))

        return R

cdef inline PairPoint3d _Poincare(PairPoint3d P) except *:
    r"""
    EXAMPLES::

    """
    cdef PairPoint3d R
    if P.x <= P.y <= P.z:
        R.x = P.x
        R.y = P.y - P.x
        R.z = P.z - P.y
        R.u = P.u + P.v + P.w
        R.v = P.v + P.w
        R.w = P.w
        R.branch = 123
    elif P.x <= P.z <= P.y:
        R.x = P.x
        R.y = P.y - P.z
        R.z = P.z - P.x
        R.u = P.u + P.v + P.w
        R.v = P.v
        R.w = P.v + P.w
        R.branch = 132
    elif P.y <= P.x <= P.z:
        R.x = P.x - P.y
        R.y = P.y
        R.z = P.z - P.x
        R.u = P.u       + P.w
        R.v = P.u + P.v + P.w
        R.w = P.w
        R.branch = 213
    elif P.z <= P.x <= P.y:
        R.x = P.x - P.z
        R.y = P.y - P.x
        R.z = P.z
        R.u = P.u + P.v
        R.v = P.v
        R.w = P.u + P.v + P.w
        R.branch = 312
    elif P.y <= P.z <= P.x:
        R.x = P.x - P.z
        R.y = P.y
        R.z = P.z - P.y
        R.u = P.u
        R.v = P.u + P.v + P.w
        R.w = P.u + P.w
        R.branch = 231
    elif P.z <= P.y <= P.x:
        R.x = P.x - P.y
        R.y = P.y - P.z
        R.z = P.z
        R.u = P.u
        R.v = P.u + P.v
        R.w = P.u + P.v + P.w
        R.branch = 321
    else:
        raise ValueError('limit case: reach set of measure zero')
    return R

cdef inline PairPoint3d _Sorted_ArnouxRauzy(PairPoint3d P):
    r"""
    EXAMPLES::

    """
    #Arnoux-Rauzy
    cdef PairPoint3d R
    R.x = P.x
    R.y = P.y
    R.z = P.z - (P.x + P.y)
    R.u = P.u + P.w
    R.v = P.v + P.w
    R.w = P.w
    R.branch = 100
    return Sort(R)

cdef inline PairPoint3d _Sorted_ArnouxRauzyMulti(PairPoint3d P):
    r"""
    EXAMPLES::

    """
    #Arnoux-Rauzy Multi
    cdef int m
    m = <int>(P.z / (P.x + P.y))
    P.z -= m * (P.x + P.y)
    P.v += m * P.w;
    P.u += m * P.w;
    P.branch = 100
    return Sort(P)

cdef inline PairPoint3d _Sorted_Poincare(PairPoint3d P):
    r"""
    EXAMPLES::

    """
    # Apply the algo
    cdef PairPoint3d R
    R.x = P.x
    R.y = P.y - P.x
    R.z = P.z - P.y
    R.u = P.u + P.v + P.w
    R.v = P.v + P.w
    R.w = P.w
    R.branch = 200
    return Sort(R)
cdef inline PairPoint3d _Sorted_Meester(PairPoint3d P):
    r"""
    EXAMPLES::

    """
    # Apply the algo
    P.y -= P.x
    P.z -= P.x
    P.u += P.v + P.w
    P.branch = 100
    return Sort(P)


cdef inline PairPoint3d Sort(PairPoint3d P):
    r"""
    EXAMPLES::

    """
    cdef double tmp
    if P.x > P.y:
        tmp = P.x
        P.x = P.y
        P.y = tmp
        tmp = P.u
        P.u = P.v
        P.v = tmp
        P.branch += 1
    if P.y > P.z:
        tmp = P.y
        P.y = P.z
        P.z = tmp
        tmp = P.v
        P.v = P.w
        P.w = tmp
        P.branch += 2
    if P.x > P.y:
        tmp = P.x
        P.x = P.y
        P.y = tmp
        tmp = P.u
        P.u = P.v
        P.v = tmp
        P.branch += 4
    return P

cdef inline (double, double) projection3to2(double x, double y, double z):
    cdef double s = -SQRT3SUR2 * x + SQRT3SUR2 * y
    cdef double t = -.5 * x -.5 * y + z
    return s,t
