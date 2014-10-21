# -*- coding: utf-8 -*-
r"""
Multidimensional Continued Fraction Algorithms

Gain d'optimisation in Cython:

    - Declare x,y,z as double instead of list : factor of 4
    - Faire les calculs dans une boucle : factor of 2
    - Do not use yield : factor of 10

TODO:

    - reuteneaour, nogueira, autres?
    - make a other class for 2d and and 1d methods
    - Read [1] and change cpdef int C[NDIVS][NDIVS][NDIVS]
    - faire une fonction sort?
    - utilise les vecteurs pour les plus grandes dimensions?

    - make use of enum of algorithms instead of char* branch

        cdef enum CheeseState:
            hard = 1
            soft = 2
            runny = 3

        cpdef int ff ():
            cdef CheeseState x = hard
            if x == 1:
                return 45
            else:
                return 88

    - try puttin arguments of inline methods as references
    - or avoid creating a new pairpoint R, the copy is already done by
      default

Question:

    - Comment factoriser le code sans utiliser les yield?
    - Comment faire un appel de fonction rapide (pour factoriser le code)

[1] https://groups.google.com/forum/?fromgroups=#!topic/sage-devel/NCBmj2KjwEM

AUTHORS:

 - Vincent Delecroix, C code, Computation of Lyapounov exponents for Brun
   algorithm, June 2013.
 - Sebastien Labbe, Invariant measures, Lyapounov exponents and natural
   extensions for a dozen of algorithms, October 2013.

"""
#*****************************************************************************
#       Copyright (C) 2014 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import itertools
import collections
from math import log
from __builtin__ import sum
from sage.rings.real_mpfr import RealField
from sage.plot.point import point
from sage.misc.prandom import random
from sage.misc.table import table

cdef struct PairPoint3d:
    double x
    double y
    double z
    double u
    double v
    double w
    int branch

cdef enum Algorithm:
    sorted_brun,
    sorted_arrevert,
    brun,
    brunmulti,
    selmer,
    meester,
    arnouxrauzy,
    arnouxrauzymulti,
    poincare,
    arp,
    arpmulti,
    arrevert,
    arrevertmulti,
    armonteil,
    delaunay,
    jacobi,
    jacobiadditif,
    jacobiadditifv2,
    sort

algo_name_to_algo_id = dict(
    sorted_brun=    sorted_brun,
    sorted_arrevert=sorted_arrevert,
    brun=           brun,
    brunmulti=      brunmulti,
    selmer=         selmer,
    meester=        meester,
    arnouxrauzy=    arnouxrauzy,
    arnouxrauzymulti = arnouxrauzymulti,
    poincare=       poincare,
    arp=            arp,
    arpmulti=       arpmulti,
    arrevert=       arrevert,
    arrevertmulti=  arrevertmulti,
    armonteil=      armonteil,
    delaunay=       delaunay,
    jacobi=         jacobi,
    jacobiadditif=  jacobiadditif,
    jacobiadditifv2=jacobiadditifv2,
    sort=           sort,
    )

cdef double SQRT3SUR2 = 0.866025403784439


class MCFAlgorithm_pyx(object):
    def __init__(self, name=None):
        self._name = name
        self._algo_id = -1
    def name(self):
        return self._name
    def algo_id(self):
        if self._algo_id == -1:
            self._algo_id = algo_name_to_algo_id[self._name]
        return self._algo_id
    def check_definition(self, int n_iterations):
        r"""
        INPUT:

        - ``n_iterations`` -- integer

        OUTPUT:

        bool

        EXAMPLES::

            sage: t = algo.arp.check_definition(10000)
            True

        """
        cdef double s,t             # temporary variables
        cdef unsigned int i         # loop counter

        cdef PairPoint3d P, R

        cdef Algorithm algo_id = self.algo_id()

        # Loop
        for i from 0 <= i < n_iterations:

            # random initial values
            P.x = random(); P.y = random(); P.z = random();
            P.u = random(); P.v = random(); P.w = random();
            s = P.x*P.u + P.y*P.v + P.z*P.w

            # Apply Algo
            R = Algo_switch(P, algo_id)
            t = R.x*R.u + R.y*R.v + R.z*R.w

            if not abs(s - t) < 0.0000001:
                m = 'Algo {} does not preserve the scalar product\n'.format(self.name())
                m += '{} != {}\n'.format(s,t)
                m += 'The problem is on branch {}\n'.format(R.branch)
                m += 'INPUT: ({}, {}, {}, {}, {}, {})\n'.format(P.x,P.y,P.z,P.u,P.v,P.w)
                m += 'OUTPUT: ({}, {}, {}, {}, {}, {})\n'.format(R.x,R.y,R.z,R.u,R.v,R.w)
                raise Exception(m)

        return True

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

        EXAMPLES

        ::

            sage: D = brun.invariant_measure_dict(4, 10, verbose=True)
            0.102209671626 0.241453007005 0.656337321369
            1 2 6
            0.134744020568 0.318309886184 0.546946093248
            1 3 5
            0.197661690901 0.335396102174 0.466942206924
            1 3 4
            0.197931587793 0.297412777066 0.504655635141
            1 2 5

        ::

            sage: D = brun.invariant_measure_dict(100000, 5)
            sage: D
            {(0, 1, 2): 21299, (0, 2, 2): 10079, (0, 0, 4): 19223, (1, 1, 1): 4912, (0, 1, 3): 19415, (0, 0, 3): 12307, (1, 1, 2): 12765}

        It is 1000 times faster using C counter instead of a python dict counter::

            sage: time D = brun.invariant_measure_dict(1000000, 10)
            Time: CPU 0.05 s, Wall: 0.05 s

        """
        # 146 works, 147 causes segmentation error!!!
        DEF NDIVS = 100
        assert ndivs <= NDIVS, "ndivs(=%s) must be less or equal to %s" % (ndivs, NDIVS)
        cdef double s
        cdef int i,j
        cdef int X,Y
        cdef Algorithm algo_id = self.algo_id()

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
            # Apply Algo
            P = Algo_switch(P, algo_id)

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

    def natural_extension(self, int n_iterations, norm_algo='sup', norm_ext='1', verbose=False):
        r"""
        INPUT:

        - ``n_iterations`` -- integer
        - ``norm_algo`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit of the algo
        - ``norm_ext`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'`` or ``'hypersurfac'``, the norm used for the orbit in the natural extension
        - ``verbose`` -- bool (default: ``False``)

        OUTPUT:

        dict, dict, dict, dict

        EXAMPLES::

            sage: t = algo.arp.natural_extension(10000, a=(1,2,3))
            sage: map(type, t)
            [collections.defaultdict,
             collections.defaultdict,
             collections.defaultdict,
             collections.defaultdict]

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
        cdef Algorithm algo_id = self.algo_id()
        P.x = x
        P.y = y
        P.z = z
        P.u = u
        P.v = v
        P.w = w

        dual_in = collections.defaultdict(list)
        dual_out = collections.defaultdict(list)
        domain_in = collections.defaultdict(list)
        domain_out = collections.defaultdict(list)

        # Loop
        for i from 0 <= i < n_iterations:

            # Apply Algo
            R = Algo_switch(P, algo_id)

            if verbose:
                print("x=%f, y=%f, z=%f" % (R.x,R.y,R.z))
                #print("u=%f, v=%f, w=%f" % (u,v,w))
                #s = x*u + y*v + z*w
                #print("scal prod <(x,y,z),(u,v,w)> = %f (after algo)" % s)

            # Normalize (xnew,ynew,znew)
            if norm_algo == '1':
                s = R.x + R.y + R.z # norm 1
            elif norm_algo == 'sup':
                s = max(R.x, R.y, R.z) # norm sup
            else:
                raise ValueError("Unknown value for norm_algo(=%s)" %norm_algo)
            R.x /= s; R.y /= s; R.z /= s

            # Normalize (unew,vnew,wnew)
            if norm_ext == '1':
                s = R.u + R.v + R.w    # norm 1
            elif norm_ext == 'sup':
                s = max(R.u, R.v, R.w) # norm sup
            elif norm_ext == 'hypersurface':
                s = R.x*R.u + R.y*R.v + R.z*R.w # hypersurface
            else:
                raise ValueError("Unknown value for norm_ext(=%s)" %norm_ext)
            R.u /= s; R.v /= s; R.w /= s

            # Projection
            if norm_algo == '1':
                s = -SQRT3SUR2 * P.x + SQRT3SUR2 * P.y
                t = -.5 * P.x -.5 * P.y + P.z
                domain_in[R.branch].append((s,t))
                s = -SQRT3SUR2 * R.x + SQRT3SUR2 * R.y
                t = -.5 * R.x -.5 * R.y + R.z
                domain_out[R.branch].append((s,t))
            elif norm_algo == 'sup':
                domain_in[R.branch].append((P.x,P.y))
                domain_out[R.branch].append((R.x,R.y))

            if norm_ext == '1':
                s = -SQRT3SUR2 * P.u + SQRT3SUR2 * P.v
                t = -.5 * P.u -.5 * P.v + P.w
                dual_in[R.branch].append((s,t))
                s = -SQRT3SUR2 * R.u + SQRT3SUR2 * R.v
                t = -.5 * R.u -.5 * R.v + R.w
                dual_out[R.branch].append((s,t))
            elif norm_ext == 'sup':
                dual_in[R.branch].append((P.u,P.v))
                dual_out[R.branch].append((R.u,R.v))

            P = R

        return domain_in, domain_out, dual_in, dual_out

    def lyapounov_exponents(self, int n_iterations=1000, verbose=False):
        r"""
        INPUT:

        - ``n_iterations`` -- integer
        - ``verbose`` -- bool (default: ``False``)

        OUTPUT:

        liapounov exponents

        theta1, theta2, theta2/theta1

        NOTE:: the code of this method was translated from C to cython. The
        C version is from Vincent Delecroix.

        EXAMPLES::

        Some benchmarks (at liafa)::

            sage: time brun.lyapounov_exponents(1000000, step=16)
            (0.3049429393152174, -0.1120652699014143, -0.367495867105725)
            Time: CPU 0.05 s, Wall: 0.06 s

        Some benchmarks (on my machine)::

            sage: %timeit brun.lyapounov_exponents(1000000, 16)
            10 loops, best of 3: 68.6 ms per loop

        Cython code on liafa is as fast as C on my machine::

            sage: time brun.lyapounov_exponents(67000000, step=16)
            (0.30452120021265766, -0.11212586210856369, -0.36820379674801734)
            Time: CPU 3.71 s, Wall: 3.71 s

        Cython code on my machine is almost as fast as C on my machine::

            sage: %time brun.lyapounov_exponents(67000000, 16)
            CPU times: user 4.56 s, sys: 0.01 s, total: 4.57 s
            Wall time: 4.58 s
            (0.30456433843239084, -0.1121770192467067, -0.36831961293987303)

        It is even faster using struct and an outside function::

            sage: %time algo.brun.lyapounov_exponents_lazy(67000000)
            CPU times: user 3.83 s, sys: 0.01 s, total: 3.84 s
            Wall time: 3.86 s
            (0.3044828923859014, -0.11213201206117002, -0.36827031949977)

        """
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
        cdef Algorithm algo_id = self.algo_id()
        P.x = x
        P.y = y
        P.z = z
        P.u = u
        P.v = v
        P.w = w

        # Loop
        for i from 0 <= i < n_iterations:

            # Apply Algo
            P = Algo_switch(P, algo_id)

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

        return theta1/n_iterations, theta2/n_iterations, theta2/theta1


    def dual_domain(self):
        r"""
        Return the dual domain of each branch.

        EXAMPLES::

            sage: algo.brun.dual_domain()
            {100: [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.5, 0.5)],
             102: [(1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.5, 0.5)],
             106: [(0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.5, 0.0, 0.5)]}
            sage: algo.arp.dual_domain()
            {100: [(1.0, 0.0, 0.0),
              (0.0, 1.0, 0.0),
              (0.3333333333333333, 0.3333333333333333, 0.3333333333333333)],
             102: [(1.0, 0.0, 0.0),
              (0.0, 0.0, 1.0),
              (0.3333333333333333, 0.3333333333333333, 0.3333333333333333)],
             106: [(0.0, 1.0, 0.0),
              (0.0, 0.0, 1.0),
              (0.3333333333333333, 0.3333333333333333, 0.3333333333333333)],
             203: [(0.0, 0.0, 1.0),
              (0.5, 0.0, 0.5),
              (0.3333333333333333, 0.3333333333333333, 0.3333333333333333)],
             206: [(0.0, 1.0, 0.0),
              (0.0, 0.5, 0.5),
              (0.3333333333333333, 0.3333333333333333, 0.3333333333333333)],
             207: [(0.0, 0.0, 1.0),
              (0.0, 0.5, 0.5),
              (0.3333333333333333, 0.3333333333333333, 0.3333333333333333)]}
        """
        cdef double x,y,z           # vector (x,y,z)
        cdef double s
        cdef PairPoint3d P, R0, R1, R2
        cdef Algorithm algo_id = self.algo_id()
        cdef int n_iterations = 100

        D = {}
        for i from 0 <= i < n_iterations:
            # random initial values
            x = random(); y = random(); z = random();

            # Order (x,y,z)
            if y > z: z,y = y,z
            if x > z: x,y,z = y,z,x
            elif x > y: x,y = y,x

            P.x,P.y,P.z = x,y,z

            # Apply Algo
            P.u, P.v, P.w = 1,0,0
            R0 = Sort(Algo_switch(P, algo_id))
            s = abs(R0.u) + abs(R0.v) + abs(R0.w);
            R0.u /= s; R0.v /= s; R0.w /= s

            # Apply Algo
            P.u, P.v, P.w = 0,1,0
            R1 = Sort(Algo_switch(P, algo_id))
            s = abs(R1.u) + abs(R1.v) + abs(R1.w);
            R1.u /= s; R1.v /= s; R1.w /= s

            # Apply Algo
            P.u, P.v, P.w = 0,0,1
            R2 = Sort(Algo_switch(P, algo_id))
            s = abs(R2.u) + abs(R2.v) + abs(R2.w);
            R2.u /= s; R2.v /= s; R2.w /= s

            assert x==P.x
            assert y==P.y
            assert z==P.z
            assert R0.branch == R1.branch == R2.branch, "problem distinct branch"

            if R0.branch not in D:
                D[R0.branch] = []
                D[R0.branch].append((R0.u,R0.v,R0.w))
                D[R0.branch].append((R1.u,R1.v,R1.w))
                D[R0.branch].append((R2.u,R2.v,R2.w))

        return D

cpdef inline PairPoint3d Algo_switch(PairPoint3d P, Algorithm algo):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, algo_name_to_algo_id['brun'])
        {'branch': 'brun1', 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.3, 'z': 0.5}

    ::

        sage: D = {'x':.7,'y':.6,'z':.5,'u':.1,'v':.2,'w':.3,'branch':999}
        sage: Algo_switch(D, algo_name_to_algo_id['sort'])
        {'branch': 'ad', 'u': 0.3, 'v': 0.2, 'w': 0.1, 'x': 0.5, 'y': 0.6, 'z': 0.7}
    """
    if algo == sorted_brun:
        return Sorted_Brun(P)
    elif algo == sorted_arrevert:
        return Sorted_ARrevert(P)
    elif algo == brun:
        return Brun(P)
    elif algo == brunmulti:
        return Sorted_BrunMulti(P)
    elif algo == selmer:
        return Sorted_Selmer(P)
    elif algo == meester:
        return Sorted_Meester(P)
    elif algo == poincare:
        return Sorted_Poincare(P)
    elif algo == arnouxrauzy:
        return Sorted_ArnouxRauzy(P)
    elif algo == arnouxrauzymulti:
        return Sorted_ArnouxRauzyMulti(P)
    elif algo == arp:
        return Sorted_ARP(P)
    elif algo == arpmulti:
        return Sorted_ARPMulti(P)
    elif algo == arrevert:
        return ARrevert(P)
    elif algo == arrevertmulti:
        return Sorted_ARrevertMulti(P)
    elif algo == armonteil:
        return Sorted_ARMonteil(P)
    elif algo == delaunay:
        return Sorted_Delaunay(P)
    elif algo == jacobi:
        return JacobiPerron(P)
    elif algo == jacobiadditif:
        return JacobiPerronAdditif(P)
    elif algo == jacobiadditifv2:
        return JacobiPerronAdditifv2(P)
    elif algo == sort:
        return Sort(P)
    else:
        raise ValueError("Unknown algo(=%s)" % algo)

cdef inline PairPoint3d Brun(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, algo_name_to_algo_id['brun'])
        {'branch': 100, 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.6, 'z': 0.20000000000000007}
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
    return P

cdef inline PairPoint3d ARrevert(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, algo_name_to_algo_id['arrevert'])
        {'branch': 100, 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.6, 'z': 0.20000000000000007}
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

cdef inline PairPoint3d Sorted_Brun(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, algo_name_to_algo_id['brun'])
        {'branch': 106, 'u': 0.3, 'v': 0.2, 'w': 0.6, 'x': 0.20000000000000007, 'y': 0.3, 'z': 0.6}
        sage: D = {'x':.3,'y':.45,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, algo_name_to_algo_id['brun'])
        {'branch': 102, 'u': 0.2, 'v': 0.3, 'w': 0.6, 'x': 0.3, 'y': 0.35000000000000003, 'z': 0.45}
        sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, algo_name_to_algo_id['brun'])
        {'branch': 100, 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.3, 'z': 0.5}
    """
    P.z -= P.y
    P.v += P.w
    P.branch = 100
    return Sort(P)

cdef inline PairPoint3d Sorted_BrunMulti(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, 'brunmulti')
        {'branch': 'brun1', 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.3, 'z': 0.5}
    """
    cdef int m = <int>(P.z / P.y)
    P.z -= m*P.y
    P.v += m*P.w
    P.branch = 100
    return Sort(P)

cdef inline PairPoint3d Sorted_Selmer(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.2,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, algo_name_to_algo_id['selmer'])
        {'branch': 100, 'u': 0.5, 'v': 0.3, 'w': 0.3, 'x': 0.2, 'y': 0.3, 'z': 0.6000000000000001}
    """
    P.z -= P.x
    P.u += P.w
    P.branch = 100
    return Sort(P)

cdef inline PairPoint3d Sorted_Meester(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.5,'y':.6,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, algo_name_to_algo_id['meester'])
        {'branch': 103, 'u': 0.3, 'v': 0.3, 'w': 0.8, 'x': 0.09999999999999998, 'y': 0.30000000000000004, 'z': 0.5}
    """
    # Apply the algo
    P.y -= P.x
    P.z -= P.x
    P.u += P.v + P.w
    P.branch = 100
    return Sort(P)

cdef inline PairPoint3d Sorted_ArnouxRauzy(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = dict(x=.3,y=.4,z=.8,u=.2,v=.3,w=.4,branch=999)
        sage: Algo_switch(D, algo_name_to_algo_id['arnouxrauzy'])
        {'branch': 106, 'u': 0.4, 'v': 0.6000000000000001, 'w': 0.7, 'x': 0.10000000000000009, 'y': 0.3, 'z': 0.4}

        ::

        sage: D = dict(x=.3,y=.7,z=.8,u=.2,v=.3,w=.4,branch=999)
        sage: Algo_switch(D, algo_name_to_algo_id['arnouxrauzy'])
        {'branch': 206, 'u': 0.4, 'v': 0.8999999999999999, 'w': 0.7, 'x': 0.10000000000000009, 'y': 0.3, 'z': 0.39999999999999997}

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

cdef inline PairPoint3d Sorted_ArnouxRauzyMulti(PairPoint3d P):
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

cdef inline PairPoint3d Sorted_ARP(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = dict(x=.3,y=.4,z=.8,u=.2,v=.3,w=.4,branch=999)
        sage: Algo_switch(D, algo_name_to_algo_id['arp'])
        {'branch': 106, 'u': 0.4, 'v': 0.6000000000000001, 'w': 0.7, 'x': 0.10000000000000009, 'y': 0.3, 'z': 0.4}

        ::

        sage: D = dict(x=.3,y=.7,z=.8,u=.2,v=.3,w=.4,branch=999)
        sage: Algo_switch(D, algo_name_to_algo_id['arp'])
        {'branch': 206, 'u': 0.4, 'v': 0.8999999999999999, 'w': 0.7, 'x': 0.10000000000000009, 'y': 0.3, 'z': 0.39999999999999997}
    """
    # Apply the algo
    if P.z > P.x + P.y:
        return Sorted_ArnouxRauzy(P)
    else:
        return Sorted_Poincare(P)

cdef inline PairPoint3d Sorted_ARPMulti(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.5,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, algo_name_to_algo_id['arpmulti'])
        {'branch': 201, 'u': 0.6, 'v': 0.8, 'w': 0.3, 'x': 0.2, 'y': 0.3, 'z': 0.30000000000000004}
    """
    # Apply the algo
    if P.z > P.x + P.y:
        return Sorted_ArnouxRauzyMulti(P)
    else:
        return Sorted_Poincare(P)

cdef inline PairPoint3d Sorted_Poincare(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = dict(x=.3,y=.7,z=.8,u=.2,v=.3,w=.4,branch=999)
        sage: Algo_switch(D, algo_name_to_algo_id['poincare'])
        {'branch': 206, 'u': 0.4, 'v': 0.9, 'w': 0.7, 'x': 0.10000000000000009, 'y': 0.3, 'z': 0.39999999999999997}

    ::

        sage: Algo_switch(D, algo_name_to_algo_id['poincare'])
        {'branch': 206, 'u': 0.10000000000000003, 'v': 1.8, 'w': 0.09999999999999998, 'x': 0.10000000000000009, 'y': 0.3, 'z': 0.39999999999999997}

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

cdef inline PairPoint3d Sorted_ARrevert(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, 'arrevert')
        {'branch': 'brun1', 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.3, 'z': 0.5}
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

cdef inline PairPoint3d Sorted_ARrevertMulti(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, 'arrevertmulti')
        {'branch': 'brun1', 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.3, 'z': 0.5}
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

cdef inline PairPoint3d Sorted_ARMonteil(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, 'armonteil')
        {'branch': 'brun1', 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.3, 'z': 0.5}
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
cdef inline PairPoint3d Sorted_Delaunay(PairPoint3d P):
    r"""
    Donné par Xavier Provençal (inspiré de en fait) le 3 février 2014.

    EXAMPLES::

        sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, 'delaunay')
        {'branch': 'brun1', 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.3, 'z': 0.5}
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
        return Sorted_Meester(P)
cdef inline PairPoint3d JacobiPerron(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, 'jacobi')
        {'branch': 'brun1', 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.3, 'z': 0.5}
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
cdef inline PairPoint3d JacobiPerronAdditif(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, 'jacobiadditif')
        {'branch': 'brun1', 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.3, 'z': 0.5}
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
cdef inline PairPoint3d JacobiPerronAdditifv2(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, 'jacobiadditifv2')
        {'branch': 'brun1', 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.3, 'z': 0.5}
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

cdef inline PairPoint3d Sort(PairPoint3d P):
    r"""
    EXAMPLES::

        sage: D = {'x':.3,'y':.3,'z':.8,'u':.2,'v':.3,'w':.3,'branch':999}
        sage: Algo_switch(D, 'brun')
        {'branch': 'brun1', 'u': 0.2, 'v': 0.6, 'w': 0.3, 'x': 0.3, 'y': 0.3, 'z': 0.5}
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

print 35

