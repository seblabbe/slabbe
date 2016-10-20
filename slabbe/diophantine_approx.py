# -*- coding: utf-8 -*-
r"""
Simultaneous diophantine approximation

EXAMPLES::

    sage: from slabbe.diophantine_approx import best_simultaneous_convergents
    sage: it = best_simultaneous_convergents([e, pi])
    sage: dirichlet10 = [next(it) for _ in range(10)]
    sage: dirichlet10
    [(3, 3, 1),
     (19, 22, 7),
     (1843, 2130, 678),
     (51892, 59973, 19090),
     (113018, 130618, 41577),
     (114861, 132748, 42255),
     (166753, 192721, 61345),
     (446524, 516060, 164267),
     (1174662, 1357589, 432134),
     (3970510, 4588827, 1460669)]

The above Dirichlet simultaneous diophantine approximations appear as columns of
matrices generated by multidimensional continued fraction algorithms::

    sage: from slabbe.diophantine_approx import mult_cont_frac_vs_dirichlet
    sage: from slabbe.mult_cont_frac import Brun,ARP,Reverse,Selmer,Cassaigne
    sage: algos = [Brun(), ARP(), Reverse(), Selmer(),Cassaigne()]
    sage: mult_cont_frac_vs_dirichlet([e,pi], dirichlet10, algos)
      Dirichlet                     Brun       ARP        Reverse   Selmer     Cassaigne
    +-----------------------------+----------+----------+---------+----------+-----------+
      (3, 3, 1)                     [4, 5]     [3]        []        [8, 12]    []
      (19, 22, 7)                   [9, 16]    [6, 11]    [7, 33]   [32]       [15, 25]
      (1843, 2130, 678)             [22, 27]   [16, 17]   []        [44, 48]   [36, 39]
      (51892, 59973, 19090)         []         []         []        [56]       []
      (113018, 130618, 41577)       []         []         []        []         []
      (114861, 132748, 42255)       [33, 35]   [22, 24]   []        [62, 66]   [51, 53]
      (166753, 192721, 61345)       []         []         []        []         []
      (446524, 516060, 164267)      [36]       [25]       []        []         [56, 57]
      (1174662, 1357589, 432134)    [39, 44]   [26, 29]   []        [68]       [61, 66]
      (3970510, 4588827, 1460669)   []         [28]       []        []         []

The indices in the table are the i-th matrices. For example, the first 3
Dirichlet approximations appear in the convergents of ARP algorithm, but not the
4th neither the 5th)::

    sage: algo = ARP()
    sage: algo.n_matrix((e,pi,1), 3)
    [3 3 2]
    [3 4 2]
    [1 1 1]
    sage: algo.n_matrix((e,pi,1), 6)
    [33 19  8]
    [38 22  9]
    [12  7  3]
    sage: algo.n_matrix((e,pi,1), 16)
    [1631 2498 1843]
    [1885 2887 2130]
    [ 600  919  678]
    sage: algo.n_matrix((e,pi,1), 22)
    [114861 101941  64812]
    [132748 117816  74905]
    [ 42255  37502  23843]
    sage: algo.n_matrix((e,pi,1), 25)
    [446524 331663 842999]
    [516060 383312 974277]
    [164267 122012 310122]
    sage: algo.n_matrix((e,pi,1), 26)
    [1621186  331663 1174662]
    [1873649  383312 1357589]
    [ 596401  122012  432134]
    sage: algo.n_matrix((e,pi,1), 28)
    [3970510 2680987 1174662]
    [4588827 3098490 1357589]
    [1460669  986280  432134]

The Dirichlet vectors can be computed from the precedent ones. So one could
expect a perfect MCF algorithm based on 3x3 matrices::

    sage: from slabbe.diophantine_approx import dirichlet_convergents_dependance
    sage: dirichlet_convergents_dependance([e,pi], 8)
      i   vi                         lin. rec.     remainder
    +---+--------------------------+-------------+-----------+
      0   (3, 3, 1)                  []            (3, 3, 1)
      1   (19, 22, 7)                [6]           (1, 4, 1)
      2   (1843, 2130, 678)          [96, 6]       (1, 0, 0)
      3   (51892, 59973, 19090)      [28, 15, 1]   (0, 0, 0)
      4   (113018, 130618, 41577)    [2, 5, 1]     (0, 0, 0)
      5   (114861, 132748, 42255)    [1, 0, 1]     (0, 0, 0)
      6   (166753, 192721, 61345)    [1, 0, 1]     (0, 0, 0)
      7   (446524, 516060, 164267)   [2, 0, 1]     (0, 0, 0)

BENCHMARKS::

    sage: it = best_simultaneous_convergents([e,pi])
    sage: %time L = [next(it) for _ in range(15)]   # not tested (4.66s)
    sage: it = best_simultaneous_convergents([e,pi])
    sage: %time L = [next(it) for _ in range(18)]   # not tested (52s)

AUTHORS:

- Sébastien Labbé, September 22, 2016
- Sébastien Labbé, October 10, 2016
- Sébastien Labbé, October 19, 2016, Cython improvements (10000x faster)
- Sébastien Labbé, October 21, 2016, Parallelization of computations

TODO:

- In general, how many of the dirichlet approximations are found by the MCF algos?
"""
#*****************************************************************************
#       Copyright (C) 2016 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import itertools
from collections import defaultdict
from bisect import bisect
from sage.functions.other import floor
from sage.misc.functional import round
from sage.modules.free_module_element import vector

def frac(x):
    r"""
    Return the fractional part of real number x.

    Not always perfect...

    EXAMPLES::

        sage: from slabbe.diophantine_approx import frac
        sage: frac(3.2)
        0.200000000000000
        sage: frac(-3.2)
        0.800000000000000
        sage: frac(pi)
        pi - 3
        sage: frac(pi).n()
        0.141592653589793

    This looks suspicious...::

        sage: frac(pi*10**15).n()
        0.000000000000000
    """
    return x - floor(x)

def distance_to_nearest_integer(x):
    raise NotImplementedError

def _best_simultaneous_convergents_upto(v, Q, start=1, verbose=False):
    r"""
    Return the smallest best simultaneous diophantine approximation of vector
    ``v`` at distance ``1/Q``.

    .. WARNING::

        This is a slow function (Python) used for comparison only. Use the
        cython function of the same name in the cython module.

    INPUT:

    - ``v`` -- list of real numbers
    - ``Q`` -- real number, Q>1
    - ``start`` -- integer (default: ``1``), starting value to check
    - ``verbose`` -- boolean (default: ``False``)

    OUTPUT:

    - A tuple (u, R) where u=(p_1, ..., p_n, q) is a vector, q is an
      integer, R is a real number such that coordinates of vector p/q are
      at most 1/R from the coordinates of v.

    EXAMPLES::

        sage: from slabbe.diophantine_approx import _best_simultaneous_convergents_upto
        sage: _best_simultaneous_convergents_upto([e,pi], 2)
        ((3, 3, 1), 3.54964677830384)
        sage: _best_simultaneous_convergents_upto([e,pi], 4)
        ((19, 22, 7), 35.7490143326079)
        sage: _best_simultaneous_convergents_upto([e,pi], 35)
        ((19, 22, 7), 35.7490143326079)
        sage: _best_simultaneous_convergents_upto([e,pi], 36)
        ((1843, 2130, 678), 203.239442934072)
        sage: _best_simultaneous_convergents_upto([e,pi], 203)   # long time (1s)
        ((1843, 2130, 678), 203.239442934072)

    We can start the next computation at step 678::

        sage: _best_simultaneous_convergents_upto([e,pi], 204, start=678) # not tested (25s)
        ((51892, 59973, 19090), 266.167750949912)

    TESTS::

        sage: _best_simultaneous_convergents_upto([1,e,pi], 1)
        Traceback (most recent call last):
        ...
        ValueError: argument Q(=1) must be > 1
    """
    if not Q > 1:
        raise ValueError("argument Q(={}) must be > 1".format(Q))
    v = vector(v)
    d = len(v)
    Qinv = 1. / Q
    un_moins_Qinv = 1 - Qinv
    for q in range(start, Q**d):
        q_v = q*v
        frac_q_v = map(frac, q_v)
        if verbose:
            print q,[a.n() for a in frac_q_v]
        if all(a <= Qinv or un_moins_Qinv <= a for a in frac_q_v):
            p = map(round, q_v)
            p.append(q)
            error = max((a if a < .5 else 1-a) for a in frac_q_v)
            return tuple(p), ~error.n()
    else:
        raise RuntimeError('Did not find diophantine approximation of vector '
                'v={} with parameter Q={}'.format(v, Q))

def best_simultaneous_convergents_upto(v, Q, start=1):
    r"""
    Return a list of all best simultaneous diophantine approximations p,q of vector
    ``v`` at distance ``|qv-p|<=1/Q`` such that `1<=q<Q^d`.

    INPUT:

    - ``v`` -- list of real numbers
    - ``Q`` -- real number, Q>1
    - ``start`` -- integer (default: ``1``), starting value to check

    EXAMPLES::

        sage: from slabbe.diophantine_approx import best_simultaneous_convergents_upto
        sage: best_simultaneous_convergents_upto([e,pi], 2)
        [((3, 3, 1), 3.55)]
        sage: best_simultaneous_convergents_upto([e,pi], 4)
        [((19, 22, 7), 35.75)]
        sage: best_simultaneous_convergents_upto([e,pi], 36, start=4**2)
        [((1843, 2130, 678), 203.24)]
        sage: best_simultaneous_convergents_upto([e,pi], 204, start=36**2)
        [((51892, 59973, 19090), 266.17), ((113018, 130618, 41577), 279.19)]
        sage: best_simultaneous_convergents_upto([e,pi], 280, start=204**2)
        [((114861, 132748, 42255), 412.79), ((166753, 192721, 61345), 749.36)]
        sage: best_simultaneous_convergents_upto([e,pi], 750, start=280**2)
        [((446524, 516060, 164267), 896.47), ((1174662, 1357589, 432134), 2935.31)]
        sage: best_simultaneous_convergents_upto([e,pi], 2936, start=750**2)
        [((3970510, 4588827, 1460669), 3654.3),
         ((21640489, 25010505, 7961091), 6257.09)]
    """
    from slabbe.diophantine_approx_pyx import good_simultaneous_convergents_upto
    from sage.parallel.decorate import parallel
    from sage.parallel.ncpus import ncpus
    step = ncpus()
    @parallel
    def F(shift):
        return good_simultaneous_convergents_upto(v, Q, start+shift, step=step)
    shifts = range(step)
    goods = []
    for (arg, kwds), output in F(shifts):
        goods.extend(output)
    if not goods:
        raise ValueError("Did not find an approximation p,q to v={} s.t.  |p-qv|<=1/Q"
            " where Q={}".format(v,Q))
    goods.sort()
    sol = p, best_error_inv = goods[0]
    bests = [sol]
    # keep only the best ones
    for p,error_inv in goods[1:]:
        if error_inv > best_error_inv:
            bests.append((p,error_inv))
            best_error_inv = error_inv
    return bests

def best_simultaneous_convergents(v):  
    r"""
    Return an iterator of best convergents to a vector of real number according
    to Dirichlet theorem on simultaneous approximations.

    INPUT:

    - ``v`` -- list of real numbers

    OUTPUT:

    - iterator

    EXAMPLES::

        sage: from slabbe.diophantine_approx import best_simultaneous_convergents
        sage: it = best_simultaneous_convergents([e, pi])
        sage: [next(it) for _ in range(5)]
        [(3, 3, 1),
         (19, 22, 7),
         (1843, 2130, 678),
         (51892, 59973, 19090),
         (113018, 130618, 41577)]

    TESTS:

    Correspondance with continued fraction when d=1::

        sage: it = best_simultaneous_convergents([e])
        sage: [next(it) for _ in range(10)]
        [(3, 1),
         (8, 3),
         (11, 4),
         (19, 7),
         (87, 32),
         (106, 39),
         (193, 71),
         (1264, 465),
         (1457, 536),
         (2721, 1001)]
        sage: continued_fraction(e).convergents()[:11].list()
        [2, 3, 8/3, 11/4, 19/7, 87/32, 106/39, 193/71, 1264/465, 1457/536, 2721/1001]
    """
    Q = 2.
    start = 1
    d = len(v)
    while True:
        bests = best_simultaneous_convergents_upto(v, Q, start)
        start = Q**d
        for u,Q in bests:
            yield u
        Q += 0.001 # make sure we do not get the same again

def dirichlet_convergents_dependance(v, n, verbose=False):
    r"""
    INPUT:

    - ``v`` -- list of real numbers
    - ``n`` -- integer, number of iterations
    - ``verbose`` -- bool (default: ``False``),

    OUTPUT:

    - table of linear combinaisons of dirichlet approximations in terms of
      previous dirichlet approximations

    EXAMPLES::

        sage: from slabbe.diophantine_approx import dirichlet_convergents_dependance
        sage: dirichlet_convergents_dependance([e,pi], 4)
          i   vi                      lin. rec.     remainder
        +---+-----------------------+-------------+-----------+
          0   (3, 3, 1)               []            (3, 3, 1)
          1   (19, 22, 7)             [6]           (1, 4, 1)
          2   (1843, 2130, 678)       [96, 6]       (1, 0, 0)
          3   (51892, 59973, 19090)   [28, 15, 1]   (0, 0, 0)

    The last 3 seems enough::

        sage: dirichlet_convergents_dependance([e,pi], 8)
          i   vi                         lin. rec.     remainder
        +---+--------------------------+-------------+-----------+
          0   (3, 3, 1)                  []            (3, 3, 1)
          1   (19, 22, 7)                [6]           (1, 4, 1)
          2   (1843, 2130, 678)          [96, 6]       (1, 0, 0)
          3   (51892, 59973, 19090)      [28, 15, 1]   (0, 0, 0)
          4   (113018, 130618, 41577)    [2, 5, 1]     (0, 0, 0)
          5   (114861, 132748, 42255)    [1, 0, 1]     (0, 0, 0)
          6   (166753, 192721, 61345)    [1, 0, 1]     (0, 0, 0)
          7   (446524, 516060, 164267)   [2, 0, 1]     (0, 0, 0)

    But not in this case::

        sage: dirichlet_convergents_dependance([pi,sqrt(3)], 12)
          i    vi                      lin. rec.                  remainder
        +----+-----------------------+--------------------------+-----------+
          0    (3, 2, 1)               []                         (3, 2, 1)
          1    (22, 12, 7)             [6]                        (4, 0, 1)
          2    (47, 26, 15)            [2, 1]                     (0, 0, 0)
          3    (69, 38, 22)            [1, 1]                     (0, 0, 0)
          4    (176, 97, 56)           [2, 0, 1, 4]               (4, 1, 1)
          5    (223, 123, 71)          [1, 0, 1]                  (0, 0, 0)
          6    (399, 220, 127)         [1, 1]                     (0, 0, 0)
          7    (1442, 795, 459)        [3, 1, 0, 0, 0, 1]         (0, 0, 0)
          8    (6390, 3523, 2034)      [4, 1, 1]                  (0, 0, 0)
          9    (26603, 14667, 8468)    [4, 0, 2, 1, 0, 0, 0, 1]   (0, 0, 0)
          10   (32993, 18190, 10502)   [1, 1]                     (0, 0, 0)
          11   (40825, 22508, 12995)   [1, 0, 1, 1]               (0, 0, 0)

    The v4 is not a lin. comb. of the previous four::

        sage: dirichlet_convergents_dependance([e,pi,sqrt(3)], 5)
          i   vi                                lin. rec.        remainder
        +---+---------------------------------+----------------+--------------+
          0   (3, 3, 2, 1)                      []               (3, 3, 2, 1)
          1   (19, 22, 12, 7)                   [6]              (1, 4, 0, 1)
          2   (193, 223, 123, 71)               [10, 1]          (0, 0, 1, 0)
          3   (5529, 6390, 3523, 2034)          [28, 6, 3]       (2, 5, 1, 1)
          4   (163067, 188461, 103904, 59989)   [29, 14, 1, 1]   (2, 4, 1, 1)
    """
    L = []
    it = best_simultaneous_convergents(v)
    rows = []
    for i in range(n):
        vi = vector(next(it))
        t = vi
        M = []
        for u in reversed(L):
            m = floor(min(a/b for a,b in zip(t,u)))
            M.append(m)
            t -= m*u
            if t == 0:
                if verbose:
                    c = ','.join("v{}".format(len(L)-j) for j in range(len(M)))
                    print "v{} = {} = <{}>.<{}>".format(i, vi, M, c)
                break
        else:
            if verbose:
                print "v{} = {} = <{}>.<v{}, ..., v0> + {}".format(i, vi, M, i-1, t)
        L.append(vi)
        row = [i, vi, M, t]
        rows.append(row)
    header_row = ['i', 'vi', 'lin. rec.', 'remainder']
    from sage.misc.table import table
    return table(rows=rows, header_row=header_row)

def mult_cont_frac_vs_dirichlet_dict(v, dirichlet, algos):
    r"""
    INPUT:

    - ``v`` -- list of real numbers
    - ``dirichlet`` -- list, first dirichlet approximations
    - ``algos`` -- list, list of mult. cont. frac. algorithms

    OUTPUT:

        dict

    EXAMPLES::

        sage: from slabbe.diophantine_approx import best_simultaneous_convergents
        sage: from slabbe.diophantine_approx import mult_cont_frac_vs_dirichlet_dict
        sage: from slabbe.mult_cont_frac import ARP, Brun
        sage: v = [e, pi]
        sage: it = best_simultaneous_convergents(v)
        sage: dirichlet = [next(it) for _ in range(3)]
        sage: mult_cont_frac_vs_dirichlet_dict([e,pi], dirichlet, [Brun(), ARP()])
        {Arnoux-Rauzy-Poincar\'e 3-dimensional continued fraction algorithm:
         defaultdict(<type 'list'>, {(19, 22, 7): [6, 7, 8, 9, 10, 11],
         (3, 3, 1): [3], (1843, 2130, 678): [16, 17]}),
         Brun 3-dimensional continued fraction algorithm:
         defaultdict(<type 'list'>, {(19, 22, 7): [9, 10, 11, 12, 13,
         14, 15, 16], (3, 3, 1): [4, 5], (1843, 2130, 678): [22, 23,
         24, 25, 26, 27]})}

     Or from precomputed Dirichlet approximations::

        sage: dirichlet = [(3, 3, 1), (19, 22, 7), (1843, 2130, 678), (51892, 59973, 19090)]
        sage: mult_cont_frac_vs_dirichlet_dict([e,pi], dirichlet, [Brun(), ARP()])
        {Arnoux-Rauzy-Poincar\'e 3-dimensional continued fraction algorithm:
         defaultdict(<type 'list'>, {(19, 22, 7): [6, 7, 8, 9, 10, 11],
         (3, 3, 1): [3], (1843, 2130, 678): [16, 17]}),
         Brun 3-dimensional continued fraction algorithm:
         defaultdict(<type 'list'>, {(19, 22, 7): [9, 10, 11, 12, 13,
         14, 15, 16], (3, 3, 1): [4, 5], (1843, 2130, 678): [22, 23,
         24, 25, 26, 27]})}
    """
    dirichlet_set = set(tuple(u) for u in dirichlet)
    MAX = max(u[-1] for u in dirichlet)
    vv = list(v) + [1]
    D = {}
    for algo in algos:
        D[algo] = defaultdict(list)
        for i in itertools.count():
            columns = algo.n_matrix(vv, i).columns()
            if min(col[-1] for col in columns) > MAX:
                break
            for col in columns:
                col = tuple(col)
                if col in dirichlet_set:
                    D[algo][col].append(i)
    return D

def mult_cont_frac_vs_dirichlet(v, dirichlet, algos):
    r"""
    Returns the indices i such that dirichlet approximations appears as columns
    of the i-th matrix obtained from mult. dim. cont. frac. algorithms.

    INPUT:

    - ``v`` -- list of real numbers
    - ``dirichlet`` -- list, first dirichlet approximations
    - ``algos`` -- list, list of mult. cont. frac. algorithms

    OUTPUT:

    table

    EXAMPLES::

        sage: from slabbe.diophantine_approx import best_simultaneous_convergents
        sage: from slabbe.diophantine_approx import mult_cont_frac_vs_dirichlet
        sage: from slabbe.mult_cont_frac import Brun,ARP,Reverse,Selmer,Cassaigne
        sage: v = [e, pi]
        sage: it = best_simultaneous_convergents(v)
        sage: dirichlet = [next(it) for _ in range(10)]
        sage: algos = [Brun(), ARP(), Reverse(), Selmer(),Cassaigne()]
        sage: mult_cont_frac_vs_dirichlet(v, dirichlet, algos)
          Dirichlet                     Brun       ARP        Reverse   Selmer     Cassaigne
        +-----------------------------+----------+----------+---------+----------+-----------+
          (3, 3, 1)                     [4, 5]     [3]        []        [8, 12]    []
          (19, 22, 7)                   [9, 16]    [6, 11]    [7, 33]   [32]       [15, 25]
          (1843, 2130, 678)             [22, 27]   [16, 17]   []        [44, 48]   [36, 39]
          (51892, 59973, 19090)         []         []         []        [56]       []
          (113018, 130618, 41577)       []         []         []        []         []
          (114861, 132748, 42255)       [33, 35]   [22, 24]   []        [62, 66]   [51, 53]
          (166753, 192721, 61345)       []         []         []        []         []
          (446524, 516060, 164267)      [36]       [25]       []        []         [56, 57]
          (1174662, 1357589, 432134)    [39, 44]   [26, 29]   []        [68]       [61, 66]
          (3970510, 4588827, 1460669)   []         [28]       []        []         []
    """
    D = mult_cont_frac_vs_dirichlet_dict(v, dirichlet, algos)
    rows = []
    for v in dirichlet:
        v = tuple(v)
        row = [v] 
        for algo in algos:
            indices = D[algo][v]
            if len(indices) > 1:
                row.append([min(indices), max(indices)])
            else:
                row.append(indices)
        rows.append(row)
    header_row = ['Dirichlet'] + [algo.class_name() for algo in algos]
    from sage.misc.table import table
    return table(rows=rows, header_row=header_row)




