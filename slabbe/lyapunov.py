# -*- coding: utf-8 -*-
r"""
Lyapunov parallel computation for MCF algorithms

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
from __future__ import absolute_import, print_function

from sage.misc.prandom import random
from sage.parallel.decorate import parallel
from sage.misc.table import table

def lyapunov_sample(algo, n_orbits, n_iterations=1000, ncpus=2,
        verbose=False):
    r"""
    Return lists of values for theta1, theta2 and 1-theta2/theta1 computed
    on many orbits.

    This is computed in parallel.

    INPUT:

    - ``algo`` -- MCF algorithm
    - ``n_orbits`` -- integer, number of orbits
    - ``n_iterations`` -- integer, length of each orbit
    - ``ncpus`` -- integer (default:``2``), number of cpus to use
    - ``verbose`` -- bool (default: ``False``)

    OUTPUT:

        tuple of three lists

    EXAMPLES::

        sage: from slabbe.lyapunov import lyapunov_sample
        sage: from slabbe.mult_cont_frac import Brun
        sage: lyapunov_sample(Brun(), 5, 1000000) # abs tol 0.01
        [(0.3027620661266397,
          0.3033468535021702,
          0.3044950176856005,
          0.3030531162480779,
          0.30601169862996064),
         (-0.11116236859835525,
          -0.11165563059874498,
          -0.1122595926203868,
          -0.11190323336181864,
          -0.11255687513610782),
         (1.367160820443926,
          1.3680790794750939,
          1.3686746452327765,
          1.3692528714016428,
          1.3678188632657973)]

    Works for higher dimensional algorithms::

        sage: lyapunov_sample(Brun(8), 5, 10^5) # abs tol 0.01
        [(0.24494574466175367,
          0.24492293068699247,
          0.24494468166245503,
          0.2447894172680422,
          0.2452400265773239),
         (-0.012028663893597698,
          -0.012174084683987084,
          -0.012268030566904951,
          -0.012286468841900471,
          -0.012166894088285442),
         (1.049107462185996,
          1.049705777445336,
          1.0500849027773989,
          1.0501919935061854,
          1.0496121871217023)]

    """
    dim = algo.dimension()
    S = [tuple(random() for _ in range(dim)) for _ in range(n_orbits)]
    @parallel(ncpus=ncpus)
    def compute_exponents(i):
        try:
            # When start=None, a random starting point is chosen
            # but it is always the same(!) because the random seed is
            # reinitialize to the same value each time in @parallel. 
            # This is why we need the list S of random starting points
            return algo.lyapunov_exponents(start=S[i], n_iterations=n_iterations) 
        except Exception as err:
            return "{}: {}".format(err.__class__.__name__, err)
    inputs = list(range(n_orbits))
    L = [v for _,v in compute_exponents(inputs)]
    L_filtered = [v for v in L if isinstance(v, tuple)]
    if verbose:
        L_error_msg = [v for v in L if not isinstance(v, tuple)]
        print(L_error_msg)
    return list(zip(*L_filtered))

def lyapunov_table(algo, n_orbits, n_iterations=1000, ncpus=2):
    r"""
    Return a table of values of Lyapunov exponents for this algorithm.

    INPUT:

    - ``algo`` -- MCF algorithm
    - ``n_orbits`` -- integer, number of orbits
    - ``n_iterations`` -- integer, length of each orbit
    - ``ncpus`` -- integer (default:``2``), number of cpus to use

    OUTPUT:

        table of liapounov exponents

    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Brun
        sage: from slabbe.lyapunov import lyapunov_table
        sage: lyapunov_table(Brun(), 10, 1000000) # random
          10 succesful orbits     min       mean      max       std
        +-----------------------+---------+---------+---------+---------+
          $\theta_1$              0.303     0.305     0.307     0.0013
          $\theta_2$              -0.1131   -0.1124   -0.1115   0.00051
          $1-\theta_2/\theta_1$   1.3678    1.3687    1.3691    0.00043

    Works for higher dimensional algorithms::

        sage: lyapunov_table(Brun(8), 10, 10^6) # random
          10 succesful orbits     min        mean       max        std
        +-----------------------+----------+----------+----------+----------+
          $\theta_1$              0.24491    0.24500    0.24506    0.000041
          $\theta_2$              -0.01230   -0.01211   -0.01198   0.000096
          $1-\theta_2/\theta_1$   1.0489     1.0494     1.0502     0.00040

    """
    import numpy as np
    from sage.misc.functional import numerical_approx
    from sage.rings.real_mpfr import RR
    from sage.misc.table import table
    rep = lyapunov_sample(algo, n_orbits, n_iterations, ncpus=ncpus)
    def floor_log(number):
        return RR(number).abs().log10().floor()
    def my_rounded(number, s):
        m = floor_log(number)
        return numerical_approx(number, digits=m-s+1)
    names = [r"$\theta_1$", r"$\theta_2$", r"$1-\theta_2/\theta_1$"]
    rows = []
    for i, data in enumerate(rep):
        data = np.array(data)
        s = floor_log(data.std())
        row = [names[i]]
        row.append(my_rounded(data.min(),s))
        row.append(my_rounded(data.mean(),s))
        row.append(my_rounded(data.max(),s))
        row.append(my_rounded(data.std(),s))
        rows.append(row)
    header = ['{} succesful orbits'.format(len(rep[0])), 'min','mean','max','std']
    return table(rows=rows,header_row=header)

def _lyapunov_row(algo, n_orbits, n_iterations=1000, ncpus=2):
    r"""
    Return a row of values of Lyapunov exponents.

    This method is used for making a larger table comparing many
    algorithms.

    INPUT:

    - ``algo`` -- MCF algorithm
    - ``n_orbits`` -- integer, number of orbits
    - ``n_iterations`` -- integer, length of each orbit
    - ``ncpus`` -- integer (default:``2``), number of cpus to use

    OUTPUT:

        row of liapounov exponents: theta1, theta2, 1-theta2/theta1

    EXAMPLES::

        sage: from slabbe.lyapunov import _lyapunov_row
        sage: from slabbe.mult_cont_frac import Brun
        sage: _lyapunov_row(Brun(), 10, 10^6) # abs tol 0.01
        ['Brun (d=3)', 10, '0.305 (0.0013)', '-0.1121 (0.00061)', '1.368 (0.0010)']

    Works for higher dimensional algorithms::

        sage: _lyapunov_row(Brun(8), 10, 10^6) # abs tol 0.01
        ['Brun (d=8)',
         10,
         '0.24502 (0.000058)',
         '-0.01204 (0.000069)',
         '1.0491 (0.00027)']

    """
    import numpy as np
    from sage.misc.functional import numerical_approx
    from sage.rings.real_mpfr import RR
    rep = lyapunov_sample(algo, n_orbits, n_iterations, ncpus=ncpus)
    def floor_log(number):
        return RR(number).abs().log10().floor()
    def my_rounded(number, s):
        m = floor_log(number)
        return numerical_approx(number, digits=m-s+1)
    row = []
    row.append("{} (d={})".format(algo.name(),algo.dimension()))
    row.append(len(rep[0]))
    for data in rep:
        data = np.array(data)
        s = floor_log(data.std())
        val = my_rounded(data.mean(), s)
        std = my_rounded(data.std(), s)
        row.append("{} ({})".format(val, std))
    return row

def lyapunov_comparison_table(L, n_orbits=100, n_iterations=10000, ncpus=2):
    r"""
    Return a table of values of Lyapunov exponents for many algorithm.

    INPUT:

    - ``L`` -- list of algorithms
    - ``n_orbits`` -- integer
    - ``n_iterations`` -- integer
    - ``ncpus`` -- integer (default:``2``), number of cpus to use

    OUTPUT:

        table

    EXAMPLES::

        sage: import slabbe.mult_cont_frac as mcf
        sage: from slabbe.lyapunov import lyapunov_comparison_table
        sage: algos = [mcf.Brun(), mcf.ARP()]
        sage: lyapunov_comparison_table(algos)    # abs tol 0.01
          Algorithm                       \#Orbits   $\theta_1$ (std)   $\theta_2$ (std)   $1-\theta_2/\theta_1$ (std)
        +-------------------------------+----------+------------------+------------------+-----------------------------+
          Arnoux-Rauzy-Poincar\'e (d=3)   100        0.44 (0.014)       -0.173 (0.0060)    1.389 (0.0051)
          Brun (d=3)                      100        0.305 (0.0085)     -0.112 (0.0042)    1.368 (0.0073)


    Works for higher dimensional algorithms::

        sage: algos = [mcf.Brun(a) for a in range(3,6)]
        sage: lyapunov_comparison_table(algos)    # abs tol 0.01
          Algorithm    \#Orbits   $\theta_1$ (std)   $\theta_2$ (std)   $1-\theta_2/\theta_1$ (std)
        +------------+----------+------------------+------------------+-----------------------------+
          Brun (d=3)   100        0.304 (0.0083)     -0.112 (0.0035)    1.369 (0.0068)
          Brun (d=4)   100        0.326 (0.0023)     -0.072 (0.0018)    1.221 (0.0049)
          Brun (d=5)   100        0.309 (0.0010)     -0.046 (0.0012)    1.150 (0.0037)

    """
    rows = []
    for algo in L:
        try:
            row = _lyapunov_row(algo, n_orbits, n_iterations, ncpus=ncpus)
        except Exception as err:
            s = "{}: {}".format(err.__class__.__name__, err)
            print("Ignoring {} in Lyapunov table. {}".format(algo.class_name(), s))
        else:
            rows.append(row)
    rows.sort(key=lambda d:d[4], reverse=True)
    header = ("Algorithm", r"\#Orbits", r"$\theta_1$ (std)", r"$\theta_2$ (std)",
               r"$1-\theta_2/\theta_1$ (std)")
    return table(rows=rows, header_row=header)

