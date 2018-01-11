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

def lyapunov_sample(algo, n_orbits, n_iterations=1000, verbose=False):
    r"""
    Return lists of values for theta1, theta2 and 1-theta2/theta1 computed
    on many orbits.

    This is computed in parallel.

    INPUT:

    - ``n_orbits`` -- integer, number of orbits
    - ``n_iterations`` -- integer, length of each orbit
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
    """
    S = [(random(), random(), random()) for _ in range(n_orbits)]
    @parallel
    def compute_exponents(i):
        try:
            return algo.lyapunov_exponents(start=S[i], n_iterations=n_iterations) 
        except Exception as err:
            return "{}: {}".format(err.__class__.__name__, err)
    L = [v for _,v in compute_exponents(range(n_orbits))]
    L_filtered = [v for v in L if isinstance(v, tuple)]
    if verbose:
        L_error_msg = [v for v in L if not isinstance(v, tuple)]
        print(L_error_msg)
    return zip(*L_filtered)

def lyapunov_table(algo, n_orbits, n_iterations=1000):
    r"""
    Return a table of values of Lyapunov exponents for this algorithm.

    INPUT:

    - ``n_orbits`` -- integer, number of orbits
    - ``n_iterations`` -- integer, length of each orbit

    OUTPUT:

        table of liapounov exponents

    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Brun
        sage: from slabbe.lyapunov import lyapunov_table
        sage: lyapunov_table(Brun(), 10, 1000000) # random
          10 succesfull orbits    min       mean      max       std
        +-----------------------+---------+---------+---------+---------+
          $\theta_1$              0.303     0.305     0.307     0.0013
          $\theta_2$              -0.1131   -0.1124   -0.1115   0.00051
          $1-\theta_2/\theta_1$   1.3678    1.3687    1.3691    0.00043
    """
    import numpy as np
    from sage.misc.functional import numerical_approx
    from sage.rings.real_mpfr import RR
    from sage.misc.table import table
    rep = lyapunov_sample(algo, n_orbits, n_iterations)
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
    header = ['{} succesfull orbits'.format(len(rep[0])), 'min','mean','max','std']
    return table(rows=rows,header_row=header)

def _lyapunov_row(algo, n_orbits, n_iterations=1000):
    r"""
    Return a row of values of Lyapunov exponents.

    This method is used for making a larger table comparing many
    algorithms.

    INPUT:

    - ``n_orbits`` -- integer, number of orbits
    - ``n_iterations`` -- integer, length of each orbit

    OUTPUT:

        row of liapounov exponents: theta1, theta2, 1-theta2/theta1

    EXAMPLES::

        sage: from slabbe.lyapunov import _lyapunov_row
        sage: from slabbe.mult_cont_frac import Brun
        sage: _lyapunov_row(Brun(), 10, 100000) # abs tol 0.01
        ['Brun', 10, '0.303 (0.0038)', '-0.112 (0.0019)', '1.368 (0.0026)']
    """
    import numpy as np
    from sage.misc.functional import numerical_approx
    from sage.rings.real_mpfr import RR
    rep = lyapunov_sample(algo, n_orbits, n_iterations)
    def floor_log(number):
        return RR(number).abs().log10().floor()
    def my_rounded(number, s):
        m = floor_log(number)
        return numerical_approx(number, digits=m-s+1)
    row = []
    row.append(algo.name())
    row.append(len(rep[0]))
    for data in rep:
        data = np.array(data)
        s = floor_log(data.std())
        val = my_rounded(data.mean(), s)
        std = my_rounded(data.std(), s)
        row.append("{} ({})".format(val, std))
    return row

def lyapunov_comparison_table(L, n_orbits=100, n_iterations=10000):
    r"""
    Return a table of values of Lyapunov exponents for many algorithm.

    INPUT:

    - ``L`` -- list of algorithms
    - ``n_orbits`` -- integer
    - ``n_iterations`` -- integer

    OUTPUT:

        table

    EXAMPLES::

        sage: import slabbe.mult_cont_frac as mcf
        sage: from slabbe.lyapunov import lyapunov_comparison_table
        sage: algos = [mcf.Brun(), mcf.ARP()]
        sage: lyapunov_comparison_table(algos)    # abs tol 0.01
          Algorithm                 \#Orbits   $\theta_1$ (std)   $\theta_2$ (std)   $1-\theta_2/\theta_1$ (std)
        +-------------------------+----------+------------------+------------------+-----------------------------+
          Arnoux-Rauzy-Poincar\'e   100        0.44 (0.012)       -0.172 (0.0060)    1.388 (0.0054)
          Brun                      100        0.30 (0.011)       -0.113 (0.0049)    1.370 (0.0070)
    """
    rows = []
    for algo in L:
        try:
            row = _lyapunov_row(algo, n_orbits, n_iterations)
        except Exception as err:
            s = "{}: {}".format(err.__class__.__name__, err)
            print("Ignoring {} in Lyapunov table. {}".format(algo.class_name(), s))
        else:
            rows.append(row)
    rows.sort(key=lambda d:d[4], reverse=True)
    header = ("Algorithm", r"\#Orbits", r"$\theta_1$ (std)", r"$\theta_2$ (std)",
               r"$1-\theta_2/\theta_1$ (std)")
    return table(rows=rows, header_row=header)

