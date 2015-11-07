# -*- coding: utf-8 -*-
r"""
Comparison of Multidimensional Continued Fraction Algorithms

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
from sage.misc.table import table

def lyapunov_exponents_table(L, n_orbits=100, n_iterations=10000):
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
        sage: algos = [mcf.Brun(), mcf.ARP()]
        sage: compare_algos_for_lyapunov(algos)
          Algorithm   $\theta_1$ (std)   $\theta_2$ (std)   $1-\theta_2/\theta_1$
          ARP         0.44 (0.013)       -0.172 (0.0057)    1.388 (0.0048)
          Brun        0.30 (0.017)       -0.111 (0.0062)    1.368 (0.0068)
    """
    rows = []
    for algo in L:
        try:
            row = algo._lyapunov_exponents_row(n_orbits, n_iterations)
        except Exception as err:
            s = "{}: {}".format(err.__class__.__name__, err)
            print "Ignoring {} in Lyapunov table. {}".format(algo.name(), s)
        else:
            rows.append(row)
    rows.sort(key=lambda d:d[3], reverse=True)
    header = ("Algorithm", r"$\theta_1$ (std)", r"$\theta_2$ (std)",
               r"$1-\theta_2/\theta_1$ (std)")
    return table(rows=rows, header_row=header)

