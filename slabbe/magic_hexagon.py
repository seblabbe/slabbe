# -*- coding: utf-8 -*-
r"""
Magic hexagon solver

We reduce the problem to an integer linear program.

EXAMPLES:

It allows to compare the efficiency of LP solvers::

    sage: from slabbe.magic_hexagon import solve_magic_hexagon
    sage: _ = solve_magic_hexagon(solver='Coin')   # not tested (>5min?)
    sage: _ = solve_magic_hexagon(solver='GLPK')   # not tested (90s)
    sage: _ = solve_magic_hexagon(solver='Gurobi') # not tested (<1s)

See https://en.wikipedia.org/wiki/Magic_hexagon
"""
#*****************************************************************************
#       Copyright (C) 2018 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
from sage.numerical.mip import MixedIntegerLinearProgram

def solve_magic_hexagon(solver=None):
    r"""
    Solves the magic hexagon problem

    We use the following convention for the positions::

           1   2  3
         4  5   6   7
       8   9  10  11  12
        13  14  15  16
          17  18  19

    INPUT:

    - ``solver`` -- string (default:``None``)

    EXAMPLES::

        sage: from slabbe.magic_hexagon import solve_magic_hexagon
        sage: solve_magic_hexagon() # long time (90s if GLPK, <1s if Gurobi)
        [15, 14, 9, 13, 8, 6, 11, 10, 4, 5, 1, 18, 12, 2, 7, 17, 16, 19, 3]

    """
    p = MixedIntegerLinearProgram(solver=solver)
    x = p.new_variable(binary=True)

    A = range(1,20)

    # exactly one tile at each position pos
    for pos in A:
        S = p.sum(x[pos,tile] for tile in A)
        name = "one tile at {}".format(pos)
        p.add_constraint(S==1, name=name)

    # each tile used exactly once
    for tile in A:
        S = p.sum(x[pos,tile] for pos in A)
        name = "tile {} used once".format(tile)
        p.add_constraint(S==1, name=name)

    lines = [(1,2,3), (4,5,6,7), (8,9,10,11,12), (13,14,15,16), (17,18,19),
             (1,4,8), (2,5,9,13), (3,6,10,14,17), (7,11,15,18), (12,16,19),
             (8,13,17), (4,9,14,18), (1,5,10,15,19), (2,6,11,16), (3,7,12)]

    # the sum is 38 on each line
    for line in lines:
        S = p.sum(tile*x[pos,tile] for tile in A for pos in line) 
        name = "sum of line {} must be 38".format(line)
        p.add_constraint(S==38, name=name)

    p.solve()
    soln = p.get_values(x)
    nonzero = sorted(key for key in soln if soln[key]!=0)
    return [tile for (pos,tile) in nonzero]


