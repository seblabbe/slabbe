# -*- coding: utf-8 -*-
r"""
The code to construct the partitions in [Lab2019]_.

REFERENCES:

.. [Lab2019] S. Labbé.  Markov partitions for toral
    `\mathbb{Z}^2`-rotations featuring Jeandel-Rao Wang shift and model
    sets, https://arxiv.org/abs/1903.06137, April 2020 (v3).

.. [Lab2018] S. Labbé. A self-similar aperiodic set of 19 Wang
    tiles. Geom. Dedicata, 2018.
    https://doi.org/10.1007/s10711-018-0384-8.

EXAMPLES:

The partition associated to Jeandel-Rao Wang shift::

    sage: from slabbe.arXiv_1903_06137 import jeandel_rao_wang_shift_partition
    sage: P0 = jeandel_rao_wang_shift_partition()
    sage: P0
    Polyhedron partition of 24 atoms with 11 letters

The partition associated to the self-similar Wang shift `\Omega_{\mathcal{U}}`::

    sage: from slabbe.arXiv_1903_06137 import self_similar_19_atoms_partition
    sage: PU = self_similar_19_atoms_partition()
    sage: PU
    Polyhedron partition of 19 atoms with 19 letters

"""
#*****************************************************************************
#       Copyright (C) 2019-2020 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field import NumberField
from sage.geometry.polyhedron.constructor import Polyhedron

from slabbe import PolyhedronPartition

def jeandel_rao_wang_shift_partition(backend=None):
    r"""
    This construct the polygon partition associated to Jeandel-Rao
    tilings introduced in [Lab2019]_.

    INPUT:

    - ``backend`` -- string, polyhedron backend

    EXAMPLES::

        sage: from slabbe.arXiv_1903_06137 import jeandel_rao_wang_shift_partition
        sage: P0 = jeandel_rao_wang_shift_partition()
        sage: P0.is_pairwise_disjoint()
        True
        sage: P0.volume()
        4*phi + 1

    The volume is consistent with::

        sage: z = polygen(QQ, 'z')
        sage: K = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
        sage: phi = K.gen()
        sage: phi * (phi + 3)
        4*phi + 1

    """
    # the golden mean
    z = polygen(QQ, 'z')
    K = NumberField(z**2-z-1, 'phi', embedding=QQ(1.6))
    phi = K.gen()

    # x and y coordinates
    xcoords = [0, phi**-2, phi**-1, 1, phi]
    ycoords = [0, 1, 2, phi+1, phi+2, phi+3]

    # the vertices
    A = [(x,0) for x in xcoords]
    B = [(x,1) for x in xcoords]
    C = [(x,2) for x in xcoords]
    D = [(x,phi+1) for x in xcoords]
    E = [(x,phi+2) for x in xcoords]
    F = [(x,phi+3) for x in xcoords]

    # atoms corresponding to Jeandel-Rao tile numbers
    L = [
        (0, (A[0], A[2], B[2])), 
        (0, (A[2], A[3], B[3])),
        (0, (A[3], A[4], B[4])),
        (1, (A[0], B[0], B[2])),
        (1, (A[2], B[2], B[3])),
        (1, (A[3], B[3], B[4])),
        (2, (C[0], D[0], E[2])),
        (3, (B[3], C[3], E[4], C[4])),
        (4, (C[0], E[2], E[3])),
        (4, (E[2], E[3], F[3])),
        (5, (D[0], E[0], E[2])),
        (5, (E[0], E[1], F[1])),
        (5, (E[1], E[2], F[3])),
        (6, (E[0], F[0], F[1])),
        (6, (E[1], F[1], F[3])),
        (6, (E[3], F[3], F[4])),
        (7, (D[3], E[3], E[4])),
        (7, (E[3], E[4], F[4])),
        (7, (B[0], C[0], E[3])),
        (8, (B[2], C[2], E[4])),
        (9, (B[0], B[2], C[2])),
        (9, (B[2], B[3], C[3])),
        (9, (B[3], B[4], C[4])),
        (10, (B[0], D[3], E[3])),
        ]
    L = [(key, Polyhedron(vertices, base_ring=K, backend=backend)) 
         for (key,vertices) in L]
    return PolyhedronPartition(L)

def self_similar_19_atoms_partition(backend=None):
    r"""
    This construct the polygon partition introduced in [Lab2019]_
    associated to the self-similar 19 Wang tiles [Lab2018]_.

    INPUT:

    - ``backend`` -- string, polyhedron backend

    EXAMPLES::

        sage: from slabbe.arXiv_1903_06137 import self_similar_19_atoms_partition
        sage: PU = self_similar_19_atoms_partition()
        sage: PU.is_pairwise_disjoint()
        True
        sage: PU.volume()
        1

    """
    # the golden mean
    z = polygen(QQ, 'z')
    K = NumberField(z**2-z-1, 'phi', embedding=QQ(1.6))
    phi = K.gen()

    # the partition vertices
    L = [
        (0, [(phi - 1, -2*phi + 4), (phi - 1, phi - 1), (-2*phi + 4, phi - 1)]),
        (1,
        [(phi - 1, -2*phi + 4),
        (phi - 1, 1),
        (1, 1),
        (-2*phi + 4, phi - 1),
        (1, phi - 1)]),
        (2, [(0, 1), (-phi + 2, phi - 1), (2*phi - 3, phi - 1)]),
        (3, [(0, 1), (0, phi - 1), (2*phi - 3, phi - 1)]),
        (4, [(phi - 1, phi - 1), (-phi + 2, 1), (phi - 1, -2*phi + 4)]),
        (5, [(phi - 1, 1), (-phi + 2, 1), (phi - 1, -2*phi + 4)]),
        (6, [(0, 1), (-phi + 2, 1), (-phi + 2, phi - 1)]),
        (7, [(-phi + 2, 1), (-phi + 2, phi - 1), (phi - 1, phi - 1)]),
        (8, [(1, 0), (phi - 1, 0), (phi - 1, -phi + 2)]),
        (9, [(-2*phi + 4, phi - 1), (phi - 1, phi - 1), (1, -phi + 2), (1, 0)]),
        (10, [(1, -phi + 2), (-2*phi + 4, phi - 1), (1, phi - 1)]),
        (11, [(1, 0), (phi - 1, -phi + 2), (phi - 1, phi - 1)]),
        (12, [(-phi + 2, -phi + 2), (-phi + 2, phi - 1), (2*phi - 3, phi - 1)]),
        (13,
        [(-phi + 2, phi - 1),
        (phi - 1, -phi + 2),
        (-phi + 2, -phi + 2),
        (phi - 1, 0)]),
        (14,
        [(2*phi - 3, phi - 1), (0, phi - 1), (-phi + 2, -phi + 2), (-phi + 2, 0)]),
        (15, [(phi - 1, 0), (-phi + 2, 0), (-phi + 2, -phi + 2)]),
        (16, [(0, 0), (-phi + 2, 0), (0, -phi + 2)]),
        (17, [(0, -phi + 2), (-phi + 2, 0), (0, phi - 1)]),
        (18, [(phi - 1, -phi + 2), (-phi + 2, phi - 1), (phi - 1, phi - 1)])
        ]

    L = [(key, Polyhedron(vertices, base_ring=K, backend=backend)) 
         for (key,vertices) in L]
    return PolyhedronPartition(L)

