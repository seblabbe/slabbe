# -*- coding: utf-8 -*-
r"""
The code contained in `arXiv:1906.01104`__

S. Labbé, Rauzy induction of polygon partitions and toral $\Z^2$-rotations
arXiv:1906.01104, June 2019, 36 p. (v2: revised in May 2020).

__ https://arxiv.org/abs/1906.01104

EXAMPLES::

    sage: from slabbe.arXiv_1906_01104 import beta8,beta9,tau
    sage: beta8*beta9*tau
    Substitution 2d: {0: [[17]], 1: [[12]], 2: [[16, 10]], 3: [[16, 9]], 
    4: [[17, 7]], 5: [[12, 7]], 6: [[16], [2]], 7: [[14], [4]], 
    8: [[17], [2]], 9: [[13], [3]], 10: [[13], [2]], 11: [[12], [2]], 
    12: [[15, 11], [5, 1]], 13: [[18, 10], [4, 1]], 14: [[16, 10], [3, 1]],
    15: [[16, 9], [2, 0]], 16: [[14, 6], [4, 1]], 17: [[14, 8], [4, 1]], 
    18: [[13, 6], [3, 1]]}

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
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.geometry.polyhedron.library import polytopes

# Polyhedron exchange transformations and induction of polyhedron partitions are
# implemented in Sage with the optional package slabbe
from slabbe import PolyhedronExchangeTransformation as PET
from slabbe import Substitution2d

# We construct the golden mean as a element of a quadratic number field
# because it is more efficient for arithmetic operations and comparisons:
z = polygen(QQ, 'z')
K = NumberField(z**2-z-1, 'phi', embedding=QQ(1.6))
phi = K.gen()

# The initial torus and Z^2-actions
Gamma0 = matrix.column([(phi,0), (1,phi+3)])
fundamental_domain = polytopes.parallelotope([(phi,0), (0,phi+3)])
R0e1 = PET.toral_translation(Gamma0, vector((1,0)), fundamental_domain)
R0e2 = PET.toral_translation(Gamma0, vector((0,1)), fundamental_domain)

# The initial partition P0
from slabbe.arXiv_1903_06137 import jeandel_rao_wang_shift_partition
P0 = jeandel_rao_wang_shift_partition()

# Compute P1 and beta0
y_le_1 = [1, 0, -1]   # syntax for the inequality y <= 1
P1,beta0 = R0e2.induced_partition(y_le_1, P0, substitution_type='column')
R1e1,_ = R0e1.induced_transformation(y_le_1)
R1e2,_ = R0e2.induced_transformation(y_le_1)

# The lattices (we don't need them, it's just to show how to construct them)
Gamma2 = Gamma1 = matrix.column([(phi,0), (0,1)])

# Change of base (shearing the action)
P2 = P1
R2e1 = R1e1
R2e2 = (R1e1 * R1e2).merge_atoms_with_same_translation()

# Compute P3 and beta2
x_le_1 = [1, -1, 0]   # syntax for x <= 1
P3,beta2 = R2e1.induced_partition(x_le_1, P2, substitution_type='row')
R3e1,_ = R2e1.induced_transformation(x_le_1)
R3e2,_ = R2e2.induced_transformation(x_le_1)

# Compute P4 and beta3
x_le_phi_inv = [phi**-1, -1, 0]    # syntax for x <= phi^-1
P4,beta3 = R3e1.induced_partition(x_le_phi_inv, P3, substitution_type='row')
R4e1,_ = R3e1.induced_transformation(x_le_phi_inv)
R4e2,_ = R3e2.induced_transformation(x_le_phi_inv)

# Compute P5 and beta4
y_le_phi_inv = [phi**-1, 0, -1]    # syntax for y <= phi^-1
P5,beta4 = R4e2.induced_partition(y_le_phi_inv, P4, substitution_type='column')
R5e1,_ = R4e1.induced_transformation(y_le_phi_inv)
R5e2,_ = R4e2.induced_transformation(y_le_phi_inv)

# Renormalization
P5_scaled = (-phi*P5).translate((1,1))
R5e1_scaled = (-phi*R5e1).translate_domain((1,1))
R5e2_scaled = (-phi*R5e2).translate_domain((1,1))

# Compute P6 and beta5
P6,beta5 = R5e1_scaled.induced_partition(x_le_phi_inv, P5_scaled, substitution_type='row')
R6e1,_ = R5e1_scaled.induced_transformation(x_le_phi_inv)
R6e2,_ = R5e2_scaled.induced_transformation(x_le_phi_inv)

# Compute P7 and beta6
P7,beta6 = R6e2.induced_partition(y_le_phi_inv, P6, substitution_type='column')
R7e1,_ = R6e1.induced_transformation(y_le_phi_inv)
R7e2,_ = R6e2.induced_transformation(y_le_phi_inv)

# Renormalization
P7_scaled = (-phi*P7).translate((1,1))
R7e1_scaled = (-phi*R7e1).translate_domain((1,1))
R7e2_scaled = (-phi*R7e2).translate_domain((1,1))

# Compute P8 and beta7
P8,beta7 = R7e1_scaled.induced_partition(x_le_phi_inv, P7_scaled, substitution_type='row')
R8e1,_ = R7e1_scaled.induced_transformation(x_le_phi_inv)
R8e2,_ = R7e2_scaled.induced_transformation(x_le_phi_inv)

# Compute P9 and beta8
P9,beta8 = R8e2.induced_partition(y_le_phi_inv, P8, substitution_type='column')
R9e1,_ = R8e1.induced_transformation(y_le_phi_inv)
R9e2,_ = R8e2.induced_transformation(y_le_phi_inv)

# Renormalization
P9_scaled = (-phi*P9).translate((1,1))
R9e1_scaled = (-phi*R9e1).translate_domain((1,1))
R9e2_scaled = (-phi*R9e2).translate_domain((1,1))

# Compute P10 and beta9
P10,beta9 = R9e1_scaled.induced_partition(x_le_phi_inv, P9_scaled, substitution_type='row')
R10e1,_ = R9e1_scaled.induced_transformation(x_le_phi_inv)
R10e2,_ = R9e2_scaled.induced_transformation(x_le_phi_inv)

# the self-similarity for P8
assert P8 == P10
tau = Substitution2d.from_permutation(P8.keys_permutation(P10))
self_similarity_P8 = beta8*beta9*tau

# We may check that the self-similarity for $\Pcal_8$ satisfies 
# $\zeta^{-1}\beta_8\beta_9\tau\zeta=\beta_{\mathcal{U}}$.

zeta = Substitution2d.from_permutation({0:0, 1:1, 2:9, 3:7, 4:8, 5:11, 6:10, 
7:6, 8:2, 9:4, 10:5, 11:3, 12:18, 13:14, 14:16, 15:13, 16:12, 17:17, 18:15})

betaU = Substitution2d({0: [[17]], 1: [[16]], 2: [[15], [11]], 3: [[13], [9]], 4: [[17], [8]], 5: [[16], [8]], 6: [[15], [8]], 7: [[14], [8]], 8: [[14, 6]], 9: [[17, 3]], 10: [[16, 3]], 11: [[14, 2]], 12: [[15, 7], [11, 1]], 13: [[14, 6], [11, 1]], 14: [[13, 7], [9, 1]], 15: [[12, 6], [9, 1]], 16: [[18, 5], [10, 1]], 17: [[13, 4], [9, 1]], 18: [[14, 2], [8, 0]]})

self_sim_PU_deduced = zeta.inverse()*beta8*beta9*tau*zeta
assert betaU == self_sim_PU_deduced 

