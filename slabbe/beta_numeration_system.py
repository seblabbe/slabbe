# -*- coding: utf-8 -*-
r"""
Beta-numeration

See for instance [Rényi1957]_ or [Parry1960]_.

EXAMPLES::

    sage: from slabbe import BetaTransformation
    sage: z = polygen(QQ, 'z')
    sage: K.<phi> = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
    sage: T = BetaTransformation(phi)
    sage: T.orbit(.1, 5)
    [0.100000000000000,
     0.161803398874990,
     0.261803398874990,
     0.423606797749979,
     0.685410196624969,
     0.109016994374948]
    sage: T.greedy_beta_expansion(.1, 10)
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]

REFERENCES:

.. [Rényi1957] Rényi, A. « Representations for real numbers and their ergodic
   properties ». Acta Mathematica. Academiae Scientiarum Hungaricae 8
   (1957): 477–493. https://doi.org/10.1007/BF02020331.

.. [Parry1960] Parry, W. « On the $\beta$-expansions of real numbers ». Acta
   Mathematica. Academiae Scientiarum Hungaricae 11 (1960): 401–416.
   https://doi.org/10.1007/BF02020954.

AUTHOR:

    - Sébastien Labbé, September 8, 2020
"""
#*****************************************************************************
#       Copyright (C) 2020 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.functions.other import floor

class BetaTransformation(object):
    r"""
    INPUT:

    - ``beta`` -- real number

    EXAMPLES::

        sage: from slabbe import BetaTransformation
        sage: z = polygen(QQ, 'z')
        sage: K.<phi> = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
        sage: BetaTransformation(phi)
        β-transformation with β = phi ≈ 1.61803398874989
    """
    def __init__(self, beta):
        r"""
        EXAMPLES::

            sage: from slabbe import BetaTransformation
            sage: z = polygen(QQ, 'z')
            sage: K.<phi> = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: BetaTransformation(4)
            β-transformation with β = 4 ≈ 4.00000000000000

        """
        self._beta = beta

    def __repr__(self):
        r"""
        Return the string representation.

        EXAMPLES::

            sage: from slabbe import BetaTransformation
            sage: z = polygen(QQ, 'z')
            sage: K.<phi> = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: BetaTransformation(phi)
            β-transformation with β = phi ≈ 1.61803398874989

        """
        return f"β-transformation with β = {self._beta} ≈ {self._beta.n()}"

    def __call__(self, x):
        r"""
        Return the image under the beta-transformation

        INPUT:

        - ``x`` -- real number in [0,1]

        EXAMPLES::

            sage: from slabbe import BetaTransformation
            sage: z = polygen(QQ, 'z')
            sage: K.<phi> = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: T = BetaTransformation(phi)
            sage: T(.1)
            0.161803398874990
            sage: T(.4)
            0.647213595499958
            sage: T(.9)
            0.456230589874906

        """
        beta_x = self._beta * x
        return beta_x - floor(beta_x)

    def orbit(self, x, n):
        r"""
        Return the orbit of length n under the beta transformation starting
        at x.

        INPUT:

        - ``x`` -- real number in [0,1]
        - ``n`` -- postive integer

        OUTPUT:

        list

        EXAMPLES::

            sage: from slabbe import BetaTransformation
            sage: z = polygen(QQ, 'z')
            sage: K.<phi> = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: T = BetaTransformation(phi)
            sage: T.orbit(.1, 5)
            [0.100000000000000,
             0.161803398874990,
             0.261803398874990,
             0.423606797749979,
             0.685410196624969,
             0.109016994374948]
            sage: T.orbit(pi-3, 5)
            [pi - 3,
             1/2*(pi - 3)*(sqrt(5) + 1),
             1/4*(pi - 3)*(sqrt(5) + 1)^2,
             1/8*(pi - 3)*(sqrt(5) + 1)^3,
             1/16*(pi - 3)*(sqrt(5) + 1)^4,
             1/32*(pi - 3)*(sqrt(5) + 1)^5 - 1]

        """
        L = [x]
        for _ in range(n):
            L.append(self(L[-1]))
        return L

    def plot_me(self):
        r"""
        Return the plot of the beta-transformation on interval [0,1].

        EXAMPLES::

            sage: from slabbe import BetaTransformation
            sage: z = polygen(QQ, 'z')
            sage: K.<phi> = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: T = BetaTransformation(phi)
            sage: _ = T.plot_me()

        """
        from sage.plot.plot import plot
        return plot(self, xmin=0, xmax=1)

    def greedy_beta_expansion(self, x, n):
        r"""
        Return the beta expansion of ``x``.

        INPUT:

        - ``x`` -- real number in [0,1]
        - ``n`` -- positive integer

        OUTPUT:

        list

        EXAMPLES::

            sage: from slabbe import BetaTransformation
            sage: z = polygen(QQ, 'z')
            sage: K.<phi> = NumberField(z**2-z-1, 'phi', embedding=RR(1.6))
            sage: T = BetaTransformation(phi)
            sage: T.greedy_beta_expansion(.1, 10)
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]
            sage: T.greedy_beta_expansion(pi-3, 10)
            [0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0]

        """
        L = self.orbit(x, n)
        return [floor(self._beta*a) for a in L]

