# -*- coding: utf-8 -*-
r"""
Kolakoski Word

The classical Kolakoski sequence [K65]_ was first studied by Oldenburger
[O39]_, where it appears as the unique solution to the problem of a
trajectory on the alphabet `\{1,2\}` which is identical to its exponent
trajectory.

See http://en.wikipedia.org/wiki/Kolakoski_sequence

It uses cython implementation inspired from the 10 lines of C code written
by Dominique Bernardi and shared during Sage Days 28 in Orsay, France, in
January 2011. The alphabet must be ``(1,2)``.

EXAMPLES::

    sage: from slabbe import KolakoskiWord
    sage: K = KolakoskiWord()
    sage: K
    word: 1221121221221121122121121221121121221221...

The cython implementation is much faster than the python one::

    sage: K = KolakoskiWord()
    sage: K[10^6]      # takes 0.02 seconds
    2
    sage: K = words.KolakoskiWord()
    sage: K[10^6]      # not tested : takes too long
    2

TESTS:

We make sure both implementation correspond for the prefix of
length 100::

    sage: Kb = KolakoskiWord()
    sage: Kp = words.KolakoskiWord()
    sage: Kp[:100] == Kb[:100]
    True

REFERENCES:

.. [K65] William Kolakoski, proposal 5304, American Mathematical Monthly
   72 (1965), 674; for a partial solution, see "Self Generating Runs,"
   by Necdet Üçoluk, Amer. Math. Mon. 73 (1966), 681-2.

.. [O39] R. Oldenburger, Exponent trajectories in dynamical systems,
   Trans. Amer. Math. Soc. 46 (1939), 453–466.

AUTHORS:

- Sébastien Labbé (February 2011) - Added a Cython implementation
  which was easily translated from the 10 lines of C code written
  by Dominique Bernardi and shared during Sage Days 28 in Orsay.
  Needs review at #13346 since too long time...
"""
#*****************************************************************************
#       Copyright (C) 2011-2014 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function
from sage.combinat.words.infinite_word import InfiniteWord_class
from .kolakoski_word_pyx import WordDatatype_Kolakoski
 
class KolakoskiWord(WordDatatype_Kolakoski, InfiniteWord_class):
    def __init__(self, parent=None):
        r"""
        Constructor. See documentation of WordDatatype_Kolakoski for more
        details.

        EXAMPLES::

            sage: from slabbe import KolakoskiWord
            sage: K = KolakoskiWord()
            sage: K
            word: 1221121221221121122121121221121121221221...

        TESTS:

        Pickle is supported::

            sage: K = KolakoskiWord()
            sage: loads(dumps(K))
            word: 1221121221221121122121121221121121221221...
        """
        if parent is None:
            from sage.combinat.words.words import Words
            parent = Words([1,2])
        WordDatatype_Kolakoski.__init__(self, parent)

