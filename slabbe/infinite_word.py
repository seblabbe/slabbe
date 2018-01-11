# -*- coding: utf-8 -*-
r"""
Infinite words

Methods that are not in Sage (for now!)

AUTHORS:

- Sébastien Labbé, 2016

EXAMPLES::

    sage: from slabbe.infinite_word import derived_sequence
    sage: w = words.ThueMorseWord()
    sage: derived_sequence(w, w[:1])
    word: 0120210121020120210201210120210121020121...
"""
#*****************************************************************************
#       Copyright (C) 2016 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
from sage.combinat.words.words import InfiniteWords
from sage.rings.semirings.non_negative_integer_semiring import NN

def derived_sequence(self, u, coding=False):
    r"""
    Return the derived sequence of according to the return words to a factor of
    self.

    INPUT:

    - ``u`` -- finite word, the length of the nonempty prefix
    - ``coding`` -- boolean (default: ``False``), whether to
      include the return word coding dictionnary

    EXAMPLES::

        sage: from slabbe.infinite_word import derived_sequence
        sage: w = words.ThueMorseWord()
        sage: derived_sequence(w, w[:1])
        word: 0120210121020120210201210120210121020121...
        sage: derived_sequence(w, w[:2])
        word: 0123013201232013012301320130123201230132...
        sage: derived_sequence(w, w[:3])
        word: 0123013201232013012301320130123201230132...

    With the return word coding::

        sage: w = words.ThueMorseWord()
        sage: derived, D = derived_sequence(w, w[:1], True)
        sage: derived
        word: 0120210121020120210201210120210121020121...
        sage: D
        {word: 0: 2, word: 01: 1, word: 011: 0}

    It gets into a cycle of length 1::

        sage: words.ThueMorseWord()
        word: 0110100110010110100101100110100110010110...
        sage: derived_sequence(_, _[:1])
        word: 0120210121020120210201210120210121020121...
        sage: derived_sequence(_, _[:1])
        word: 0123013201232013012301320130123201230132...
        sage: derived_sequence(_, _[:1])
        word: 0123013201232013012301320130123201230132...

    .. NOTE::

        Note that method ``return_words_derivate`` of finite words in Sage does
        the same for finite words but without returning the translation
        dictionnary::

            sage: w = words.ThueMorseWord()
            sage: prefix = w[:1000]
            sage: prefix.return_words_derivate(prefix[:1])
            word: 1231321232131231321312321231321232131232...
    """
    D = {}
    it = (D.setdefault(w, len(D)) for w in self.return_words_iterator(u))
    W = InfiniteWords(alphabet=NN)
    w = W(it, datatype='iter')
    if coding:
        return w, D
    else:
        return w

