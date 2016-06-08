# -*- coding: utf-8 -*-
r"""
Infinite words

Methods that are not in Sage (for now!)

AUTHORS:

- Sébastien Labbé, 2016

EXAMPLES::

    sage: from slabbe.infinite_word import derived_sequence
    sage: w = words.ThueMorseWord()
    sage: derived_sequence(w, 1)
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
import itertools
import collections
from sage.combinat.words.words import InfiniteWords
from sage.rings.semirings.non_negative_integer_semiring import NN

def derived_sequence(self, k, include_translation_dict=False):
    r"""
    Return the derived sequence of according to the return words to the prefix
    of length k of self.

    INPUT:

    - ``k`` -- integer, the length of the nonempty prefix
    - ``include_translation_dict`` -- boolean (default: ``False``), whether to
      include translation dict

    EXAMPLES::

        sage: from slabbe.infinite_word import derived_sequence
        sage: w = words.ThueMorseWord()
        sage: derived_sequence(w, 1)
        word: 0120210121020120210201210120210121020121...
        sage: derived_sequence(w, 2)
        word: 0123013201232013012301320130123201230132...
        sage: derived_sequence(w, 3)
        word: 0123013201232013012301320130123201230132...

    With the translation dict::

        sage: w = words.ThueMorseWord()
        sage: derived, D = derived_sequence(w, 1, True)
        sage: derived
        word: 0120210121020120210201210120210121020121...
        sage: D
        defaultdict(<function <lambda> at ...>, {word: 0: 2, word: 01: 1, word: 011: 0})

    It gets into a cycle of length 1::

        sage: words.ThueMorseWord()
        word: 0110100110010110100101100110100110010110...
        sage: derived_sequence(_, 1)
        word: 0120210121020120210201210120210121020121...
        sage: derived_sequence(_, 1)
        word: 0123013201232013012301320130123201230132...
        sage: derived_sequence(_, 1)
        word: 0123013201232013012301320130123201230132...
    """
    u = self[:k]
    c = itertools.count()
    D = collections.defaultdict(lambda:next(c))
    it = (D[w] for w in self.return_words_iterator(u))
    W = InfiniteWords(alphabet=NN)
    w = W(it, datatype='iter')
    if include_translation_dict:
        return w, D
    else:
        return w
        
