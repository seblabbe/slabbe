r"""
Kolakoski word (datatype)
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
import itertools

cdef class WordDatatype_Kolakoski(object):
    r"""
    Word datatype class for the Kolakoski word.

    This is a cython implementation inspired from the 10 lines of C code
    written by Dominique Bernardi and shared during Sage Days 28 in Orsay,
    France, in January 2011.

    INPUT:

    - ``parent`` - the parent

    EXAMPLES::

        sage: from slabbe.kolakoski_word_pyx import WordDatatype_Kolakoski
        sage: parent = Words([1,2])
        sage: K = WordDatatype_Kolakoski(parent)
        sage: K
        <slabbe.kolakoski_word_pyx.WordDatatype_Kolakoski object at ...>
        sage: K[0]
        1
        sage: K[1]
        2
    """
    cdef public _parent
    cdef public _hash
    cdef public _limit
    def __init__(self, parent):
        r"""
        Constructor. See documentation of WordDatatype_Kolakoski for more
        details.

        EXAMPLES::

            sage: from slabbe.kolakoski_word_pyx import WordDatatype_Kolakoski
            sage: parent = Words([1,2])
            sage: WordDatatype_Kolakoski(parent)
            <slabbe.kolakoski_word_pyx.WordDatatype_Kolakoski object at ...>
        """
        self._parent = parent
        self._limit = 365583569409

    def __getitem__(self, n):
        r"""
        Return the n-th letter of the Kolakoski word.

        This is a cython implementation inspired from the 10 lines of C
        code written by Dominique Bernardi and shared during Sage Days 28
        in Orsay, France, in January 2011.

        INPUT:

        - ``n`` - integer such that `n < 365583569409` or slice

        OUPUT:

        n-th letter of the Kolakoski Word or a substring

        .. NOTE::

            Assuming ``sizeof(unsigned long long)`` is `8`, i.e. 64 bit,
            when the variable ``i`` is equal to ``365583569408`` then ``f``
            is equal to ``2^64 - 1`` so that the line ``g = f + 1`` in the
            loop will bust the capacity.  Hence ``__getitem__`` returns bad
            values beyond ``365583569409``. Maybe the limit is smaller
            because of other operations in the loop. This still has to be
            verified.

        EXAMPLES::

            sage: from slabbe.kolakoski_word_pyx import WordDatatype_Kolakoski
            sage: parent = Words([1,2])
            sage: K = WordDatatype_Kolakoski(parent)
            sage: K
            <slabbe.kolakoski_word_pyx.WordDatatype_Kolakoski object at ...>
            sage: K[0]
            1
            sage: K[1]
            2
            sage: [K[i] for i in range(20)]
            [1, 2, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1]

        Slices works too::

            sage: K[:]
            word: 1221121221221121122121121221121121221221...
            sage: K[:20]
            word: 12211212212211211221
            sage: K[2:20:2]
            word: 211221212
            sage: K[2:20:-2]
            word:
            sage: K[20:2:-2]
            word: 221212211

        TESTS:

        Not yet implemented for too large integers ::

            sage: K[365583569410]
            Traceback (most recent call last):
            ...
            NotImplementedError: when n is larger then 365583569409
        """
        cdef unsigned long long e = 0, f = 0, g, m, i
        if isinstance(n, slice):
            key = n
            from sage.rings.all import Infinity
            from math import ceil
            if not(key.start is None) and key.start < 0 or \
                    not(key.stop is None) and key.stop < 0:
                raise ValueError, "for infinite words, start and stop values cannot be negative"
            step = 1 if key.step is None else key.step
            if step >= 0:
                start = 0 if key.start is None else key.start
                if key.stop is None:
                    length = Infinity
                else: # key.stop > 0
                    length = int(max(0,ceil((key.stop-start)/float(step))))
                data = itertools.islice(self, start, key.stop, step)
                return self._parent(data, length=length)
            else:
                if key.start is None or key.start < 0:
                    raise ValueError, "start value must be nonnegative for negative step values"
                start = key.start
                stop = 0 if key.stop is None else key.stop
                length = int(max(0,ceil((key.stop-start)/float(step))))
                data = list(itertools.islice(self, start+1))[key]
                return self._parent(data, length=length)
        elif n > self._limit:
            raise NotImplementedError("when n is larger then %s" % self._limit)
        else:
            if n == 0:
                e = 1
            elif n == 1:
                e = 0
            else:
                for i from 0 <= i < n-1:
                    g = f + 1
                    m = f ^ g
                    e ^= m / 2
                    f = g + (e & m) / 2
            return int(2 - (e & 1))

    def __iter__(self):
        r"""
        Iterator over the Kolakoski sequence.

        This is a cython implementation inspired from the 10 lines of C
        code written by Dominique Bernardi and shared during Sage Days 28
        in Orsay, France, in January 2011.

        .. NOTE::

            This iterator is defined not further than `365583569409` terms
            if ``sizeof(unsigned long long) == 8`` and less in fact (around
            `10^11`).

        TODO: Define this iterator further then the limit.

        EXAMPLES::

            sage: from slabbe.kolakoski_word_pyx import WordDatatype_Kolakoski
            sage: parent = Words([1,2])
            sage: K = WordDatatype_Kolakoski(parent)
            sage: it = iter(K)
            sage: [it.next() for _ in range(20)]
            [1, 2, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1]
        """
        yield 1
        yield 2
        cdef unsigned long long e = 0, f = 0, g, m
        while True:
            g = f + 1
            m = f ^ g
            e ^= m / 2
            f = g + (e & m) / 2
            yield int(2 - (e & 1))

    def __reduce__(self):
        r"""
        Pickle support

        TESTS::

            sage: from slabbe.kolakoski_word_pyx import WordDatatype_Kolakoski
            sage: parent = Words([1,2])
            sage: K = WordDatatype_Kolakoski(parent)
            sage: K.__reduce__()
            (<type 'slabbe.kolakoski_word_pyx.WordDatatype_Kolakoski'>,
            (Finite and infinite words over {1, 2},))
        """
        return self.__class__, (self._parent, )
