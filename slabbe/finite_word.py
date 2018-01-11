# -*- coding: utf-8 -*-
r"""
Finite words

Methods that are not in Sage (for now!)

AUTHORS:

- Sébastien Labbé, 2015
- Sébastien Labbé, 2017, added lexicographic Lyndon stuff

EXAMPLES::

    sage: from slabbe.finite_word import discrepancy
    sage: w = words.ChristoffelWord(5,8)
    sage: discrepancy(w)
    12/13
"""
#*****************************************************************************
#       Copyright (C) 2015,2017 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
import itertools
from collections import Counter
from random import randrange
from sage.rings.rational_field import QQ
from sage.functions.other import abs

def discrepancy(self):
    r"""
    Return the discrepancy of the word.
    
    This is a distance to the euclidean line defined in [T1980]_.

    EXAMPLES::

        sage: from slabbe.finite_word import discrepancy
        sage: w = words.ChristoffelWord(5,8)
        sage: w
        word: 0010010100101
        sage: discrepancy(w)
        12/13

    ::

        sage: for c in w.conjugates(): print (c, discrepancy(c))
        0010010100101 12/13
        0100101001010 7/13
        1001010010100 10/13
        0010100101001 10/13
        0101001010010 7/13
        1010010100100 12/13
        0100101001001 8/13
        1001010010010 9/13
        0010100100101 11/13
        0101001001010 6/13
        1010010010100 11/13
        0100100101001 9/13
        1001001010010 8/13

    REFERENCES:

    .. [T1980] R., Tijdeman. The chairman assignment problem. Discrete
       Mathematics 32, no 3 (1980): 323-30. doi:10.1016/0012-365X(80)90269-1.
    """
    length = self.length()
    freq = {a:QQ((v,length)) for (a,v) in self.evaluation_dict().iteritems()}
    C = Counter()
    M = 0
    for i,a in enumerate(self):
        C[a] += 1
        MM = max(abs(freq[b]*(i+1) - C[b]) for b in freq)
        M = max(M, MM)
    return M

def to_image(self, width=1000):
    r"""
    Creates an image from a word

    INPUT:

    - ``width`` -- integer, width of image

    EXAMPLES::

        sage: from slabbe.finite_word import to_image
        sage: t = words.ThueMorseWord()
        sage: img = to_image(t[:10000], width=100)
        sage: img
        <PIL.Image.Image image mode=RGB size=100x100 at 0x...>
        sage: img.show()    # not tested

    ::

        sage: W = FiniteWords(range(10))
        sage: d = {a:W.random_element(7) for a in range(10)}
        sage: m = WordMorphism(d, codomain=W)
        sage: w = m.periodic_points()[0][0]

    ::

        sage: s = map(int, str(pi.n(digits=40001))[2:])
        sage: len(s)
        40000
        sage: img = to_image(W(s), 200)
        sage: img.show()    # not tested
    """
    #http://stackoverflow.com/questions/434583/what-is-the-fastest-way-to-draw-an-image-from-discrete-pixel-values-in-python
    import numpy as np
    import scipy.misc as smp

    height = self.length() // width

    # Create a 1024x1024x3 array of 8 bit unsigned integers
    data = np.zeros( (height,width,3), dtype=np.uint8 )
    data += 255   # white as default color
    alphabet = self.parent().alphabet()

    color_dict={a:[randrange(256) for _ in range(3)] for a in alphabet}

    for i,a in enumerate(self):
        x = i % width
        y = i // width
        data[y,x] = color_dict[a]

    img = smp.toimage( data )       # Create a PIL image
    #img.show()                     # View in default viewer
    return img

def run_length_encoding(self):
    r"""
    EXAMPLES::

        sage: from slabbe.finite_word import run_length_encoding
        sage: L = [0, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3]
        sage: run_length_encoding(L)
        [(0, 1), (1, 6), (2, 1), (1, 1), (2, 5), (3, 6)]
    """
    return [(value, sum(1 for _ in it)) for value,it in itertools.groupby(self)]

def word_to_polyomino(self):
    r"""
    Returns the inside points of a polyomino.

    INPUT:

    - ``self`` -- list of integers in 0,1,2,3 describing a closed path

    OUTPUT:

    - list of 2d vectors

    EXAMPLES::

        sage: from slabbe.finite_word import word_to_polyomino
        sage: w = [0,0,0,1,1,1,2,2,2,3,3,3]
        sage: word_to_polyomino(w)
        [(0, 0), (0, 1), (1, 0), (2, 0), (1, 1), (0, 2), (1, 2), (2, 1), (2, 2)]
        sage: w = [1,1,1,0,0,0,3,3,3,2,2,2]
        sage: word_to_polyomino(w)
        [(0, 0), (0, 1), (1, 0), (2, 0), (1, 1), (0, 2), (1, 2), (2, 1), (2, 2)]
    """
    from sage.combinat.words.words import FiniteWords
    alphabet = [0,1,2,3]
    FW = FiniteWords(alphabet)
    w = FW(self)

    # Checking it is closed
    end = (w.count(0)-w.count(2), w.count(1) - w.count(3))
    if end != (0,0):
        raise ValueError("input word is not closed path, "
                "it finishes at {}".format(end))

    # Checking its turning number is 1 = 4/4
    w_ = w + FW([self[0]])
    turns = w_.finite_differences(mod=4)
    turning_number = turns.count(1) - turns.count(3)
    if turning_number == 4:
        pass
    elif turning_number == -4:
        self = [(a+2)%4 for a in reversed(self)]
    else:
        raise ValueError("input word is not simple path, "
                "its turning number is {}".format(turning_number))

    # Setting up directions, outside and inside points of each edge
    from sage.modules.free_module import FreeModule
    from sage.rings.integer_ring import ZZ
    F = FreeModule(ZZ, 2)
    directions = map(F, [(1,0), (0,1), (-1,0), (0,-1)])
    outside = map(F, [(0,-1), (0,0), (-1,0), (-1,-1)])
    inside = map(F, [(0,0), (-1,0), (-1,-1), (0,-1)])
    # obtained from below after a translation by (-1/2,-1/2)
    #outside = map(F, [(u/2,-u/2), (u/2,u/2), (-u/2,u/2), (-u/2,-u/2)])
    #inside = map(F, [(u/2,u/2), (-u/2,u/2), (-u/2,-u/2), (u/2,-u/2)])
    for p in directions: p.set_immutable()
    for p in outside: p.set_immutable()
    for p in inside: p.set_immutable()
    alph_to_dir = dict(zip(alphabet, directions))
    alph_to_outside = dict(zip(alphabet, outside))
    alph_to_inside = dict(zip(alphabet, inside))

    # Computing the outside points
    cur = F.zero()
    outside_points = set()
    for a in self:
        outside = cur + alph_to_outside[a]
        outside.set_immutable()
        outside_points.add(outside)
        cur += alph_to_dir[a]

    # The recursively enumerated set
    def children(p):
        possible = []
        for d in directions:
            p_d = p+d
            p_d.set_immutable()
            if p_d not in outside_points:
                possible.append(p_d)
        return possible
    seeds = [alph_to_inside[self[0]]]
    from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
    R = RecursivelyEnumeratedSet(seeds, children, structure='symmetric')
    return list(R)

def minimum_lexicographic_conjugate_reversal(self):
    r"""
    TODO: Use Lyndon factorisation to improve the time/space...

    EXAMPLES::

        sage: from slabbe.finite_word import minimum_lexicographic_conjugate_reversal
        sage: minimum_lexicographic_conjugate_reversal(Word([1,3,2,2,2]))
        word: 12223

    ::

        sage: w = words.ChristoffelWord(72452,462443)
        sage: minimum_lexicographic_conjugate_reversal(w)
        word: 0000000100000010000000100000010000001000...
        sage: _ == w
        True
        sage: minimum_lexicographic_conjugate_reversal(Word([1,3,2,2,1,1,2]))
        word: 1121322
    """
    u = minimum_lexicographic_conjugate(self)
    v = minimum_lexicographic_conjugate(self.reversal())
    return min((u,v))

def minimum_lexicographic_conjugate(self):
    r"""
    Return the conjugate word which is minimal for the lexicographic order.

    The output is a Lyndon word (or some power of).

    EXAMPLES::

        sage: from slabbe.finite_word import minimum_lexicographic_conjugate
        sage: minimum_lexicographic_conjugate(Word([1,3,2,2,2]))
        word: 13222
        sage: minimum_lexicographic_conjugate(Word([1,4,3]))
        word: 143
        sage: minimum_lexicographic_conjugate(Word([3,4,1]))
        word: 134

    The code is fast::

        sage: w = words.ChristoffelWord(72452, 462443)
        sage: minimum_lexicographic_conjugate(w)
        word: 0000000100000010000000100000010000001000...
    """
    from sage.misc.misc_c import prod
    from sage.combinat.words.word import Word
    f = self.lyndon_factorization()
    return f[-1]*prod(f[:-1], Word())

def is_lyndon_mod_reverse(self):
    r"""
    EXAMPLES::

        sage: from slabbe.finite_word import is_lyndon_mod_reverse
        sage: is_lyndon_mod_reverse(Word('111222'))
        True
        sage: is_lyndon_mod_reverse(Word('1112221'))
        False
        sage: is_lyndon_mod_reverse(Word('143'))
        False

    ::

        sage: w = words.ChristoffelWord(72452,462443)
        sage: is_lyndon_mod_reverse(w)
        True
    """
    return self.is_lyndon() and self <= minimum_lexicographic_conjugate(self.reversal())

