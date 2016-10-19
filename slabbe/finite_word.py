# -*- coding: utf-8 -*-
r"""
Finite words

Methods that are not in Sage (for now!)

AUTHORS:

- Sébastien Labbé, 2015

EXAMPLES::

    sage: from slabbe.finite_word import discrepancy
    sage: w = words.ChristoffelWord(5,8)
    sage: discrepancy(w)
    12/13
"""
#*****************************************************************************
#       Copyright (C) 2015 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from collections import Counter
from sage.rings.rational_field import QQ
from sage.functions.other import abs
from random import randrange

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

        sage: for c in w.conjugates(): print c, discrepancy(c)
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
