# -*- coding: utf-8 -*-
r"""
3d DyckWords

Generalisation of Dyck Word to Surface of cubes in the n x n x n cube above
the plane x+y+z=2n.

EXAMPLES::

    sage: from slabbe.dyck_3d import DyckBlocks3d
    sage: L = [len(DyckBlocks3d(i)) for i in range(1, 7)]    # not tested
    [1, 2, 9, 96, 2498, 161422]

::

    sage: L = [1, 2, 9, 96, 2498, 161422]
    sage: oeis.find_by_subsequence(L)                        # not tested internet   
    0: A115965: Number of planar subpartitions of size n pyramidal planar partition.

AUTHOR:

- Sébastien Labbé, 31 october 2014

"""
from __future__ import absolute_import, print_function

def Possible(n):
    r"""
    Possible stack of DyckWords inside a n x n cube.

    EXAMPLES::

        sage: from slabbe.dyck_3d import Possible
        sage: Possible(1)
        The Cartesian product of ({[1, 0]},)
        sage: Possible(2)
        The Cartesian product of ({[1, 1, 0, 0]}, {[1, 0, 1, 0], [1, 1, 0, 0]})
        sage: Possible(3).list()
        [([1, 1, 1, 0, 0, 0], [1, 1, 0, 1, 0, 0], [1, 0, 1, 0, 1, 0]),
         ([1, 1, 1, 0, 0, 0], [1, 1, 0, 1, 0, 0], [1, 0, 1, 1, 0, 0]),
         ([1, 1, 1, 0, 0, 0], [1, 1, 0, 1, 0, 0], [1, 1, 0, 0, 1, 0]),
         ([1, 1, 1, 0, 0, 0], [1, 1, 0, 1, 0, 0], [1, 1, 0, 1, 0, 0]),
         ([1, 1, 1, 0, 0, 0], [1, 1, 0, 1, 0, 0], [1, 1, 1, 0, 0, 0]),
         ([1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0], [1, 0, 1, 0, 1, 0]),
         ([1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0], [1, 0, 1, 1, 0, 0]),
         ([1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0], [1, 1, 0, 0, 1, 0]),
         ([1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0], [1, 1, 0, 1, 0, 0]),
         ([1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0])]

    """
    from sage.combinat.dyck_word import DyckWords, DyckWord
    from sage.categories.cartesian_product import cartesian_product
    L = []
    for i in range(1, n+1):
        K = []
        for w in DyckWords(i):
            w = DyckWord([1]*(n-i) + list(w) + [0]*(n-i))
            K.append(w) 
        L.append(K)
    return cartesian_product(L)

def is_larger_than(x, y):
    r"""
    EXAMPLES::

        sage: from slabbe.dyck_3d import is_larger_than
        sage: w = DyckWord([1,1,1,0,0,0])
        sage: w.heights()
        (0, 1, 2, 3, 2, 1, 0)
        sage: z = DyckWord([1,0,1,0,1,0])
        sage: is_larger_than(w,z)
        True
        sage: is_larger_than(z,w)
        False
        sage: is_larger_than(w,w)
        True

    """
    return all(g>=h for g,h in zip(x.heights(), y.heights()))

def DyckBlocks3d(n):
    r"""
    EXAMPLES::

        sage: from slabbe.dyck_3d import DyckBlocks3d
        sage: L = [len(DyckBlocks3d(i)) for i in range(1, 6)]
        sage: L
        [1, 2, 9, 96, 2498]
    """
    L = []
    for c in Possible(n):
        if all(is_larger_than(c[i], c[i+1]) for i in range(n-1)):
            L.append(c)
    return L

