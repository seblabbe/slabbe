# -*- coding: utf-8 -*-
r"""
Word morphisms methods and iterators

"""
#*****************************************************************************
#       Copyright (C) 2016 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.combinat.words.morphism import WordMorphism
from sage.combinat.words.words import FiniteWords
import itertools

def iter_palindromes(words, length):
    r"""
    EXAMPLES::

        sage: from slabbe.word_morphisms import iter_palindromes
        sage: list(iter_palindromes(Words('ab'), 2))
        [word: aa, word: bb]
        sage: list(iter_palindromes(Words('ab'), 3))
        [word: aaa, word: aba, word: bab, word: bbb]
        sage: list(iter_palindromes(Words('ab'), 4))
        [word: aaaa, word: abba, word: baab, word: bbbb]
        sage: list(iter_palindromes(Words('ab'), 5))
        [word: aaaaa,
         word: aabaa,
         word: ababa,
         word: abbba,
         word: baaab,
         word: babab,
         word: bbabb,
         word: bbbbb]
    """
    if length % 2 == 0:
        for w in words.iterate_by_length(length // 2):
            yield w * w.reversal()
    else:
        for w in words.iterate_by_length(length // 2):
            for a in words.iterate_by_length(1):
                yield w * a * w.reversal()

def iter_conjugate_classP(words, n):
    r"""
    EXAMPLES::

        sage: from slabbe.word_morphisms import iter_conjugate_classP
        sage: F = FiniteWords('ab')
        sage: list(iter_conjugate_classP(F, 2))
        [WordMorphism: a->a, b->a,
         WordMorphism: a->a, b->b,
         WordMorphism: a->b, b->a,
         WordMorphism: a->b, b->b]
        sage: list(iter_conjugate_classP(F, 3))
        [WordMorphism: a->aa, b->a,
         WordMorphism: a->aa, b->b,
         WordMorphism: a->bb, b->a,
         WordMorphism: a->bb, b->b,
         WordMorphism: a->a, b->aa,
         WordMorphism: a->a, b->bb,
         WordMorphism: a->b, b->aa,
         WordMorphism: a->b, b->bb,
         WordMorphism: a->ba, b->b,
         WordMorphism: a->ab, b->a,
         WordMorphism: a->b, b->ba,
         WordMorphism: a->a, b->ab]
    """
    alphabet = words.alphabet()
    length = alphabet.cardinality()
    # images are palindromes
    for sizes in IntegerListsLex(n=n, length=length, min_part=1):
        L = [iter_palindromes(words, size) for size in sizes]
        for pals in itertools.product(*L):
            d = dict(itertools.izip(alphabet, pals))
            yield WordMorphism(d, codomain=words)
    # images are one common letter + palindrome
    for sizes in IntegerListsLex(n=n-length, length=length, min_part=0):
        L = [iter_palindromes(words, size) for size in sizes]
        for pals in itertools.product(*L):
            for b in alphabet:
                d = {a:words([b])*p for a,p in itertools.izip(alphabet, pals)}
                if all(w.is_palindrome() for w in d.values()):
                    # already yielded above
                    continue
                yield WordMorphism(d, codomain=words)

def iter_pisot_irreductible(d=3, arg=None):
    r"""
    Return an iterator over Pisot irreductible substitutions

    INPUT:

    - ``d`` -- size of alphabet, [0,1,...,d-1]

    - "arg" -- (optional, default: None) It can be one of the
      following :

      * "None" -- then the method iterates through all morphisms.

      * tuple (a, b) of two integers  - It specifies the range
        "range(a, b)" of values to consider for the sum of the length

    EXAMPLES::

        sage: from slabbe.word_morphisms import iter_pisot_irreductible
        sage: it = iter_pisot_irreductible(3)
        sage: for _ in range(4): next(it)
        WordMorphism: 0->01, 1->2, 2->0
        WordMorphism: 0->02, 1->0, 2->1
        WordMorphism: 0->10, 1->2, 2->0
        WordMorphism: 0->12, 1->0, 2->1

    Pour linstant, avec le tuple, il y a un bogue::

        sage: it = iter_pisot_irreductible(3, (5,10))
        sage: for _ in range(4): next(it)
        WordMorphism: 0->0000001, 1->2, 2->0
        WordMorphism: 0->0000002, 1->0, 2->1
        WordMorphism: 0->0000010, 1->2, 2->0
        WordMorphism: 0->0000012, 1->0, 2->1
    """
    from slabbe.matrices import is_pisot
    W = FiniteWords(range(d))
    for m in W.iter_morphisms(arg):
        incidence_matrix = m.incidence_matrix()
        if not incidence_matrix.det() == 1:
            continue
        if not m.is_primitive(): # mathematiquement non necessaire
            continue
        if not is_pisot(incidence_matrix):
            continue
        yield m

def is_left_marked(m):
    r"""
    EXAMPLES::

        sage: from slabbe.word_morphisms import is_left_marked
        sage: m = WordMorphism('0->00001,1->00010')
        sage: is_left_marked(m)
        True
        sage: m = WordMorphism('0->00001,1->00001')
        sage: is_left_marked(m)
        False
        sage: m = WordMorphism('0->00001,1->00001,2->201')
        sage: is_left_marked(m)
        False
        sage: m = WordMorphism('0->00001,1->10001,2->201')
        sage: is_left_marked(m)
        True
        sage: m = WordMorphism('0->000001,1->010001,2->0201')
        sage: is_left_marked(m)
        True
        sage: m = WordMorphism('0->000001,1->010001,2->0101')
        sage: is_left_marked(m)
        False
    """
    images = m.images()
    N = len(images)
    L = map(len, images)
    for i in range(2*max(L)):
        s = len(set(image[i % L[a]] for (a,image) in enumerate(images)))
        if s == N:
            return True
        elif s > 1:
            return False
    # is cyclic
    return False

def is_marked(m):
    r"""
    EXAMPLES::

        sage: from slabbe.word_morphisms import is_marked
        sage: m = WordMorphism('0->00001,1->00010')
        sage: is_marked(m)
        True
    """
    return is_left_marked(m) and is_left_marked(m.reversal())

def iter_primitive_marked_classP_morphisms(words, n):
    r"""
    EXAMPLES::

        sage: from slabbe.word_morphisms import iter_primitive_marked_classP_morphisms
        sage: F = FiniteWords('ab')
        sage: it = iter_primitive_marked_classP_morphisms(F, 2)
        sage: list(it)
        []
        sage: it = iter_primitive_marked_classP_morphisms(F, 3)
        sage: list(it)
        [WordMorphism: a->ab, b->a, WordMorphism: a->b, b->ba]
        sage: it = iter_primitive_marked_classP_morphisms(F, 4)
        sage: list(it)
        [WordMorphism: a->aba, b->a,
         WordMorphism: a->bab, b->a,
         WordMorphism: a->b, b->aba,
         WordMorphism: a->b, b->bab,
         WordMorphism: a->abb, b->a,
         WordMorphism: a->ab, b->aa,
         WordMorphism: a->bb, b->ba,
         WordMorphism: a->b, b->baa]
    """
    for m in iter_conjugate_classP(words, n):
        if m.is_primitive() and is_marked(m):
            yield m

