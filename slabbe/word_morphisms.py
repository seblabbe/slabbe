# -*- coding: utf-8 -*-
r"""
Word morphisms methods and iterators

EXAMPLES::

    sage: from slabbe.word_morphisms import iter_primitive_marked_classP_morphisms
    sage: F = FiniteWords('ab')
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
#*****************************************************************************
#       Copyright (C) 2016 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
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

def desubstitute_prefix_code(self, u):
    r"""
    Return the preimage of u under self.

    INPUT:

    - ``self`` -- word morphism, a prefix code
    - ``u`` -- word

    EXAMPLES::

        sage: from slabbe.word_morphisms import desubstitute_prefix_code
        sage: s = WordMorphism({0:[0,1],1:[1,0]})
        sage: w = desubstitute_prefix_code(s, Word([0,1,0,1,1,0]))
        sage: w
        word: 001

    The result lives in the domain of the given substitution::

        sage: w.parent()
        Finite words over {0, 1}

    TESTS::

        sage: s = WordMorphism({0:[0,1],1:[1,0],2:[1,0]})
        sage: desubstitute_prefix_code(s, Word([0,1,0,1,1,0]))
        Traceback (most recent call last):
        ...
        ValueError: non unique desubstitution, m(1)=10, m(2)=10 are prefixes of u[4:]

    ::

        sage: s = WordMorphism({0:[0,1],1:[1,0]})
        sage: desubstitute_prefix_code(s, Word([0,1,0,1,1,1]))
        Traceback (most recent call last):
        ...
        ValueError: desubstitution is impossible for u[4:]
    """
    morph = self._morph
    i = 0
    len_u = len(u)
    result = []
    while i < len_u:
        ui = u[i:]
        keys = [(k,v) for k,v in morph.iteritems() if v.is_prefix(ui)]
        if len(keys) > 1:
            s = ", ".join(["m({})={}".format(k,v) for (k,v) in keys])
            msg = ("non unique desubstitution, "
                   "{} are prefixes of u[{}:] ".format(s,i))
            raise ValueError(msg)
        if len(keys) == 0:
            raise ValueError("desubstitution is impossible for u[{}:]".format(i))
        k,v = keys[0]
        result.append(k)
        i += len(v)
    W = self.domain()
    return W(result)

def desubstitute(self, u):
    r"""
    EXAMPLES:

    Unique preimage::

        sage: from slabbe.word_morphisms import desubstitute
        sage: s = WordMorphism({0:[0,1],1:[1,0]})
        sage: desubstitute(s, Word([0,1,0,1,1,0]))
        [word: 001]

    Non-unique preimage::

        sage: s = WordMorphism({0:[0,1],1:[1,0],2:[1,0]})
        sage: desubstitute(s, Word([0,1,0,1,1,0]))
        [word: 001, word: 002]

    No preimage::

        sage: s = WordMorphism({0:[0,1],1:[1,0]})
        sage: desubstitute(s, Word([0,1,0,1,1,1]))
        []

    Lot of preimages (computation is done in parallel with Florent's Hivert
    parallel map reduce code)::

        sage: s = WordMorphism({0:[0,1],1:[0,1]})
        sage: w = Word([0,1]) ^ 10
        sage: L = desubstitute(s, w)
        sage: len(L)
        1024
    """
    morph = self._morph
    len_u = len(u)
    def successor(node):
        L,i = node
        if i == len_u:
            return []
        else:
            ui = u[i:]
            return [(L+[k],i+len(v)) for k,v in morph.iteritems() if v.is_prefix(ui)]
    roots = [([],0)]
    post_process = lambda node:node if node[1]==len_u else None
    from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
    R = RecursivelyEnumeratedSet(roots, successor, structure='forest',
            post_process=post_process)
    W = self.domain()
    map_function = lambda node:[W(node[0])]
    reduce_function = lambda x,y:x+y
    reduce_init = []
    return R.map_reduce(map_function, reduce_function, reduce_init)

def return_substitution(self, u, coding=False, length=1000):
    r"""
    Return the return substitution of self according to factor u.

    INPUT:

    - ``self`` -- word morphism
    - ``u`` -- word such that u is a prefix of self(u)
    - ``coding`` -- boolean (default: ``False``), whether to
      include the return word coding morphism
    - ``length`` -- integer (default: ``1000``), compute the first 1000 letters
      of the derived sequence to make sure every return word are seen

    EXAMPLES::

        sage: from slabbe.word_morphisms import return_substitution
        sage: s = WordMorphism({0:[0,1],1:[1,0]})
        sage: return_substitution(s, Word([0]))
        WordMorphism: 0->012, 1->02, 2->1
        sage: return_substitution(s, Word([0,1]))
        WordMorphism: 0->01, 1->23, 2->013, 3->2
        sage: return_substitution(s, Word([0,1,1]))
        WordMorphism: 0->01, 1->23, 2->013, 3->2

    ::

        sage: return_substitution(s, Word([0]), True)
        (WordMorphism: 0->012, 1->02, 2->1, 
         WordMorphism: 0->011, 1->01, 2->0)
        sage: return_substitution(s, Word([0,1]), True)
        (WordMorphism: 0->01, 1->23, 2->013, 3->2,
         WordMorphism: 0->011, 1->010, 2->0110, 3->01)

    ::

        sage: s = WordMorphism({0:[0,0,1],1:[0,1]})
        sage: return_substitution(s, Word([0]))
        WordMorphism: 0->01, 1->011

    TESTS::

        sage: s = WordMorphism({0:[0,1],1:[1,0]})
        sage: sigma_u, theta_u = return_substitution(s, Word([0]), coding=True)
        sage: sigma_u
        WordMorphism: 0->012, 1->02, 2->1
        sage: theta_u
        WordMorphism: 0->011, 1->01, 2->0
        sage: theta_u*sigma_u == s*theta_u
        True
        sage: theta_u*sigma_u
        WordMorphism: 0->011010, 1->0110, 2->01
    """
    from slabbe.infinite_word import derived_sequence
    a = u[0]
    x = self.fixed_point(a)
    s, D = derived_sequence(x, u, coding=True)
    _ = s[length] # make sure that D is complete (exact value 
                  # is known by J. Leroy and F. Durand)
    code_to_return_word = WordMorphism({v:k for k,v in D.iteritems()})
    rep = {}
    for key,value in D.iteritems():
        self_key = self(key)
        L = desubstitute(code_to_return_word, self_key)
        if len(L) == 0:
            raise ValueError("desubstitution of {} by {} "
                     "is impossible ".format(self_key, code_to_return_word))
        elif len(L) > 1:
            s = "=".join(["m({})".format(u) for u in L])
            msg = ("non unique desubstitution, "
                   "{}={}".format(s,self_key))
            raise ValueError(msg)
        #print(key,value,self_key,L[0])
        preimage = L[0]
        rep[value] = preimage
    m = WordMorphism(rep) 
    if coding:
        return m, code_to_return_word
    else:
        return m

def compute_xsi(self, u):
    r"""
    EXAMPLES::

        sage: from slabbe.word_morphisms import compute_xsi
        sage: s = WordMorphism({0:[0,1],1:[1,0]})
        sage: compute_xsi(s, Word([0]))
        sigma_u= 0->012, 1->02, 2->1
        theta_u= 0->011, 1->01, 2->0
        psi= 0->(0, 0),(0, 1),(0, 2), 1->(1, 0),(1, 1), 2->(2, 0)
        psi*sigma_u= 0->(0, 0),(0, 1),(0, 2),(1, 0),(1, 1),(2, 0), 1->(0, 0),(0, 1),(0, 2),(2, 0), 2->(1, 0),(1, 1)
        Finite words over {(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (2, 0)}
        [1 0 0]
        [1 0 0]
        [1 0 0]
        [0 1 0]
        [0 1 0]
        [0 0 1]
        We want zeta such that:
        zeta((0, 0),(0, 1),(0, 2)) = (0, 0),(0, 1),(0, 2),(1, 0),(1, 1),(2, 0)
        zeta((1, 0),(1, 1)) = (0, 0),(0, 1),(0, 2),(2, 0)
        zeta((2, 0)) = (1, 0),(1, 1)
    """
    sigma_u, theta_u = return_substitution(self, u, coding=True)
    assert theta_u*sigma_u == self*theta_u, "identity is not verified"
    print("sigma_u=", sigma_u)
    print("theta_u=", theta_u)
    d = {k:[(k,i) for i in range(len(v))] for k,v in theta_u._morph.iteritems()}
    psi = WordMorphism(d)
    print("psi=", psi)
    print("psi*sigma_u=", psi*sigma_u)
    print(psi.codomain())
    print(psi.incidence_matrix())
    print("We want zeta such that:")
    for k,v in psi._morph.iteritems():
        print("zeta({}) = {}".format(v, psi(sigma_u(k))))


