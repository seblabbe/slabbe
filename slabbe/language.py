# -*- coding: utf-8 -*-
r"""
Regular languages

EXAMPLES:

Language over all finite words on an alphabet::

    sage: from slabbe.language import Language
    sage: Language(alphabet=['a', 'b'])
    Language of finite words over alphabet ['a', 'b']

Finite language::

    sage: from slabbe.language import FiniteLanguage
    sage: S = ['a', 'ab', 'aab', 'aaab']
    sage: FiniteLanguage(alphabet=['a', 'b'], words=S)
    Finite language of cardinality 4 over alphabet ['a', 'b']

Regular language::

    sage: from slabbe.language import RegularLanguage
    sage: alphabet = ['a', 'b']
    sage: trans = [(0, 1, 'a'), (1, 2, 'b'), (2, 3, 'b'), (3, 4, 'a')]
    sage: automaton = Automaton(trans, initial_states=[0], final_states=[4])
    sage: RegularLanguage(alphabet, automaton)
    Regular language over ['a', 'b']
    defined by: Automaton with 5 states

Predefined languages::

    sage: from slabbe.language import languages
    sage: languages.ARP()
    Regular language over [1, 2, 3, 123, 132, 213, 231, 312, 321]
    defined by: Automaton with 7 states

AUTHORS:

 - Sébastien Labbé, initial clean and full doctested version, October 2015
"""
#*****************************************************************************
#       Copyright (C) 2015 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
import itertools
from sage.combinat.words.words import Words
from sage.combinat.finite_state_machine import Automaton

class Language(object):
    r"""
    Language of finite words

    INPUT:

    - ``alphabet`` -- iterable of letters

    EXAMPLES::

        sage: from slabbe.language import Language
        sage: Language(alphabet=['a', 'b'])
        Language of finite words over alphabet ['a', 'b']
    """
    def __init__(self, alphabet):
        r"""
        EXAMPLES::

            sage: from slabbe.language import Language
            sage: Language(alphabet=['a', 'b'])
            Language of finite words over alphabet ['a', 'b']
        """
        self._alphabet = alphabet

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe.language import Language
            sage: Language(alphabet=['a'])
            Language of finite words over alphabet ['a']
        """
        s = "Language of finite words over alphabet {}"
        return s.format(self._alphabet)

    def words_of_length_iterator(self, length):
        r"""
        Return an iterator over words of given length.

        INPUT:

        - ``length`` -- integer

        EXAMPLES::

            sage: from slabbe.language import Language
            sage: F = Language(alphabet=['a', 'b'])
            sage: it = F.words_of_length_iterator(2)
            sage: list(it)
            [word: aa, word: ab, word: ba, word: bb]
        """
        W = Words(self._alphabet)
        return W.iterate_by_length(length)

    def complexity(self, length):
        r"""
        Returns the number of words of given length.

        .. NOTE::

            This method is defined from :func:`~words_of_length_iterator`.

        INPUT:

        - ``length`` -- integer

        EXAMPLES::

            sage: from slabbe.language import Language
            sage: F = Language(alphabet=['a', 'b'])
            sage: map(F.complexity, range(5))
            [1, 2, 4, 8, 16]
        """
        return sum(1 for _ in self.words_of_length_iterator(length))
class FiniteLanguage(Language):
    r"""
    Finite language

    INPUT:

    - ``alphabet`` -- iterable of letters
    - ``words`` -- finite iterable of words

    EXAMPLES::

        sage: from slabbe.language import FiniteLanguage
        sage: L = ['a', 'aa', 'aaa']
        sage: FiniteLanguage(alphabet=['a'], words=L)
        Finite language of cardinality 3 over alphabet ['a']
    """
    def __init__(self, alphabet, words):
        r"""
        EXAMPLES::

            sage: from slabbe.language import FiniteLanguage
            sage: L = ['a', 'aa', 'aaa']
            sage: F = FiniteLanguage(alphabet=['a'], words=L)
        """
        Language.__init__(self, alphabet)
        self._words = words
    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe.language import FiniteLanguage
            sage: L = ['a', 'ab', 'aab', 'aaab']
            sage: FiniteLanguage(alphabet=['a', 'b'], words=L)
            Finite language of cardinality 4 over alphabet ['a', 'b']
        """
        s = "Finite language of cardinality {} over alphabet {}"
        return s.format(len(self._words), self._alphabet)

    def automaton(self):
        r"""
        Return the automaton recognizing this finite language.

        EXAMPLES::

            sage: from slabbe.language import FiniteLanguage
            sage: L = ['a', 'aa', 'aaa']
            sage: F = FiniteLanguage(alphabet=['a'], words=L)
            sage: F.automaton()
            Automaton with 7 states
        """
        transitions = []
        final_states = []
        end = 1
        for w in self._words:
            start = 0
            for a in w:
                transitions.append((start, end, a))
                start, end = end, end+1
            final_states.append(start)
        return Automaton(transitions, initial_states=[0], final_states=final_states)

    def minimal_automaton(self):
        r"""
        Return the minimal automaton recognizing this finite language.

        .. NOTE:: 
        
            One of the state is not final. You may want to remove it...

        EXAMPLES::

            sage: from slabbe.language import FiniteLanguage
            sage: L = ['a', 'aa', 'aaa']
            sage: F = FiniteLanguage(alphabet=['a'], words=L)
            sage: F.minimal_automaton()
            Automaton with 5 states
        """
        return self.automaton().minimization().relabeled()

    def number_of_states(self):
        r"""
        EXAMPLES::

            sage: from slabbe.language import FiniteLanguage
            sage: L = ['a', 'aa', 'aaa']
            sage: F = FiniteLanguage(alphabet=['a'], words=L)
            sage: F.number_of_states()
            5
        """
        return len(self.minimal_automaton().states())

class RegularLanguage(Language):
    r"""
    Regular language

    INPUT:

    - ``alphabet`` -- iterable of letters
    - ``automaton`` -- finite state automaton

    EXAMPLES::

        sage: from slabbe.language import RegularLanguage
        sage: alphabet = ['a', 'b']
        sage: trans = [(0, 1, 'a'), (1, 2, 'b'), (2, 3, 'b'), (3, 4, 'a')]
        sage: automaton = Automaton(trans, initial_states=[0], final_states=[4])
        sage: RegularLanguage(alphabet, automaton)
        Regular language over ['a', 'b']
        defined by: Automaton with 5 states
    """
    def __init__(self, alphabet, automaton):
        r"""
        EXAMPLES::

            sage: from slabbe.language import RegularLanguage
            sage: alphabet = ['a', 'b']
            sage: trans = [(0, 1, 'a'), (1, 2, 'b'), (2, 3, 'b'), (3, 4, 'a')]
            sage: automaton = Automaton(trans, initial_states=[0], final_states=[4])
            sage: R = RegularLanguage(alphabet, automaton)
        """
        Language.__init__(self, alphabet)
        self._automaton = automaton


    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe.language import RegularLanguage
            sage: alphabet = ['a', 'b']
            sage: trans = [(0, 1, 'a'), (1, 2, 'b'), (2, 3, 'b'), (3, 4, 'a')]
            sage: automaton = Automaton(trans, initial_states=[0], final_states=[4])
            sage: RegularLanguage(alphabet, automaton)
            Regular language over ['a', 'b']
            defined by: Automaton with 5 states
        """
        s = "Regular language over {}\ndefined by: {}"
        return s.format(self._alphabet, self._automaton)

    def words_of_length_iterator(self, length):
        r"""
        Return an iterator over words of given length.

        INPUT:

        - ``length`` -- integer

        EXAMPLES::

            sage: from slabbe.language import RegularLanguage
            sage: alphabet = ['a', 'b']
            sage: trans = [(0, 1, 'a'), (1, 2, 'b'), (2, 3, 'b'), (3, 4, 'a')]
            sage: automaton = Automaton(trans, initial_states=[0], final_states=[4])
            sage: R = RegularLanguage(alphabet, automaton)
            sage: [list(R.words_of_length_iterator(i)) for i in range(6)]
            [[], [], [], [], [word: abba], []]
        """
        it = super(RegularLanguage, self).words_of_length_iterator(length)
        return itertools.ifilter(self._automaton, it)

#####################
# Language generators
#####################
class LanguageGenerator(object):
    def ARP(self):
        r"""
        Return the Arnoux-Rauzy-Poincaré regular language.

            sage: from slabbe.language import languages
            sage: L = languages.ARP()
            sage: L
            Regular language over [1, 2, 3, 123, 132, 213, 231, 312, 321]
            defined by: Automaton with 7 states
            sage: map(L.complexity, range(4))
            [1, 9, 57, 345]
        """
        alphabet = [1, 2, 3, 123, 132, 213, 231, 312, 321]
        automaton = self._ARP_automaton()
        return RegularLanguage(alphabet, automaton)

    def _ARP_automaton(self):
        r"""
        Return the automaton for the ARP language.

        EXAMPLES::

            sage: from slabbe.language import languages
            sage: A = languages._ARP_automaton()
            sage: A
            Automaton with 7 states

        TESTS::

            sage: A.process([1, 312, 1, 213])
            (True, 'H213')
            sage: A([1, 312, 1, 213])
            True
        """
        def H(i,j,k):
            return 'H{}{}{}'.format(i,j,k)
        def P(i,j,k):
            return int('{}{}{}'.format(i,j,k))
        def A(k):
            return int('{}'.format(k))
        D = 'Delta'
        states = [H(*p) for p in itertools.permutations((1,2,3))] + [D]
        autom = Automaton(initial_states=[D], final_states=states)
        for p in itertools.permutations((1,2,3)):
            i,j,k = p
            v = H(*p)
            autom.add_transition(v, H(k,i,j), P(k,i,j))
            autom.add_transition(v, H(j,k,i), P(j,k,i))
            autom.add_transition(v, H(k,j,i), P(k,j,i))
            autom.add_transition(v, D, A(i))
            autom.add_transition(v, v, A(j))
            autom.add_transition(D, v, P(i,j,k))
        for k in [1,2,3]:
            autom.add_transition(D, D, A(k))
        return autom

    def Brun(self):
        r"""
        Return the Brun regular language.

        EXAMPLES::

            sage: from slabbe.language import languages
            sage: L = languages.Brun()
            sage: L
            Regular language over [123, 132, 213, 231, 312, 321]
            defined by: Automaton with 6 states
            sage: map(L.complexity, range(4))
            [1, 6, 18, 54]
            sage: list(L.words_of_length_iterator(2))
            [word: 123,123,
             word: 123,132,
             word: 123,312,
             word: 132,123,
             word: 132,132,
             word: 132,213,
             word: 213,213,
             word: 213,231,
             word: 213,321,
             word: 231,123,
             word: 231,213,
             word: 231,231,
             word: 312,231,
             word: 312,312,
             word: 312,321,
             word: 321,132,
             word: 321,312,
             word: 321,321]
        """
        alphabet = [123, 132, 213, 231, 312, 321]
        automaton = self._Brun_automaton()
        return RegularLanguage(alphabet, automaton)

    def _Brun_automaton(self):
        r"""
        Return the automaton for the Brun language.

        EXAMPLES::

            sage: from slabbe.language import languages
            sage: A = languages._Brun_automaton()
            sage: A
            Automaton with 6 states

        TESTS::

            sage: A([123, 123, 132, 213])
            True
            sage: A([123, 123, 132, 213, 123])
            False
        """
        def B(i,j,k):
            return int('{}{}{}'.format(i,j,k))
        states = [B(*p) for p in itertools.permutations((1,2,3))]
        autom = Automaton(initial_states=states, final_states=states)
        for p in itertools.permutations((1,2,3)):
            i,j,k = p
            autom.add_transition(B(*p), B(i,j,k), B(i,j,k))
            autom.add_transition(B(*p), B(i,k,j), B(i,k,j))
            autom.add_transition(B(*p), B(k,i,j), B(k,i,j))
        return autom

    def Selmer(self):
        r"""
        Return the Selmer regular language.

        EXAMPLES::

            sage: from slabbe.language import languages
            sage: L = languages.Selmer()
            sage: L
            Regular language over [123, 132, 213, 231, 312, 321]
            defined by: Automaton with 6 states
            sage: map(L.complexity, range(4))
            [1, 6, 12, 24]
            sage: list(L.words_of_length_iterator(2))
            [word: 123,132,
             word: 123,312,
             word: 132,123,
             word: 132,213,
             word: 213,231,
             word: 213,321,
             word: 231,123,
             word: 231,213,
             word: 312,231,
             word: 312,321,
             word: 321,132,
             word: 321,312]
        """
        alphabet = [123, 132, 213, 231, 312, 321]
        automaton = self._Selmer_automaton()
        return RegularLanguage(alphabet, automaton)

    def _Selmer_automaton(self):
        r"""
        Return the automaton for the Selmer language.

        EXAMPLES::

            sage: from slabbe.language import languages
            sage: A = languages._Selmer_automaton()
            sage: A
            Automaton with 6 states

        TESTS::

            sage: A([123, 132, 213])
            True
            sage: A([123, 132, 213, 123])
            False
        """
        def S(i,j,k):
            return int('{}{}{}'.format(i,j,k))
        states = [S(*p) for p in itertools.permutations((1,2,3))]
        autom = Automaton(initial_states=states, final_states=states)
        for p in itertools.permutations((1,2,3)):
            i,j,k = p
            autom.add_transition(S(*p), S(i,k,j), S(i,k,j))
            autom.add_transition(S(*p), S(k,i,j), S(k,i,j))
        return autom

    def Cassaigne(self):
        r"""
        Return the Cassaigne regular language over the alphabet
        [11, 22, 122, 211, 121, 212].

        EXAMPLES::

            sage: from slabbe.language import languages
            sage: L = languages.Cassaigne()
            sage: L
            Regular language over [11, 22, 122, 211, 121, 212]
            defined by: Automaton with 1 state
            sage: map(L.complexity, range(4))
            [1, 6, 36, 216]
        """
        alphabet = [11, 22, 122, 211, 121, 212]
        automaton = self._Cassaigne_automaton()
        return RegularLanguage(alphabet, automaton)

    def _Cassaigne_automaton(self):
        r"""
        Return the automaton for the Cassaigne language over the alphabet
        [11, 22, 122, 211, 121, 212].

        EXAMPLES::

            sage: from slabbe.language import languages
            sage: A = languages._Cassaigne_automaton()
            sage: A
            Automaton with 1 state
        """
        q = 0
        states = [q]
        autom = Automaton(initial_states=states, final_states=states)
        alphabet = [11, 22, 122, 211, 121, 212]
        for i in alphabet:
            autom.add_transition(q, q, i)
        return autom

languages = LanguageGenerator()
