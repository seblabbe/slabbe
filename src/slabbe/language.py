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
    Regular language over ['A1', 'A2', 'A3', 'P12', 'P13', 'P21', 'P23', 'P31', 'P32']
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
            Regular language over ['A1', 'A2', 'A3', 'P12', 'P13', 'P21', 'P23', 'P31', 'P32']
            defined by: Automaton with 7 states
            sage: map(L.complexity, range(4))
            [1, 9, 57, 345]
        """
        alphabet = ['A1', 'A2', 'A3', 'P12', 'P13', 'P21', 'P23', 'P31', 'P32']
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

            sage: A.process(['A1', 'P12', 'A1', 'P13'])
            (True, 'H13')
            sage: A(['A1', 'P12', 'A1', 'P13'])
            True
        """
        def H(j,k):
            return 'H%s%s' % (j,k)
        def P(j,k):
            return 'P%s%s' % (j,k)
        def A(k):
            return 'A%s' % (k)
        jk = [(3,1), (2,1), (3,2), (1,2), (2,3), (1,3)]
        D = 'Delta'
        states = [H(j,k) for (j,k) in jk] + [D]
        autom = Automaton(initial_states=[D], final_states=states)
        for (j,k) in jk:
            i = 6 - j - k
            v = H(j,k)
            autom.add_transition(v, H(i,j), P(i,j))
            autom.add_transition(v, H(k,i), P(k,i))
            autom.add_transition(v, H(j,i), P(j,i))
            autom.add_transition(v, D, A(i))
            autom.add_transition(v, v, A(j))
            autom.add_transition(D, v, P(j,k))
        for k in [1,2,3]:
            autom.add_transition(D, D, A(k))
        return autom

    def Brun(self):
        r"""
        Return the Brun regular language.

            sage: from slabbe.language import languages
            sage: L = languages.Brun()
            sage: L
            Regular language over ['B12', 'B13', 'B21', 'B23', 'B31', 'B32']
            defined by: Automaton with 6 states
            sage: map(L.complexity, range(4))
            [1, 6, 18, 54]
            sage: list(L.words_of_length_iterator(2))
            [word: B12,B12,
             word: B12,B21,
             word: B12,B31,
             word: B13,B13,
             word: B13,B21,
             word: B13,B31,
             word: B21,B12,
             word: B21,B21,
             word: B21,B32,
             word: B23,B12,
             word: B23,B23,
             word: B23,B32,
             word: B31,B13,
             word: B31,B23,
             word: B31,B31,
             word: B32,B13,
             word: B32,B23,
             word: B32,B32]
        """
        alphabet = ['B12', 'B13', 'B21', 'B23', 'B31', 'B32']
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

            sage: A(['B23', 'B23', 'B32', 'B13'])
            True
            sage: A(['B23', 'B23', 'B32', 'B13', 'B23'])
            False
        """
        def B(j,k):
            return 'B%s%s' % (j,k)
        jk = [(3,1), (2,1), (3,2), (1,2), (2,3), (1,3)]
        states = [B(j,k) for (j,k) in jk]
        autom = Automaton(initial_states=states, final_states=states)
        for (j,k) in jk:
            i = 6 - j - k
            v = B(j,k)
            autom.add_transition(v, B(j,k), B(j,k))
            autom.add_transition(v, B(k,j), B(k,j))
            autom.add_transition(v, B(i,j), B(i,j))
        return autom

languages = LanguageGenerator()
