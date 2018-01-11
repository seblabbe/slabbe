# -*- coding: utf-8 -*-
r"""
Markov transformation

EXAMPLES:

    ...

TODO:

    - Remove cylinder code from matrix cocycle
    - Remove rounded_string_vector from matrix cocycle

AUTHORS:

 - Sébastien Labbé, initial version, January 2016
"""
#*****************************************************************************
#       Copyright (C) 2016 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function
import itertools
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.matrix.constructor import matrix
from sage.combinat.finite_state_machine import Automaton
from .language import RegularLanguage
from .matrices import projection_matrix

######################
# Markov Transformation
######################
class MarkovTransformation(object):
    r"""
    Markov Transformation

    INPUT:

    - ``partition`` -- dict, mapping each key to a cone (matrix)
    - ``transitions`` -- dict, mapping each key to set of keys
    - ``linear_maps`` -- dict, mapping each key to a linear map (matrix)

    EXAMPLES:

    Brun MCF algorithm is a Markov transformation::

        sage: import itertools
        sage: B12 = matrix(3, [1,0,0, 1,1,0, 0,0,1])
        sage: B13 = matrix(3, [1,0,0, 0,1,0, 1,0,1])
        sage: B21 = matrix(3, [1,1,0, 0,1,0, 0,0,1])
        sage: B23 = matrix(3, [1,0,0, 0,1,0, 0,1,1])
        sage: B31 = matrix(3, [1,0,1, 0,1,0, 0,0,1])
        sage: B32 = matrix(3, [1,0,0, 0,1,1, 0,0,1])
        sage: gens = (B23, B32, B13, B31, B12, B21)
        sage: alphabet = [123, 132, 213, 231, 312, 321]
        sage: partition = dict(zip(alphabet, gens))
        sage: def B(i,j,k): return int('{}{}{}'.format(i,j,k))
        sage: transitions = {B(i,j,k):[B(i,j,k), B(i,k,j), B(k,i,j)]
        ....:         for i,j,k in itertools.permutations((1,2,3))}
        sage: linear_maps = partition
        sage: from slabbe.markov_transformation import MarkovTransformation 
        sage: T = MarkovTransformation(partition, transitions, linear_maps)

    """
    def __init__(self, partition, transitions, linear_maps):
        self._partition = partition
        self._transitions = transitions
        self._linear_maps = linear_maps

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: T
            Markov Transformation with transitions
            {321: [321, -321], 132: [132, -132], -123: [-132, -312], 231:
            [231, -231], -312: [-321, -231], -213: [-231, -321], 213: [213,
            -213], 312: [312, -312], -231: [-213, -123], 123: [123, -123],
            -132: [-123, -213], -321: [-312, -132]}
        """
        return ("Markov Transformation with "
               "transitions\n{}".format(self._transitions))

    def automaton(self):
        r"""
        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: T.automaton()
            Automaton with 12 states
        """
        states = self._partition.keys()
        autom = Automaton(initial_states=states, final_states=states)
        for key,values in self._transitions.iteritems():
            for v in values:
                autom.add_transition(key, v, v)
        return autom

    def language(self):
        r"""
        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: T.language()
            Regular language over [321, 132, -123, 231, -312, -213, 213,
            312, -231, 123, -132, -321]
            defined by: Automaton with 12 states
        """
        alphabet = self._partition.keys()
        return RegularLanguage(alphabet, self.automaton())

    @cached_method
    def identity_matrix(self):
        r"""
        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: T.identity_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return self._linear_maps.values()[0].parent().one()

    def word_to_matrix(self, w):
        r"""
        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: T.word_to_matrix([123,321,-231])
            [1 1 1]
            [0 1 0]
            [1 1 2]

        Empty word::

            sage: T.word_to_matrix([])
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return prod((self._linear_maps[a] for a in w), z=self.identity_matrix())

    def n_words_iterator(self, n):
        r"""
        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: list( T.n_words_iterator(1))
            [word: 321, word: 132, word: -123, word: 231, word: -312, word:
            -213, word: 213, word: 312, word: -231, word: 123, word: -132,
            word: -321]

        TESTS::

            sage: list(T.n_words_iterator(0))
            [word: ]
        """
        return self.language().words_of_length_iterator(n)

    def n_matrices_iterator(self, n):
        r"""
        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: A,B = zip(*list(T.n_matrices_iterator(1)))
            sage: A
            (word: 321, word: 132, word: -123, word: 231, word: -312, word:
            -213, word: 213, word: 312, word: -231, word: 123, word: -132,
            word: -321)
            sage: B
            (
            [1 0 1]  [1 0 0]  [1 0 0]  [1 1 0]  [1 0 0]  [1 0 0]  [1 0 0]  [1 0 0]
            [0 1 0]  [1 1 0]  [0 1 0]  [0 1 0]  [0 1 1]  [0 1 0]  [0 1 0]  [0 1 1]
            [0 0 1], [0 0 1], [1 0 1], [0 0 1], [0 0 1], [0 1 1], [0 1 1], [0 0 1],
            <BLANKLINE>
            [1 1 0]  [1 0 0]  [1 0 0]  [1 0 1]
            [0 1 0]  [0 1 0]  [1 1 0]  [0 1 0]
            [0 0 1], [1 0 1], [0 0 1], [0 0 1]
            )

        TESTS::

            sage: list(T.n_matrices_iterator(0))
            [(
                    [1 0 0]
                    [0 1 0]
            word: , [0 0 1]
            )]
        """
        for w in self.n_words_iterator(n):
            yield w, self.word_to_matrix(w)

    def n_cylinders_iterator(self, n):
        r"""
        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: A,B = zip(*list(T.n_cylinders_iterator(1)))
            sage: A[:5]
            (word: 321, word: 321, word: 132, word: 132, word: -123)
            sage: B
            (
            [1 3 1]  [2 3 1]  [0 1 0]  [1 1 0]  [1 1 0]  [1 1 1]  [1 3 1]  [2 3 1]
            [0 1 1]  [1 1 1]  [1 3 1]  [2 3 1]  [1 2 1]  [1 2 1]  [0 1 0]  [1 1 0]
            [0 1 0], [1 1 0], [0 1 1], [1 1 1], [2 2 1], [2 2 1], [0 1 1], [1 1 1],
            <BLANKLINE>
            [1 2 1]  [1 2 1]  [1 2 1]  [1 2 1]  [0 1 1]  [1 1 1]  [0 1 1]  [1 1 1]
            [2 2 1]  [2 2 1]  [1 1 0]  [1 1 1]  [0 1 0]  [1 1 0]  [1 3 1]  [2 3 1]
            [1 1 0], [1 1 1], [2 2 1], [2 2 1], [1 3 1], [2 3 1], [0 1 0], [1 1 0],
            <BLANKLINE>
            [2 2 1]  [2 2 1]  [0 1 0]  [1 1 0]  [1 1 0]  [1 1 1]  [2 2 1]  [2 2 1]
            [1 1 0]  [1 1 1]  [0 1 1]  [1 1 1]  [2 2 1]  [2 2 1]  [1 2 1]  [1 2 1]
            [1 2 1], [1 2 1], [1 3 1], [2 3 1], [1 2 1], [1 2 1], [1 1 0], [1 1 1]
            )
        """
        for w,m in self.n_matrices_iterator(n):
            if w:
                parts = self._transitions[w[-1]]
            else:
                parts = self._partition.keys()
            for part in parts:
                part_matrix = self._partition[part]
                yield w, m*part_matrix

    def n_cylinders_edges(self, n):
        r"""
        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: E = T.n_cylinders_edges(1)
            sage: len(E)
            39
        """
        from sage.rings.finite_rings.integer_mod_ring import Integers
        edges = set()
        for w,cyl in self.n_cylinders_iterator(n):
            cols = cyl.columns()
            indices = Integers(len(cols))
            edges.update(frozenset((cols[i], cols[i+1])) for i in indices)
        return edges

    def plot_n_cylinders(self, n, labels=True):
        r"""
        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: G = T.plot_n_cylinders(3)

        TESTS::

            sage: G = T.plot_n_cylinders(0)
        """
        from sage.plot.graphics import Graphics
        from sage.plot.polygon import polygon
        from sage.plot.text import text
        M3to2 = projection_matrix(3, 2)
        G = Graphics()
        for w,cyl in self.n_cylinders_iterator(n):
            columns = cyl.columns()
            G += polygon((M3to2*col/col.norm(1) for col in columns), fill=False) 
            if labels:
                sum_cols = sum(columns)
                G += text("{}".format(w), M3to2*sum_cols/sum_cols.norm(1))
        return G

    def tikz_n_cylinders(self, n, labels=None, scale=1):
        r"""
        INPUT:

        - ``labels`` -- None, True or False (default: None), if None, it
          takes value True if n is 1.

        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: t = T.tikz_n_cylinders(1, labels=True, scale=4)
            sage: t
            \documentclass[tikz]{standalone}
            \usepackage{amsmath}
            \begin{document}
            \begin{tikzpicture}
            [scale=4]
            \draw (-0.8660, -0.5000) -- (-0.3464, -0.2000);
            \draw (0.0000, 0.0000) -- (-0.2165, -0.1250);
            \draw (0.0000, -0.2000) -- (0.0000, 0.0000);
            ...
            ... 56 lines not printed (2702 characters in total) ...
            ...
            \node at (0.2742, 0.0750) {$-132$};
            \node at (0.1299, -0.0083) {$-132$};
            \node at (-0.0722, -0.2750) {$-321$};
            \node at (-0.0722, -0.1083) {$-321$};
            \end{tikzpicture}
            \end{document}

        ::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.pdf')
            sage: _ = t.pdf(filename)
        """
        if labels is None:
            labels = True if n == 1 else False
        lines = []
        lines.append(r"\begin{tikzpicture}")
        lines.append("[scale={}]".format(scale))
        M3to2 = projection_matrix(3, 2)
        for (u,v) in self.n_cylinders_edges(n):
            u = rounded_string_vector(M3to2 * u / u.norm(1), digits=4)
            v = rounded_string_vector(M3to2 * v / v.norm(1), digits=4)
            lines.append(r"\draw {} -- {};".format(u,v))
        if labels:
            for w,cyl in self.n_cylinders_iterator(n):
                u = sum(c / c.norm(1) for c in cyl.columns())
                u = rounded_string_vector(M3to2 * u / u.norm(1), digits=4)
                lines.append(r"\node at {} {{${}$}};".format(u, w))
        lines.append(r"\end{tikzpicture}")
        from slabbe import TikzPicture
        return TikzPicture("\n".join(lines))

##################################
# Markov Transformation Generators
##################################

class MarkovTransformationGenerators(object):
    def Brun(self):
        B12 = matrix(3, [1,0,0, 1,1,0, 0,0,1])
        B13 = matrix(3, [1,0,0, 0,1,0, 1,0,1])
        B21 = matrix(3, [1,1,0, 0,1,0, 0,0,1])
        B23 = matrix(3, [1,0,0, 0,1,0, 0,1,1])
        B31 = matrix(3, [1,0,1, 0,1,0, 0,0,1])
        B32 = matrix(3, [1,0,0, 0,1,1, 0,0,1])
        gens = (B23, B32, B13, B31, B12, B21)
        alphabet = [123, 132, 213, 231, 312, 321]
        linear_maps = dict(zip(alphabet, gens))

        partition = linear_maps 

        def B(i,j,k): return int('{}{}{}'.format(i,j,k))
        transitions = {B(i,j,k):[B(i,j,k), B(i,k,j), B(k,i,j)]
                for i,j,k in itertools.permutations((1,2,3))}

        return MarkovTransformation(partition, transitions, linear_maps)

    def Selmer(self):
        r"""
        EXAMPLES::

            sage: from slabbe.markov_transformation import markov_transformations
            sage: T = markov_transformations.Selmer()
            sage: T
            Markov Transformation with transitions
            {321: [321, -321], 132: [132, -132], -123: [-132, -312], 231:
            [231, -231], -312: [-321, -231], -213: [-231, -321], 213: [213,
            -213], 312: [312, -312], -231: [-213, -123], 123: [123, -123],
            -132: [-123, -213], -321: [-312, -132]}
        """
        L = {}
        L[123] = L[-123] = matrix(3, [1,0,0, 0,1,0, 1,0,1])
        L[132] = L[-132] = matrix(3, [1,0,0, 1,1,0, 0,0,1])
        L[213] = L[-213] = matrix(3, [1,0,0, 0,1,0, 0,1,1])
        L[231] = L[-231] = matrix(3, [1,1,0, 0,1,0, 0,0,1])
        L[312] = L[-312] = matrix(3, [1,0,0, 0,1,1, 0,0,1])
        L[321] = L[-321] = matrix(3, [1,0,1, 0,1,0, 0,0,1])
        linear_maps = L

        partition = {}
        partition[213] = matrix.column(3, [0,0,1, 1,1,2, 1,0,1])
        partition[123] = matrix.column(3, [0,0,1, 1,1,2, 0,1,1])
        partition[132] = matrix.column(3, [0,1,0, 1,2,1, 0,1,1])
        partition[312] = matrix.column(3, [0,1,0, 1,2,1, 1,1,0])
        partition[321] = matrix.column(3, [1,0,0, 2,1,1, 1,1,0])
        partition[231] = matrix.column(3, [1,0,0, 2,1,1, 1,0,1])
        partition[-213] = matrix.column(3, [1,1,1, 1,1,2, 1,0,1])
        partition[-123] = matrix.column(3, [1,1,1, 1,1,2, 0,1,1])
        partition[-132] = matrix.column(3, [1,1,1, 1,2,1, 0,1,1])
        partition[-312] = matrix.column(3, [1,1,1, 1,2,1, 1,1,0])
        partition[-321] = matrix.column(3, [1,1,1, 2,1,1, 1,1,0])
        partition[-231] = matrix.column(3, [1,1,1, 2,1,1, 1,0,1])

        def B(i,j,k): return int('{}{}{}'.format(i,j,k))
        transitions = {}
        for i,j,k in itertools.permutations((1,2,3)):
            transitions[B(i,j,k)] = [B(i,j,k), -B(i,j,k)]
            transitions[-B(i,j,k)] = [-B(i,k,j), -B(k,i,j)]

        return MarkovTransformation(partition, transitions, linear_maps)

markov_transformations = MarkovTransformationGenerators()

#####################
# Helper functions
#####################
def rounded_string_vector(v, digits=4):
    r"""
    EXAMPLES::

        sage: from slabbe.matrix_cocycle import rounded_string_vector
        sage: v = (-0.144337567297406, 0.166666666666667)
        sage: rounded_string_vector(v)
        '(-0.1443, 0.1667)'
        sage: rounded_string_vector(v, digits=6)
        '(-0.144338, 0.166667)'
    """
    s = "{{:.{}f}}".format(digits)
    content = ", ".join(s.format(float(a)) for a in v)
    return "({})".format(content)
