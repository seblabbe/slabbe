# -*- coding: utf-8 -*-
r"""
Bispecial factors and extension types

EXAMPLES::

    sage: print "todo"
    not done

TODO:

    - add examples above
    - use __classcall_private__ stuff for ExtensionType ?
    - rename ExtensionType2to1 to ExtendedExtensionType ?
    - export tikz to pdf using view instead of tikz2pdf ?
    - fix bug of apply for ExtensionType2to1 when the word appears in the
      image of a letter
    - add the bispecial word to the attribute of extension type
    - use this to compute the factor complexity function

"""
from collections import defaultdict, Counter
import itertools
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.cachefunc import cached_method
from sage.combinat.words.word import FiniteWord_class, Word
from sage.combinat.words.morphism import WordMorphism

common_substitutions_dict = dict(
ar1=WordMorphism({1:[1],      2:[2,1],   3:[3,1]}),
ar2=WordMorphism({1:[1,2],   2:[2],     3:[3,2]}),
ar3=WordMorphism({1:[1,3],   2:[2,3],   3:[3]}),
b12=WordMorphism({1:[1,2],   2:[2],     3:[3]}),
b13=WordMorphism({1:[1,3],   2:[2],     3:[3]}),
b21=WordMorphism({1:[1],     2:[2,1],   3:[3]}),
b23=WordMorphism({1:[1],     2:[2,3],   3:[3]}),
b31=WordMorphism({1:[1],     2:[2],     3:[3,1]}),
b32=WordMorphism({1:[1],     2:[2],     3:[3,2]}),
p12=WordMorphism({1:[1,2],   2:[2],     3:[3,1,2]}),
p13=WordMorphism({1:[1,3],   2:[2,1,3], 3:[3]}),
p21=WordMorphism({1:[1],     2:[2,1],   3:[3,2,1]}),
p23=WordMorphism({1:[1,2,3], 2:[2,3],   3:[3]}),
p31=WordMorphism({1:[1],     2:[2,3,1], 3:[3,1]}),
p32=WordMorphism({1:[1,3,2], 2:[2],     3:[3,2]})
)

######################################
# Utility functions
######################################
def factors_length_2_from_morphism_and_factors_length_2(m, F):
    r"""
    Return the set of factors of lengths two in the image by a morphism of
    a set of factors of length 2.

    INPUT:

    - ``m`` - endomorphim
    - ``F`` - set of factors of length 2

    EXAMPLES::

        sage: b12 = WordMorphism({1:[1,2],2:[2],3:[3]})
        sage: sorted(factors_length_2_from_morphism_and_factors_length_2(b12, []))
        []
        sage: sorted(factors_length_2_from_morphism_and_factors_length_2(b12, [(1,1)]))
        [(1, 2), (2, 1)]
        sage: b23 = WordMorphism({1:[1],2:[2,3],3:[3]})
        sage: sorted(factors_length_2_from_morphism_and_factors_length_2(b23, [(1,1)]))
        [(1, 1)]

    """
    assert m.is_endomorphism(), "m(=%s) must be an endomorphism"
    L = []
    # externes
    letters = set()
    for a,b in F:
        L.append( (m(a)[-1], m(b)[0]) )
        letters.add(a)
        letters.add(b)
    # internes
    for a in letters:
        image = m(a)
        for i in range(len(image)-1):
            L.append( (image[i], image[i+1]) )
    return set(L)


######################################
# Extension Type
######################################
class ExtensionType(object):
    #__metaclass__ = ClasscallMetaclass
    #@staticmethod
    #def __classcall_private__(cls, *args, **kwds):
    #    if len(args) == 2:
    #        a0, a1 = args
    #        if (isinstance(a0, FiniteWord_class) 
    #            and isinstance(a1, FiniteWord_class)):
    #            return self.from_factor(a0, a1)
    #        else:
    #            raise NotImplementedError
    #    else:
    #        raise NotImplementedError

    @staticmethod
    def from_factor(bispecial, word):
        r"""
        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: ExtensionType.from_factor(W(), W([0,1,1,2,0]))
              E(w)   0   1   2
               0         X
               1         X   X
               2     X
             m(w)=-1, not ord.
        """
        W = word.parent()
        letters = list(W.iterate_by_length(1))
        L = []
        for a,b in itertools.product(letters, letters):
            if (a*bispecial*b).is_factor(word):
                L.append((a[0],b[0]))
        return ExtensionType1to1(L, W.alphabet())

    @staticmethod
    def from_factor2(bispecial, word):
        r"""
        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: ExtensionType.from_factor2(W(), W([0,1,1,2,0]))
              E(w)   0   1   2
               0         X
               1         X   X
               2     X
             m(w)=-1, not ord.
        """
        W = word.parent()
        len_bispecial = len(bispecial)
        L = []
        for f in word.factor_iterator(bispecial+2):
            if bispecial == f[1:-1]:
                L.append((f[0], f[-1]))
        return ExtensionType1to1(L, W.alphabet())

    @staticmethod
    def from_morphism(m):
        r"""
        Return the extension type of the empty word in the language defined by
        the image of the free monoid under the morphism m.

        INPUT:

        - ``m`` - endomorphim

        EXAMPLES::

            sage: ar = WordMorphism({1:[1,3],2:[2,3],3:[3]})
            sage: ExtensionType.from_morphism(ar)
              E(w)   1   2   3
               1             X
               2             X
               3     X   X   X
             m(w)=0, ordinary

        ::

            sage: p = WordMorphism({1:[1,2,3],2:[2,3],3:[3]})
            sage: ExtensionType.from_morphism(p)
              E(w)   1   2   3
               1         X    
               2             X
               3     X   X   X
             m(w)=0, not ord.

        ::

            sage: b12 = WordMorphism({1:[1,2],2:[2],3:[3]})
            sage: ExtensionType.from_morphism(b12)
              E(w)   1   2   3
               1         X    
               2     X   X   X
               3     X   X   X
             m(w)=2, not ord.

        """
        assert m.is_endomorphism(), "m(=%s) must be an endomorphism"
        L = []
        images = m.images()
        # externes
        for ma in images:
            for mb in images:
                L.append( (ma[-1],mb[0]) )
        # internes
        for image in images:
            for i in range(len(image)-1):
                L.append( (image[i], image[i+1]) )
        alphabet = m.domain().alphabet()
        return ExtensionType1to1(L, alphabet)

    def __hash__(self):
        return hash(self._pairs)
    def __iter__(self):
        return iter(self._pairs)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E
              E(w)   1   2   3
               1             X
               2             X
               3     X   X   X
              m(w)=0, ordinary

        With chignons::

            sage: E = ExtensionType1to1(L, [1,2,3], ('a','b'))
            sage: E
              E(awb)   1   2   3
                1              X
                2              X
                3      X   X   X
              m(w)=0, ordinary

        ::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:     2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: ExtensionType2to1(L, (1,2,3))
              E(w)   1   2   3
               21        X
               31        X
               12    X   X   X
               22    X
               23    X
              m(w)=0, not ord., empty
        """
        ordinary = 'ordinary' if self.is_ordinaire() else "not ord."
        empty = ', empty' if self.is_empty() else ""
        s = "\n m(w)={}, {}{}".format(self.multiplicity(), ordinary, empty)
        return self.table()._repr_() + s

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3], ('a', 'b'))
            sage: latex(E)
            \begin{tabular}{c}
            \begin{tabular}{cccc}
            E(awb) & $1$ & $2$ & $3$ \\
            $1$ &   &   & X \\
            $2$ &   &   & X \\
            $3$ & X & X & X \\
            \end{tabular}\\
            $m(w) = 0$, ordinary
            \end{tabular}

        ::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:     2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: latex(E)
            \begin{tabular}{c}
            \begin{tabular}{cccc}
            E(w) & $1$ & $2$ & $3$ \\
            $21$ &   & X &   \\
            $31$ &   & X &   \\
            $12$ & X & X & X \\
            $22$ & X &   &   \\
            $23$ & X &   &   \\
            \end{tabular}\\
            $m(w) = 0$, not ord., empty
            \end{tabular}

        """
        ordinary = 'ordinary' if self.is_ordinaire() else "not ord."
        empty = ', empty' if self.is_empty() else ""
        s = '\\begin{tabular}{c}\n'
        s += self.table()._latex_()
        s += "\\\\\n$m(w) = {}$, {}{}".format(self.multiplicity(), ordinary, empty)
        s += '\n\\end{tabular}'
        return s

    def life_graph(self, substitutions, substitutions_dict=None):
        r"""
        Return the graph of extension types generated under a sequence of
        substitutions.

        INPUT:

        - ``substitutions`` - list of substitutions keys
        - ``substitutions_dict`` - dict of substitutions, if None then it
          gets replaced by ``common_substitutions_dict`` defined in the
          module.

        EXAMPLES:

        From an ordinaire word::

            sage: e = ExtensionType1to1([(1,3),(2,3),(3,1),(3,2),(3,3)], [1,2,3])
            sage: e.life_graph(['p23'])
            Looped multi-digraph on 2 vertices
            sage: e.life_graph(['p32'])
            Looped multi-digraph on 3 vertices
            sage: e.life_graph(['p32','p13','ar2'])
            Looped multi-digraph on 5 vertices

        2to1::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:     2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.life_graph(['b12','b21','b12'])
            Looped multi-digraph on 8 vertices
        """
        if substitutions_dict is None:
            substitutions_dict = common_substitutions_dict
        multiedges = True
        loops = True
        G = DiGraph(multiedges=multiedges,loops=loops)
        L,newL = [self],[]
        for key in substitutions:
            s = substitutions_dict[key]
            for e in L:
                for new in e.apply(s):
                    G.add_edge((e,new,key))
                    newL.append(new)
            L,newL = newL,[]
        return G

    def life_graph_save_tikz(self, filename, substitutions, **kwds):
        r"""
        INPUT:

        - "filename" - string
        - "format" - string, default: 'tkz_graph' -- either 'dot2tex' or
          'tkz_graph'.

        If format is 'dot2tex', then all the LaTeX generation will be
        delegated to "dot2tex" (which must be installed).

        For the 'dot2tex' format, the possible option names and associated
        values are given below:

        - "prog" -- the program used for the layout. It must be a string
          corresponding to one of the software of the graphviz suite:
          'dot', 'neato', 'twopi', 'circo' or 'fdp'.

        See this for more options::

            sage: g = graphs.PetersenGraph()
            sage: opts = g.latex_options()
            sage: opts.set_option?        # not tested

        EXAMPLES::

            sage: e = ExtensionType1to1([(1,3),(2,3),(3,1),(3,2),(3,3)], [1,2,3])
            sage: e.life_graph_save_tikz('a.tikz', ['p32','p13','ar2'])    # not tested

        ::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:    2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.life_graph_save_tikz('b.tikz', ['b12','b21','b12']) # not tested
            Creation of file b.tikz
            Using template '/Users/slabbe/.tikz2pdf.tex'.
            tikz2pdf: calling pdflatex...
            tikz2pdf: Output written to 'b.pdf'.

        """
        g = self.life_graph(substitutions)
        default_kwds = dict(format='dot2tex', edge_labels=True, color_by_label=False)
        default_kwds.update(kwds)
        g.latex_options().set_options(**default_kwds)
        tikz = latex(g)
        if tikz:
            with open(filename, 'w') as f:
                f.write(tikz)
                print "Creation of file %s" % filename
            import os
            os.system("tikz2pdf %s" % filename)

    def equivalence_class(self):
        r"""
        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:    2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: len(E.equivalence_class())
            6
        """
        L = []
        for i,j,k in itertools.permutations((1,2,3)):
            perm = WordMorphism({1:[i],2:[j],3:[k]})
            e, = [e for e in self.apply(perm) if e.is_chignons_empty()]
            L.append(e)
        return L

    def is_equivalent(self, other):
        r"""
        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:    2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.is_equivalent(E)
            True
        """
        if not isinstance(other, ExtensionType):
            return False
        else:
            return other in self.equivalence_class()

    def is_empty(self):
        if hasattr(self, '_empty'):
            return self._empty
        else:
            return False

    def is_bispecial(self):
        r"""
        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:    2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.is_bispecial()
            True
        """
        return self.right_valence() > 1 and self.left_valence() > 1
    def is_ordinaire(self):
        raise NotImplementedError

    def is_neutral(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E.is_neutral()
            True
        """
        return self.multiplicity() == 0

    def image(self, m):
        r"""

        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:          2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: b23 = WordMorphism({1:[1],2:[2,3],3:[3]})
            sage: E.image(b23)
              E(w)   1   2   3
               31        X
               12            X
               32            X
               23    X   X   X
               33    X
             m(w)=0, not ord., empty
        """
        L = []
        for e in self.apply(m):
            if e.is_chignons_empty():
                L.append(e)
        assert len(L) == 1, "len of L should be 1"
        return L[0]

    def cardinality(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E.cardinality()
            5
        """
        return len(self._pairs)
    def left_right_projection(self):
        r"""

        EXAMPLES::

            sage: L = [(1,2), (2,2), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, alphabet=(1,2,3))
            sage: E.left_right_projection()
            (Counter({3: 3, 1: 1, 2: 1}), Counter({2: 3, 1: 1, 3: 1}))

        ::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:          2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: left, right = E.left_right_projection()
            sage: sorted(left.iteritems())
            [(word: 12, 3), (word: 21, 1), (word: 22, 1), (word: 23, 1), (word: 31, 1)]
            sage: sorted(right.iteritems())
            [(word: 1, 3), (word: 2, 3), (word: 3, 1)]
        """
        left_projection = Counter()
        right_projection = Counter()
        for a,b in self:
            left_projection[a] += 1
            right_projection[b] += 1
        return left_projection, right_projection

    def left_extensions(self):
        raise NotImplementedError
    def right_extensions(self):
        raise NotImplementedError
    def left_valence(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E.left_valence()
            3

        """
        return len(self.left_extensions())
    def right_valence(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E.right_valence()
            3

        """
        return len(self.right_extensions())
    def multiplicity(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E.multiplicity()
            0

        ::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:    2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.multiplicity()
            0
        """
        return self.cardinality() - self.left_valence() - self.right_valence() + 1
class ExtensionType1to1(ExtensionType):
    r"""
    INPUT:

    - ``L`` - list of pairs of letters
    - ``alphabet`` - the alphabet
    - ``chignons`` - optional (default: None), pair of words added to the
      left  and to the right of the image of the previous bispecial

    EXAMPLES::

        sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
        sage: E = ExtensionType1to1(L, alphabet=(1,2,3))
        sage: E
          E(w)   1   2   3
           1             X
           2             X
           3     X   X   X
         m(w)=0, ordinary

    With chignons::

        sage: E = ExtensionType1to1(L, [1,2,3], ('a','b'))
        sage: E
          E(awb)   1   2   3
            1              X
            2              X
            3      X   X   X
         m(w)=0, ordinary
    """
    def __init__(self, L, alphabet, chignons=('','')):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E
              E(w)   1   2   3
               1             X
               2             X
               3     X   X   X
             m(w)=0, ordinary
        """
        self._pairs = frozenset(L)
        self._alphabet = alphabet
        self._chignons = tuple(chignons)

    def table(self):
        r"""
        return a table representation of self.

        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, alphabet=(1,2,3))
            sage: E.table()
              E(w)   1   2   3
               1             X
               2             X
               3     X   X   X

        ::

            sage: E = ExtensionType1to1(L, alphabet=(1,2,3), chignons=('a', 'b'))
            sage: E.table()
              E(awb)   1   2   3
                1              X
                2              X
                3      X   X   X
        """
        lines = []
        L = R = sorted(self._alphabet)
        for a in L:
            line = []
            for b in R:
                if (a,b) in self._pairs:
                    line.append('X')
                else:
                    line.append(' ')
            lines.append(line)
        if self._chignons != ('',''):
            Ew = "E(%sw%s)" % self._chignons
        else:
            Ew = "E(w)"
        t = table(rows=lines, header_row=R, header_column=[Ew]+L)
        t.options(header_column=False,header_row=False,align='center')
        return t

    def _repr_old(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, alphabet=(1,2,3))
            sage: print E._repr_old()
            _123
            1OOX
            2OOX
            3XXX

        With chignons::

            sage: E = ExtensionType1to1(L, [1,2,3], ('a','b'))
            sage: print E._repr_old()
            awb
            _123
            1OOX
            2OOX
            3XXX
        """
        lines = []
        if self._chignons != ('',''):
            line = "%sw%s" % self._chignons
            lines.append(line)
        line = '_' + "".join(map(str, self._alphabet))
        lines.append(line)
        for a in self._alphabet:
            line = str(a)
            for b in self._alphabet:
                if (a,b) in self._pairs:
                    line += 'X'
                else:
                    line += 'O'
            lines.append(line)
        return '\n'.join(lines)

    def _latex_old(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: print E._latex_old()
            \begin{array}{r|rrr}
            v & 1 & 2 & 3\\
            \hline
            1 &   &   & \times\\
            2 &   &   & \times\\
            3 & \times & \times & \times
            \end{array}

        ::

            sage: E = ExtensionType1to1(L, [1,2,3], ('a','b'))
            sage: print E._latex_old()
            \begin{array}{r|rrr}
            awb & 1 & 2 & 3\\
            \hline
            1 &   &   & \times\\
            2 &   &   & \times\\
            3 & \times & \times & \times
            \end{array}
        """
        lines = []
        lines.append(r"\begin{array}{r|rrr}")
        chignons = "%sw%s" % self._chignons if self._chignons != ('','') else 'v'
        lines.append(" & ".join([chignons] + map(str, self._alphabet)) + r'\\')
        lines.append(r"\hline")
        for a in self._alphabet:
            line = [str(a)]
            for b in self._alphabet:
                if (a,b) in self._pairs:
                    line.append(r'\times')
                else:
                    line += ' '
            lines.append(' & '.join(line) + r'\\')
        assert lines[-1][-2:] == r'\\'
        lines[-1] = lines[-1][:-2]
        lines.append(r"\end{array}")
        return '\n'.join(lines)

    def __eq__(self, other):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: ar = WordMorphism({1:[1,3],2:[2,3],3:[3]})
            sage: E.apply(ar) == E
            False
            sage: F = ExtensionType1to1(L, [1,2,3])
            sage: E == F
            True
        """
        if not isinstance(other, ExtensionType1to1):
            return False
        else:
            return self._pairs == other._pairs

    def chignons_multiplicity_tuple(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3], ('a', 'b'))
            sage: E.chignons_multiplicity_tuple()
            ('a', 'b', 0)
        """
        return self._chignons + (self.multiplicity(),)

    def apply(self, m):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E
              E(w)   1   2   3
               1             X
               2             X
               3     X   X   X
             m(w)=0, ordinary
            sage: ar = WordMorphism({1:[1,3],2:[2,3],3:[3]})
            sage: E.apply(ar)
            (  E(3w)   1   2   3
                1             X
                2             X
                3     X   X   X
             m(w)=0, ordinary,)

        ::

            sage: ar = WordMorphism({1:[3,1],2:[3,2],3:[3]})
            sage: E.apply(ar)
            (  E(w3)   1   2   3
                1             X
                2             X
                3     X   X   X
             m(w)=0, ordinary,)

        Creation of a pair of ordinaire bispecial words from an
        **ordinaire** word::

            sage: e = ExtensionType1to1([(1,3),(2,3),(3,1),(3,2),(3,3)], [1,2,3])
            sage: p0 = WordMorphism({1:[1,2,3],2:[2,3],3:[3]})
            sage: e.apply(p0)
            (  E(3w)   1   2   3
                1
                2             X
                3     X   X   X
             m(w)=0, ordinary,)
            sage: p3 = WordMorphism({1:[1,3,2],2:[2],3:[3,2]})
            sage: e.apply(p3)
            (  E(2w)   1   2   3
                1
                2             X
                3     X   X   X
             m(w)=0, ordinary,
               E(32w)   1   2   3
                1              X
                2      X   X   X
                3
             m(w)=0, ordinary)

        Creation of a strong-weak pair of bispecial words from a neutral
        **not ordinaire** word::

            sage: p0 = WordMorphism({1:[1,2,3],2:[2,3],3:[3]})
            sage: e = ExtensionType1to1([(1,2),(2,3),(3,1),(3,2),(3,3)], [1,2,3])
            sage: e.apply(p0)
            (  E(3w)   1   2   3
                1
                2         X   X
                3     X   X   X
             m(w)=1, not ord.,
               E(23w)   1   2   3
                1          X
                2
                3              X
             m(w)=-1, not ord.)

        Creation of a pair of ordinaire bispecial words from an **not
        ordinaire** word::

            sage: p1 = WordMorphism({1:[1,2],2:[2],3:[3,1,2]})
            sage: e = ExtensionType1to1([(1,2),(2,3),(3,1),(3,2),(3,3)], [1,2,3])
            sage: e.apply(p1)
            (  E(2w)   1   2   3
                1     X   X   X
                2             X
                3
             m(w)=0, ordinary,
               E(12w)   1   2   3
                1
                2          X
                3      X   X   X
             m(w)=0, ordinary)

        This result is now fixed::

            sage: e = ExtensionType1to1([(1,2), (3,3)], [1,2,3])
            sage: p3 = WordMorphism({1:[1,3,2],2:[2],3:[3,2]})
            sage: e.apply(p3)
            (  E(32w)   1   2   3
                1          X
                2              X
                3
             m(w)=-1, not ord.,)

        ::

            sage: e = ExtensionType1to1([(2,2),(2,3),(3,1),(3,2),(3,3)], [1,2,3])
            sage: e.apply(p3)
            (  E(2w)   1   2   3
                1
                2         X   X
                3     X   X   X
             m(w)=1, not ord.,)

        This result is now fixed::

            sage: e = ExtensionType1to1([(2,2),(2,3),(3,1),(3,2),(3,3)], [1,2,3])
            sage: p2 = WordMorphism({1:[1],2:[2,3,1],3:[3,1]})
            sage: e.apply(p2)
            (  E(31w)   1   2   3
                1      X   X   X
                2          X   X
                3
             m(w)=1, not ord.,)

        ::

            sage: e = ExtensionType1to1([(1,2),(3,3)], [1,2,3])
            sage: e.apply(p2)
            (  E(1w)   1   2   3
                1         X
                2
                3             X
             m(w)=-1, not ord.,)

        TESTS::

            sage: L = [(1,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E
              E(w)   1   2   3
               1             X
               2
               3
             m(w)=0, not ord.
            sage: ar = WordMorphism({1:[1,3],2:[2,3],3:[3]})
            sage: E.apply(ar)
            ()

        POSSIBLE BUG::

            sage: b23 = WordMorphism({1:[1],2:[2,3],3:[3]})
            sage: b13 = WordMorphism({1:[1,3],2:[2],3:[3]})
            sage: b31 = WordMorphism({1:[1],2:[2],3:[3,1]})
            sage: e = ExtensionType.from_morphism(b23)
            sage: r = e.apply(b23)[0]
            sage: r.apply(b13)
            ()
            sage: r.apply(b31)
            ()
        """
        images = m.images()
        common_suffix = images[0]
        common_prefix = images[0]
        for image in images:
            common_suffix = image.longest_common_suffix(common_suffix)
            common_prefix = image.longest_common_prefix(common_prefix)
        word_pairs = []
        for a,b in self:
            left = common_suffix * m(a)
            right = m(b) * common_prefix
            word_pairs.append( (left, right) )
        A, B = zip(*word_pairs)
        left_nb = max(map(len, A))
        right_nb = max(map(len, B))
        #print left_nb, right_nb
        L = []
        for i in range(1, left_nb+1):
            for j in range(right_nb):
                extensions = defaultdict(list)
                for left,right in word_pairs:
                    if len(left) < i or len(right) <= j:
                        continue
                    chignons = left[len(left)-i+1:], right[:j]
                    extensions[chignons].append( (left[-i], right[j]) )
                for chignons, extension in extensions.iteritems():
                    e = ExtensionType1to1(extension, self._alphabet, chignons)
                    if e.is_bispecial():
                        L.append(e)
        return tuple(L)

    def is_ordinaire(self):
        r"""
        EXAMPLES:

        ordinary::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, alphabet=(1,2,3))
            sage: E
              E(w)   1   2   3
               1             X
               2             X
               3     X   X   X
             m(w)=0, ordinary
            sage: E.is_ordinaire()
            True

        strong::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3), (1,1)]
            sage: E = ExtensionType1to1(L, alphabet=(1,2,3))
            sage: E.is_ordinaire()
            False

        neutral but not ordinary::

            sage: L = [(1,1), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, alphabet=(1,2,3))
            sage: E
              E(w)   1   2   3
               1     X
               2             X
               3     X   X   X
             m(w)=0, not ord.
            sage: E.is_neutral()
            True
            sage: E.is_ordinaire()
            False

        not neutral, not ordinaire::

            sage: L = [(1,1), (2,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, alphabet=(1,2,3))
            sage: E
              E(w)   1   2   3
               1     X
               2     X
               3         X   X
             m(w)=-1, not ord.
            sage: E.is_neutral()
            False
            sage: E.is_ordinaire()
            False
        """

        if not self.is_bispecial():
            return False
        elif not self.is_neutral():
            return False
        left_projection, right_projection = self.left_right_projection()  
        left_most = [a for (a,v) in left_projection.iteritems() if v > 1]
        if len(left_most) != 1: 
            return False
        right_most = [b for (b,v) in right_projection.iteritems() if v > 1]
        if len(right_most) != 1: 
            return False
        #print left_most,right_most
        if (left_most[0],right_most[0]) not in self._pairs:
            return False
        return True

    def cardinality(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E.cardinality()
            5
        """
        return len(self._pairs)
    def left_extensions(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E.left_extensions()
            set([1, 2, 3])

        """
        return set(a for a,b in self)
    def right_extensions(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: E.right_extensions()
            set([1, 2, 3])
        """
        return set(b for a,b in self)
class ExtensionType2to1(ExtensionType):
    r"""
    Generalized to words.

    INPUT:

    - ``L`` - list of pairs of *words*
    - ``alphabet`` - the alphabet
    - ``chignons`` - optional (default: None), pair of words added to the
      left  and to the right of the image of the previous bispecial
    - ``factors_length_2`` - list of factors of length 2. If None, they are
      computed from the provided extension assuming the bispecial factor is
      *empty*.
    - ``empty`` - bool, (optional, default: None), if None, then it is
      computed from the chignons and takes value True iff the chignons are
      empyt.

    EXAMPLES::

        sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
        ....:          2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
        sage: E = ExtensionType2to1(L, (1,2,3))
        sage: E
          E(w)   1   2   3
           21        X
           31        X
           12    X   X   X
           22    X
           23    X
        m(w)=0, not ord., empty

    """
    def __init__(self, L, alphabet, chignons=('',''),
            factors_length_2=None, empty=None):
        r"""
        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:          2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
        """
        self._pairs = frozenset((Word(a),Word(b)) for a,b in L)
        self._alphabet = alphabet
        self._chignons = tuple(chignons)
        self._factors_length_2 = factors_length_2
        if empty is None:
            self._empty = self.is_chignons_empty()
        else:
            self._empty = empty

    def is_chignons_empty(self):
        return map(len, self._chignons) == [0,0]

    def factors_length_2(self):
        r"""
        Returns the set of factors of length 2 of the language.

        This is computed from the extension type if it was not provided at
        the construction.

        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:    2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.factors_length_2()
            set([(1, 2), (3, 1), (2, 3), (2, 1), (2, 2)])
        """
        if self._factors_length_2 is None:
            # We suppose here the factor is the empty word...
            self._factors_length_2 = set((a[-1],b[0]) for a,b in self._pairs)
        return self._factors_length_2

    def is_valid(self):
        if any(len(a)==0 or len(b)==0 for a,b in self._pairs):
            return False
        return True

    def left_word_extensions(self):
        r"""
        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:          2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: sorted(E.left_word_extensions())
            [word: 12, word: 21, word: 22, word: 23, word: 31]

        """
        return set(a for a,b in self)
    def right_word_extensions(self):
        r"""
        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:          2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: sorted(E.right_word_extensions())
            [word: 1, word: 2, word: 3]
        """
        return set(b for a,b in self)
    def table(self):
        r"""
        return a table representation of self.

        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:          2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E
              E(w)   1   2   3
               21        X
               31        X
               12    X   X   X
               22    X
               23    X
             m(w)=0, not ord., empty
        """
        lines = []
        L = sorted(self.left_word_extensions(), key=lambda w:w.reversal())
        R = sorted(self.right_word_extensions())
        for a in L:
            line = []
            for b in R:
                if (a,b) in self._pairs:
                    line.append('X')
                else:
                    line.append(' ')
            lines.append(line)
        if self._chignons != ('',''):
            Ew = "E(%sw%s)" % self._chignons
        else:
            Ew = "E(w)"
        t = table(rows=lines, header_row=R, header_column=[Ew]+L)
        t.options(header_column=False,header_row=False,align='center')
        return t

    def __eq__(self, other):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3])
            sage: ar = WordMorphism({1:[1,3],2:[2,3],3:[3]})
            sage: E.apply(ar) == E
            False
            sage: F = ExtensionType1to1(L, [1,2,3])
            sage: E == F
            True
        """
        if not isinstance(other, ExtensionType2to1):
            return False
        else:
            return self._pairs == other._pairs

    def chignons_multiplicity_tuple(self):
        r"""
        EXAMPLES::

            sage: L = [(1,3), (2,3), (3,1), (3,2), (3,3)]
            sage: E = ExtensionType1to1(L, [1,2,3], ('a', 'b'))
            sage: E.chignons_multiplicity_tuple()
            ('a', 'b', 0)
        """
        return self._chignons + (self.multiplicity(),)

    def letters_before_and_after(self):
        r"""
        Returns a pair of dict giving the possible letters that goes before
        or after a letter.

        Computed from the set of factors of length 2 of the language.

        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:      2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.factors_length_2()
            set([(1, 2), (3, 1), (2, 3), (2, 1), (2, 2)])
            sage: E.letters_before_and_after()
            (defaultdict(<type 'set'>, {1: set([2, 3]), 2: set([1, 2]), 3: set([2])}),
            defaultdict(<type 'set'>, {1: set([2]), 2: set([1, 2, 3]), 3: set([1])}))
        """
        possible_before  = defaultdict(set)
        possible_after = defaultdict(set)
        for x,y in self.factors_length_2():
            possible_before[y].add(x)
            possible_after[x].add(y)
        return possible_before, possible_after

    def apply(self, m, l=2, r=1):
        r"""
        The code works for Brun here because we take length 2 on the left
        and length 1 on the right.

        On utilise les facteurs de longueur 2 pour completer l'info qui
        peut manquer.

        TODO: bien corriger les facteurs de longueurs 2 de l'image!!!

        INPUT:

        - ``m`` - substitution
        - ``l`` - integer, length of left extension
        - ``r`` - integer, length of right extension

        OUTPUT:

        list of Extension type of the bispecial images

        POSSIBLE BUG::

            sage: b23 = WordMorphism({1:[1],2:[2,3],3:[3]})
            sage: b13 = WordMorphism({1:[1,3],2:[2],3:[3]})
            sage: b31 = WordMorphism({1:[1],2:[2],3:[3,1]})
            sage: e = ExtensionType.from_morphism(b23)
            sage: r = e.apply(b23)[0]
            sage: r.apply(b13)
            ()
            sage: r.apply(b31)
            ()

        On a le meme bug (ca se corrige avec de plus grandes extensions a gauche)::

            sage: E = ExtensionType2to1((([a],[b]) for a,b in e), (1,2,3))
            sage: E.apply(b23)[0].apply(b13)
            (  E(w)   1   2   3
               31            X
               32            X
               3     X   X   X
               13    X   X   X
               23            X
             m(w)=0, ordinary, empty,)
            sage: E.apply(b23)[0].apply(b31)
            (  E(1w)   1   2   3
                1     X   X   X
                3     X   X   X
               23             X
             m(w)=2, not ord.,   
               E(w)   1   2   3
               11    X   X   X
               31    X   X   X
               12            X
               13    X
               23    X
             m(w)=0, not ord., empty)

        EXAMPLES:

        On imagine qu'on vient de faire b12::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:      2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E
              E(w)   1   2   3
               21        X
               31        X
               12    X   X   X
               22    X
               23    X
             m(w)=0, not ord., empty
            sage: b12 = WordMorphism({1:[1,2],2:[2],3:[3]})
            sage: E.apply(b12)
            (  E(2w)   1   2   3
               21         X
               31         X
               12     X   X   X
               22     X
            m(w)=0, ordinary,
               E(w)   1   2   3
               21        X
               31        X
               12        X
               22    X   X   X
               23    X
            m(w)=0, not ord., empty)

        ::

            sage: b21 = WordMorphism({1:[1],2:[2,1],3:[3]})
            sage: E.apply(b21)
            (  E(1w)   1   2   3
               21         X
               12     X   X   X
               13         X
             m(w)=0, ordinary,
               E(w)   1   2   3
               11        X
               21    X   X   X
               31        X
               12    X
               13    X
             m(w)=0, ordinary, empty)
            sage: b23 = WordMorphism({1:[1],2:[2,3],3:[3]})
            sage: E.apply(b23)
            (  E(3w)   1   2   3
               12     X   X   X
               32     X
               23     X
             m(w)=0, ordinary,
               E(23w)   1   2   3
                31     X   X   X
                23     X
             m(w)=0, ordinary,
               E(w)   1   2   3
               31        X
               12            X
               32            X
               23    X   X   X
               33    X
             m(w)=0, not ord., empty)

        """
        letters_before, letters_after = self.letters_before_and_after()
        #print letters_before, letters_after
        word_before = defaultdict(Word)
        word_after = defaultdict(Word)
        for key,value in letters_before.iteritems():
            word_before[key] = longest_common_suffix(map(m, value))
        for key,value in letters_after.iteritems():
            word_after[key] = longest_common_prefix(map(m, value))
        #print word_before, word_after

        extensions = defaultdict(list)
        for a,b in self:
            left = word_before[a[0]] * m(a)
            right = m(b) * word_after[b[-1]]
            length_image_last_a = len(m(a[-1]))
            length_image_first_b = len(m(b[0]))
            for i in range(length_image_last_a+1):
                for j in range(length_image_first_b+1):
                    chignons = left[len(left)-i:], right[:j]
                    new_ext  = left[:len(left)-i][-l:], right[j:j+r]
                    extensions[chignons].append( new_ext )

        # The empty word (as image of the empty word) occurs in every place...
        if self.is_empty():
            chignons = Word(),Word()
            assert chignons in extensions
            for a,b in self:
                left = word_before[a[0]] * m(a)
                right = m(b) * word_after[b[-1]]
                word = left * right
                for f in word.factor_iterator(l+r):
                    new_ext = f[:l], f[-r:]
                    extensions[chignons].append( new_ext )

        F = factors_length_2_from_morphism_and_factors_length_2(m,
                self.factors_length_2())
        L = []
        for chignons, extension in extensions.iteritems():
            empty = self.is_empty() and map(len, chignons) == [0,0]
            e = ExtensionType2to1(extension, self._alphabet, chignons, F, empty)
            if e.is_valid() and e.is_bispecial():
                L.append(e)
        return tuple(L)

    @cached_method
    def extension_type_1to1(self):
        r"""
        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:    2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.extension_type_1to1()
              E(w)   1   2   3
               1         X
               2     X   X   X
               3     X
             m(w)=0, not ord.

        """
        pairs = set((a[-1],b[0]) for a,b in self)
        return ExtensionType1to1(pairs, alphabet=self._alphabet,
                chignons=self._chignons)

    def cardinality(self):
        r"""
        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:    2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.cardinality()
            5

        """
        return self.extension_type_1to1().cardinality()
    def is_ordinaire(self):
        r"""
        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:          2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.extension_type_1to1()
              E(w)   1   2   3
               1         X
               2     X   X   X
               3     X
             m(w)=0, not ord.
            sage: E.is_ordinaire()
            False
        """
        return self.extension_type_1to1().is_ordinaire()

    def left_extensions(self):
        r"""
        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:    2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.left_extensions()
            set([1, 2, 3])

        """
        return self.extension_type_1to1().left_extensions()
        #return set(a[-1:] for a in self.left_word_extensions())
    def right_extensions(self):
        r"""
        EXAMPLES::

            sage: L = [((2, 2), (1,)), ((2, 3), (1,)), ((2, 1), (2,)), ((1,
            ....:    2), (1,)), ((1, 2), (2,)), ((1, 2), (3,)), ((3, 1), (2,))]
            sage: E = ExtensionType2to1(L, (1,2,3))
            sage: E.right_extensions()
            set([1, 2, 3])
        """
        return self.extension_type_1to1().right_extensions()
        #return set(b[:1] for b in self.right_word_extensions())
######################################
# methods that should be in Sage
######################################
def longest_common_prefix(L):
    r"""
    Return the longest common prefix of a list of words.

    EXAMPLES::

        sage: longest_common_prefix((Word('ab'), Word('abc'), Word('abd')))
        word: ab
    """
    common = L[0]
    for w in L:
        common = w.longest_common_prefix(common)
    return common
def longest_common_suffix(L):
    r"""
    Return the longest common suffix of a list of words.
    """
    common = L[0]
    for w in L:
        common = w.longest_common_suffix(common)
    return common
