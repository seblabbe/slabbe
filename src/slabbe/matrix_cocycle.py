# -*- coding: utf-8 -*-
r"""
Multidimensional Continued Fraction Algorithm's cocyles

EXAMPLES:

The 1-cylinders of ARP transformation as matrices::

    sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
    sage: ARP = mcf_cocycles.ARP()
    sage: zip(*ARP.n_cylinders_iterator(1))
    [(word: A1,
      word: A2,
      word: A3,
      word: P12,
      word: P13,
      word: P21,
      word: P23,
      word: P31,
      word: P32),
     (
    [1 1 1]  [1 0 0]  [1 0 0]  [1 1 1]  [1 1 1]  [2 1 1]  [1 0 1]  [2 1 1]
    [0 1 0]  [1 1 1]  [0 1 0]  [1 2 1]  [0 1 1]  [1 1 1]  [1 1 1]  [1 1 0]
    [0 0 1], [0 0 1], [1 1 1], [0 1 1], [1 1 2], [1 0 1], [1 1 2], [1 1 1],
    <BLANKLINE>
    [1 1 0]
    [1 2 1]
    [1 1 1]
    )]

Ces calculs illustrent le bounded distorsion de ratio=4 pour ARP
multiplicatif (2 avril 2014)::

    sage: T = mcf_cocycles.Sorted_ARPMulti(2)
    sage: T.distorsion_max(1, p=oo)
    5
    sage: T.distorsion_max(2, p=oo)
    7
    sage: T.distorsion_max(3, p=oo)
    22/3
    sage: T.distorsion_max(4, p=oo)    # long time (4s)
    62/17

::

    sage: T = mcf_cocycles.Sorted_ARPMulti(3)
    sage: T.distorsion_max(1, p=oo)
    7
    sage: T.distorsion_max(2, p=oo)
    9
    sage: T.distorsion_max(3, p=oo)
    19/2
    sage: T.distorsion_max(4, p=oo)  # long time (47s)
    161/43

"""
import itertools
from collections import Counter
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from language import Language, FiniteLanguage

######################
# class MCF_Cocycle
######################
class MCF_Cocycle(object):
    r"""
    Multidimensional continued fraction algorithm cocycle

    INPUT:

    - ``gens`` -- list, tuple or dict; the matrices. Keys 0,...,n-1 are
      used for list and tuple.
    - ``cone`` -- dict or matrix; the cone for each matrix generators.
      If it is a matrix, then it serves as the cone for all matrices. The
      cone is defined by the columns of the matrix.
    - ``language`` -- regular language or None; if None, the language is
      the full shift.

    EXAMPLES::

        sage: from slabbe.mult_cont_frac_cocycle import MCF_Cocycle
        sage: B1 = matrix(3, [1,0,0, 0,1,0, 0,1,1])
        sage: B2 = matrix(3, [1,0,0, 0,0,1, 0,1,1])
        sage: B3 = matrix(3, [0,1,0, 0,0,1, 1,0,1])
        sage: gens = {'B1':B1, 'B2':B2, 'B3':B3}
        sage: cone = matrix(3, [1,1,1,0,1,1,0,0,1])
        sage: MCF_Cocycle(gens, cone)
        Cocycle with 3 gens over Language of finite words over alphabet ['B1', 'B2', 'B3']
    """
    def __init__(self, gens, cone, language=None):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import MCF_Cocycle
            sage: gens = {'A':matrix(3, [1,0,0, 0,1,0, 0,1,1])}
            sage: cone = identity_matrix(3)
            sage: MCF_Cocycle(gens, cone)
            Cocycle with 1 gens over Language of finite words over alphabet ['A']
        """
        if isinstance(gens, dict):
            self._gens = gens
        elif isinstance(gens, (list, tuple)):
            self._gens = dict(enumerate(gens))
        else:
            raise ValueError("gens must be a list, tuple or a dict")
        if isinstance(cone, dict):
            self._cone_dict = cone
        else:
            self._cone_dict = {letter:cone for letter in self._gens.keys()}
        if language is None:
            self._language = Language(sorted(self._gens.keys()))
        else:
            self._language = language

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import MCF_Cocycle
            sage: gens = {'A':matrix(3, [1,0,0, 0,1,0, 0,1,1])}
            sage: cone = identity_matrix(3)
            sage: MCF_Cocycle(gens, cone)
            Cocycle with 1 gens over Language of finite words over alphabet ['A']
        """
        s = "Cocycle with {} gens over {}"
        return s.format(len(self._gens), self._language)

    def gens(self):
        return self._gens
    def cone_dict(self):
        return self._cone_dict
    def cone(self, key):
        return self._cone_dict[key]
    def language(self):
        return self._language

    @cached_method
    def identity_matrix(self):
        return self._gens.values()[0].parent().one()

    def word_to_matrix(self, w):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: C = mcf_cocycles.Sorted_ARP()
            sage: C.word_to_matrix(Word())
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return prod((self._gens[a] for a in w), z=self.identity_matrix())

    def n_words_iterator(self, n):
        r"""
        EXAMPLES::
            
            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: ARP = mcf_cocycles.Sorted_ARP()
            sage: list(ARP.n_words_iterator(1))
            [word: A1, word: A2, word: A3, word: P1, word: P2, word: P3]
        """
        return self._language.words_of_length_iterator(n)

    def n_matrices_iterator(self, n):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: ARP = mcf_cocycles.Sorted_ARP()
            sage: A,B = zip(*list(ARP.n_matrices_iterator(1)))
            sage: A
            (word: A1, word: A2, word: A3, word: P1, word: P2, word: P3)
            sage: B
            (
            [1 0 0]  [1 0 0]  [0 1 0]  [0 1 0]  [0 0 1]  [0 0 1]
            [0 1 0]  [0 0 1]  [0 0 1]  [0 1 1]  [1 0 1]  [0 1 1]
            [1 1 1], [1 1 1], [1 1 1], [1 1 1], [1 1 1], [1 1 1]
            )
        """
        for w in self.n_words_iterator(n):
            yield w, self.word_to_matrix(w)

    def n_matrices_eigenvalues_iterator(self,n):
        r"""
        Return the eigenvalues of the matrices of level n.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: ARP = mcf_cocycles.ARP()
            sage: list(ARP.n_matrices_eigenvalues_iterator(1))
            [(word: A1, [1, 1, 1]),
             (word: A2, [1, 1, 1]),
             (word: A3, [1, 1, 1]),
             (word: P12, [1, 1, 1]),
             (word: P13, [1, 1, 1]),
             (word: P21, [1, 1, 1]),
             (word: P23, [1, 1, 1]),
             (word: P31, [1, 1, 1]),
             (word: P32, [1, 1, 1])]

        ::

            sage: B = mcf_cocycles.Sorted_Brun()
            sage: list(B.n_matrices_eigenvalues_iterator(1))
            [(word: B1, [1, 1, 1]),
             (word: B2, [1, -0.618033988749895?, 1.618033988749895?]),
             (word: B3, [1.465571231876768?, -0.2327856159383841? -
                         0.7925519925154479?*I, -0.2327856159383841? +
                         0.7925519925154479?*I])]
        """
        for w,m in self.n_matrices_iterator(n):
            yield w, m.eigenvalues()

    def n_matrices_eigenvectors(self,n, verbose=False):
        r"""
        Return the left and right eigenvectors of the matrices of level n.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: C = mcf_cocycles.ARP()
            sage: C.n_matrices_eigenvectors(1)
            [(word: A1, (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
             (word: A2, (0.0, 1.0, 0.0), (1.0, 0.0, 0.0)),
             (word: A3, (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
             (word: P12, (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
             (word: P13, (0.0, 0.0, 1.0), (0.0, 1.0, 0.0)),
             (word: P21, (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
             (word: P23, (0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
             (word: P31, (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
             (word: P32, (0.0, 1.0, 0.0), (1.0, 0.0, 0.0))]
        """
        R = []
        for w,m in self.n_matrices_iterator(n):
            try:
                a,v_right = perron_right_eigenvector(m)
                b,v_left = perron_right_eigenvector(m.transpose())
            except ValueError:
                print "problem with :\n",m
            else:
                R.append((w, v_right,v_left))
                if verbose:
                    print "indices of matrices:", w
                    print m
                    print "eigenvectors:", v_right, v_left
        return R
    @cached_method
    def n_matrices_non_pisot(self,n, verbose=False):
        r"""
        Return the list of non pisot matrices (as list of indices of base
        matrices).

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: ARP = mcf_cocycles.Sorted_ARP()
            sage: ARP.n_matrices_non_pisot(1)
            [word: A1, word: A2]
            sage: ARP.n_matrices_non_pisot(2)   # long time (1s)
            [word: A1,A1, word: A1,A2, word: A2,A1, word: A2,A2]
            sage: ARP.n_matrices_non_pisot(3)   # long time (11s)
            [word: A1,A1,A1,
             word: A1,A1,A2,
             word: A1,A2,A1,
             word: A1,A2,A2,
             word: A2,A1,A1,
             word: A2,A1,A2,
             word: A2,A2,A1,
             word: A2,A2,A2]
            sage: len(ARP.n_matrices_non_pisot(4))  # long time
            16

        ::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: B = mcf_cocycles.Sorted_Brun()
            sage: B.n_matrices_non_pisot(2)
            [word: B1,B1, word: B1,B2, word: B2,B1, word: B2,B2]
            sage: B.n_matrices_non_pisot(3)
            [word: B1,B1,B1,
             word: B1,B1,B2,
             word: B1,B2,B1,
             word: B1,B2,B2,
             word: B2,B1,B1,
             word: B2,B1,B2,
             word: B2,B2,B1,
             word: B2,B2,B2]
        """
        return [w for w in self.n_words_iterator(n) if not self.is_pisot(w)]

    def n_matrices_semi_norm_iterator(self, n, p=2):
        r"""
        EXAMPLES:

        For the 1-norm, all matrices contracts the hyperplane::
            
            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: C = mcf_cocycles.ARP()
            sage: it = C.n_matrices_semi_norm_iterator(1, p=1)
            sage: for _ in range(5): print next(it) # tolerance 0.00001
            (word: A1, 1.0, False)
            (word: A2, 1.0, False)
            (word: A3, 1.0, False)
            (word: P12, 0.9999943584528624, False)
            (word: P13, 0.9999937410436893, False)

        For the 2-norm, AR matrices do not contract::

            sage: it = C.n_matrices_semi_norm_iterator(1, p=2)
            sage: for w,s,b in it: print w,s,b   # long time (6s)
            A1 1.30656296488 False
            A2 1.30656296486 False
            A3 1.30656296475 False
            P12 0.99999999996 False
            P13 0.999999999967 False
            P21 0.999999999967 False
            P23 0.999999999997 False
            P31 0.999999999769 False
            P32 0.999999999839 False

        When, the 1-norm is < 1, the product is pisot::

            sage: it = C.n_matrices_semi_norm_iterator(2, p=1)
            sage: for w,s,b in it: print w,s,b   # long time
            A1,A1 1.0 False
            A1,A2 1.0 False
            A1,A3 1.0 False
            A1,P12 0.999998922557 False
            A1,P13 0.999997464905 False
            A1,P21 0.999993244882 False
            A1,P23 0.999999150973 True
            A1,P31 0.999994030522 False
            A1,P32 0.999998046513 True
            A2,A1 1.0 False
            A2,A2 1.0 False
            A2,A3 1.0 False
            A2,P12 0.99999375291 False
            A2,P13 0.999995591588 True
            ...
            P31,A3 0.999988326888 False
            P31,P12 0.749998931902 True
            P31,P23 0.799999157344 True
            P31,P32 0.749993104833 True
            P32,A1 0.999997170005 True
            P32,A3 0.99999420509 False
            P32,P13 0.666665046248 True
            P32,P21 0.666665629351 True
            P32,P31 0.666664488371 True
        """
        if n == 0:
            raise NotImplementedError
        for w,m in self.n_matrices_iterator(n):
            cone = m*self.cone(w[-1])
            yield w, semi_norm_cone(m.transpose(), cone, p=p), self.is_pisot(w)

    def n_matrices_distorsion_iterator(self, n, p=1):
        r"""
        Return the the distorsion of the n-cylinders.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: T = mcf_cocycles.Sorted_ARP()
            sage: it =T.n_matrices_distorsion_iterator(1)
            sage: list(it)
            [(word: A1, 2),
             (word: A2, 2),
             (word: A3, 2),
             (word: P1, 3),
             (word: P2, 3),
             (word: P3, 3)]
        """
        for w,m in self.n_matrices_iterator(n):
            yield w, distorsion(m, p=p)

    def n_cylinders_iterator(self, n):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: C = mcf_cocycles.ARP()
            sage: it = C.n_cylinders_iterator(1)
            sage: for w,cyl in it: print "{}\n{}".format(w,cyl)
            A1
            [1 1 1]
            [0 1 0]
            [0 0 1]
            A2
            [1 0 0]
            [1 1 1]
            [0 0 1]
            A3
            [1 0 0]
            [0 1 0]
            [1 1 1]
            P12
            [1 1 1]
            [1 2 1]
            [0 1 1]
            P13
            [1 1 1]
            [0 1 1]
            [1 1 2]
            P21
            [2 1 1]
            [1 1 1]
            [1 0 1]
            P23
            [1 0 1]
            [1 1 1]
            [1 1 2]
            P31
            [2 1 1]
            [1 1 0]
            [1 1 1]
            P32
            [1 1 0]
            [1 2 1]
            [1 1 1]
        """
        if n == 0:
            raise NotImplementedError
        for w,m in self.n_matrices_iterator(n):
            yield w, m*self.cone(w[-1])

    def n_cylinders_edges(self, n):
        r"""
        Return the set of edges of the n-cylinders.

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: ARP = mcf_cocycles.ARP()
            sage: ARP.n_cylinders_edges(1)
            {frozenset({(1, 1, 0), (1, 1, 1)}),
             frozenset({(0, 1, 0), (1, 1, 0)}),
             frozenset({(1, 1, 1), (2, 1, 1)}),
             frozenset({(0, 0, 1), (1, 0, 1)}),
             frozenset({(0, 1, 0), (0, 1, 1)}),
             frozenset({(0, 1, 1), (1, 0, 1)}),
             frozenset({(1, 0, 0), (1, 1, 0)}),
             frozenset({(0, 1, 1), (1, 2, 1)}),
             frozenset({(0, 1, 1), (1, 1, 2)}),
             frozenset({(1, 1, 1), (1, 1, 2)}),
             frozenset({(1, 0, 1), (1, 1, 2)}),
             frozenset({(1, 0, 1), (2, 1, 1)}),
             frozenset({(0, 0, 1), (0, 1, 1)}),
             frozenset({(1, 1, 0), (1, 2, 1)}),
             frozenset({(1, 1, 0), (2, 1, 1)}),
             frozenset({(1, 0, 0), (1, 0, 1)}),
             frozenset({(1, 1, 1), (1, 2, 1)}),
             frozenset({(1, 0, 1), (1, 1, 0)}),
             frozenset({(0, 1, 1), (1, 1, 1)}),
             frozenset({(0, 1, 1), (1, 1, 0)}),
             frozenset({(1, 0, 1), (1, 1, 1)})}
        """
        from sage.rings.finite_rings.integer_mod_ring import Integers
        edges = set()
        for w,cyl in self.n_cylinders_iterator(n):
            cols = cyl.columns()
            indices = Integers(len(cols))
            edges.update(frozenset((cols[i], cols[i+1])) for i in indices)
        return edges

    def is_pisot(self, w):
        r"""
        """
        m = self.word_to_matrix(w)
        S = sorted((abs(e) for e in m.eigenvalues()), reverse=True)
        return S[0] > 1 and S[1] < 1

    def non_pisot_automaton(self, n):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: C = mcf_cocycles.ARP()
            sage: A = C.non_pisot_automaton(2)
            sage: A
            Automaton with 2 states
            sage: A.graph().plot(edge_labels=True)   # not tested
        """
        L = []
        for i in range(n):
            L.extend(self.n_matrices_non_pisot(i))
        alphabet = self._language._alphabet
        F = FiniteLanguage(alphabet, L)
        A = F.minimal_automaton()  
        return A
        #G = A.graph()
        #to_remove = set(A.states()) - set(A.final_states())
        #G.delete_vertices(to_remove)
        #return G

    def distorsion_max(self, n, p=1):
        r"""
        EXAMPLES:

        Non borné::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: T = mcf_cocycles.Sorted_ARP()
            sage: T.distorsion_max(1, p=oo)
            1
            sage: T.distorsion_max(2, p=oo)
            3
            sage: T.distorsion_max(3, p=oo)
            5
            sage: T.distorsion_max(4, p=oo)
            7
        """
        return max(d for (w,d) in self.n_matrices_distorsion_iterator(n, p=p))

    def distorsion_argmax(self, n, p=1):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: ARP = mcf_cocycles.Sorted_ARP()
            sage: ARP.distorsion_argmax(1)
            (
                      [1 0 0]
                      [1 1 0]
            word: A1, [3 2 1]
            )
        """
        it = self.n_cylinders_iterator(n)
        key = lambda (w,m):distorsion(m, p=p)
        return max(it, key=key)

    def plot_n_cylinders(self, n, labels=True):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: C = mcf_cocycles.Sorted_ARP()
            sage: G = C.plot_n_cylinders(3)
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

    def plot_n_matrices_eigenvectors(self, n, side='right', color_index=0, draw_line=False):
        r"""
        INPUT:

        - ``n` -- integer, length
        - ``side`` -- ``'left'`` or ``'right'``, drawing left or right
          eigenvectors
        - ``color_index`` -- 0 for first letter, -1 for last letter
        - ``draw_line`` -- boolean

        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: ARP = mcf_cocycles.ARP()
            sage: G = ARP.plot_n_matrices_eigenvectors(2)
        """
        from sage.plot.graphics import Graphics
        from sage.plot.point import point
        from sage.plot.line import line
        from sage.plot.text import text
        from sage.plot.colors import hue
        from sage.modules.free_module_element import vector
        M3to2 = projection_matrix(3, 2)
        R = self.n_matrices_eigenvectors(n)
        L = [(w, M3to2*(a/sum(a)), M3to2*(b/sum(b))) for (w,a,b) in R]
        G = Graphics()
        alphabet = self._language._alphabet
        color_ = dict( (letter, hue(i/float(len(alphabet)))) for i,letter in
                enumerate(alphabet))
        for letter in alphabet:
            L_filtered = [(w,p1,p2) for (w,p1,p2) in L if w[color_index] == letter]
            words,rights,lefts = zip(*L_filtered)
            if side == 'right':
                G += point(rights, color=color_[letter], legend_label=letter)
            elif side == 'left':
                G += point(lefts,  color=color_[letter], legend_label=letter)
            else:
                raise ValueError("side(=%s) should be left or right" % side)

        if draw_line:
            for (a,b) in L:
                G += line([a,b], color='black', linestyle=":")
        G += line([M3to2*vector(a) for a in [(1,0,0), (0,1,0), (0,0,1), (1,0,0)]]) 
        title = "%s eigenvectors, colored by letter w[%s] of cylinder w" % (side, color_index)
        G += text(title, (0.5, 1.05), axis_coords=True)
        G.axes(False)
        return G

    def plot_pisot_conjugates(self, n):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: B = mcf_cocycles.Sorted_Brun()
            sage: G = B.plot_pisot_conjugates(5)   # long time (8s)

        Image envoyee a Timo (6 mai 2014)::

            sage: G = sum(B.plot_pisot_conjugates(i) for i in [1..6])  #not tested
        """
        from sage.plot.point import points
        Lreal = []
        Limag = []
        for w,s in self.n_matrices_eigenvalues_iterator(n):
            a,b,c = sorted(s, key=abs)
            if a.imag() == 0 and b.imag() == 0:
                Lreal.append((a,b))
            else:
                Limag.append((a.real(),a.imag()))
                Limag.append((b.real(),b.imag()))
        return points(Lreal) + points(Limag, color='red')

    def tikz_n_cylinders(self, n, labels=False, scale=1):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: ARP = mcf_cocycles.ARP()
            sage: s = ARP.tikz_n_cylinders(1, labels=True, scale=4)
            sage: view(s, tightpage=True)   # not tested
        """
        from sage.misc.latex import LatexExpr
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
        return LatexExpr("\n".join(lines))

####################
# quick construction
####################
class MCFCocycleGenerator(object):
    def ARP(self):
        from sage.matrix.constructor import identity_matrix
        A1 = matrix(3, [1,-1,-1, 0,1,0, 0,0,1]).inverse()
        A2 = matrix(3, [1,0,0, -1,1,-1, 0,0,1]).inverse()
        A3 = matrix(3, [1,0,0, 0,1,0, -1,-1,1]).inverse()
        P12 = matrix(3, [1, 0, 1, 1, 1, 1, 0, 0, 1])
        P13 = matrix(3, [1, 1, 0, 0, 1, 0, 1, 1, 1])
        P23 = matrix(3, [1, 0, 0, 1, 1, 0, 1, 1, 1])
        P21 = matrix(3, [1, 1, 1, 0, 1, 1, 0, 0, 1])
        P31 = matrix(3, [1, 1, 1, 0, 1, 0, 0, 1, 1])
        P32 = matrix(3, [1, 0, 0, 1, 1, 1, 1, 0, 1])
        gens = (A1, A2, A3, P12, P13, P23, P21, P31, P32)
        alphabet = ['A1', 'A2', 'A3', 'P12', 'P13', 'P23', 'P21', 'P31', 'P32']
        gens = dict(zip(alphabet, gens))

        pairs = [(1,2), (1,3), (2,3), (2,1), (3,2), (3,1)]
        cone = {"P{}{}".format(j,k):self._ARP_H_matrices(j,k) for (j,k) in pairs}
        cone.update({"A{}".format(k):identity_matrix(3) for k in [1,2,3]})

        from language import languages
        return MCF_Cocycle(gens, cone, language=languages.ARP())

    def _ARP_H_matrices(self, j,k, normalize=False):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac_cocycle import mcf_cocycles
            sage: possible = [(3,1), (2,1), (3,2), (1,2), (2,3), (1,3)]
            sage: [mcf_cocycles._ARP_H_matrices(j,k) for (j,k) in possible]
            [
            [1 0 0]  [1 0 0]  [1 1 0]  [1 0 0]  [1 0 1]  [1 0 0]
            [1 1 0]  [0 1 0]  [0 1 0]  [0 1 0]  [0 1 0]  [0 1 1]
            [0 0 1], [1 0 1], [0 0 1], [0 1 1], [0 0 1], [0 0 1]
            ]
        """
        if normalize:
            a = 1/2
        else:
            a = 1
        if j==3 and k==1:
            return matrix([(a,a,0), (0,1,0), (0,0,1)]).transpose()
        elif j==2 and k==1:
            return matrix([(a,0,a), (0,1,0), (0,0,1)]).transpose()
        elif j==3 and k==2:
            return matrix([(1,0,0), (a,a,0), (0,0,1)]).transpose()
        elif j==1 and k==2:
            return matrix([(1,0,0), (0,a,a), (0,0,1)]).transpose()
        elif j==2 and k==3:
            return matrix([(1,0,0), (0,1,0), (a,0,a)]).transpose()
        elif j==1 and k==3:
            return matrix([(1,0,0), (0,1,0), (0,a,a)]).transpose()
        else:
            raise ValueError("invalid j(=%s) and k(=%s)" % (j,k))


    def ArnouxRauzy(self):
        A1 = matrix(3, [1,-1,-1, 0,1,0, 0,0,1]).inverse()
        A2 = matrix(3, [1,0,0, -1,1,-1, 0,0,1]).inverse()
        A3 = matrix(3, [1,0,0, 0,1,0, -1,-1,1]).inverse()
        gens = (A1, A2, A3)
        alphabet = ['A1', 'A2', 'A3']
        gens = dict(zip(alphabet, gens))
        cone = matrix(3, [1,0,0,0,1,0,0,0,1])
        return MCF_Cocycle(gens, cone)

    def Sorted_Brun(self):
        B1 = matrix(3, [1,0,0, 0,1,0, 0,-1,1]).inverse()
        B2 = matrix(3, [1,0,0, 0,-1,1, 0,1,0]).inverse()
        B3 = matrix(3, [0,-1,1, 1,0,0, 0,1,0]).inverse()
        gens = (B1, B2, B3)
        alphabet = ['B1', 'B2', 'B3']
        gens = dict(zip(alphabet, gens))
        cone = matrix(3, [1,1,1,0,1,1,0,0,1])
        return MCF_Cocycle(gens, cone)

    def Sorted_ARP(self):
        A1 = matrix(3, [1,0,0,0,1,0,-1,-1,1]).inverse()
        A2 = matrix(3, [1,0,0,-1,-1,1,0,1,0]).inverse()
        A3 = matrix(3, [-1,-1,1,1,0,0,0,1,0]).inverse()
        P1 = matrix(3, [0,-1,1,1,0,0,-1,1,0]).inverse()
        P2 = matrix(3, [-1,1,0,0,-1,1,1,0,0]).inverse()
        P3 = matrix(3, [0,-1,1,-1,1,0,1,0,0]).inverse()
        gens = (A1, A2, A3, P1, P2, P3)
        alphabet = ['A1', 'A2', 'A3', 'P1', 'P2', 'P3']
        gens = dict(zip(alphabet, gens))
        cone = matrix(3, [1,0,0,1,1,0,1,1,1])
        return MCF_Cocycle(gens, cone)

    def Sorted_ARPMulti(self, order):
        A1 = matrix(3, [1,0,0,0,1,0,-1,-1,1]).inverse()
        P1 = matrix(3, [0,-1,1,1,0,0,-1,1,0]).inverse()
        P2 = matrix(3, [-1,1,0,0,-1,1,1,0,0]).inverse()
        P3 = matrix(3, [0,-1,1,-1,1,0,1,0,0]).inverse()
        cone = matrix(3, [1,0,0,1,1,0,1,1,1])
        gens = [P1,P2,P3]
        t12 = matrix(3, [1,0,0,0,0,1,0,1,0])
        t132 = matrix(3, [0,1,0,0,0,1,1,0,0])
        for i in range(1, order+1):
            A1_i = A1**i
            gens.append(A1_i * P1)
            gens.append(A1_i * P2)
            gens.append(A1_i * P3)
            gens.append(A1_i * t12)
            gens.append(A1_i * t132)
        return MCF_Cocycle(gens, cone)

mcf_cocycles = MCFCocycleGenerator()

#####################
# Helper functions
#####################
def rounded_string_vector(v, digits=4):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_cocycle import rounded_string_vector
        sage: v = (-0.144337567297406, 0.166666666666667)
        sage: rounded_string_vector(v)
        '(-0.1443, 0.1667)'
        sage: rounded_string_vector(v, digits=6)
        '(-0.144338, 0.166667)'
    """
    s = "{{:.{}f}}".format(digits)
    content = ", ".join(s.format(float(a)) for a in v)
    return "({})".format(content)

def projection_matrix(dim_from=3, dim_to=2):
    r"""
    Return a projection matrix from R^d to R^l.

    INPUT:

    - ``dim_from` -- integer (default: ``3``)
    - ``dim_to` -- integer (default: ``2``)

    OUTPUT:

        matrix

    EXAMPLES::

        sage: from slabbe.mult_cont_frac_cocycle import projection_matrix
        sage: projection_matrix(3,2)
        [-0.866025403784439  0.866025403784439  0.000000000000000]
        [-0.500000000000000 -0.500000000000000   1.00000000000000]
        sage: projection_matrix(2,3)
        [-0.577350269189626 -0.333333333333333]
        [ 0.577350269189626 -0.333333333333333]
        [ 0.000000000000000  0.666666666666667]
    """
    from math import sqrt
    from sage.rings.real_mpfr import RR
    sqrt3 = sqrt(3)
    if dim_from == 3 and dim_to == 2:
        return matrix(2,[-sqrt3,sqrt3,0,-1,-1,2],ring=RR)/2
    elif dim_from == 2 and dim_to == 3:
        return matrix(3,[-sqrt3,-1,sqrt3,-1,0,2],ring=RR)/3
    else:
        s = "for input dim_from={} and dim_to={}"
        raise NotImplementedError(s.format(dim_from, dim_to))

def distorsion(M, p=1):
    r"""
    1 Avril 2014. L'ancien ratio n'était pas le bon. Je n'utilisais pas les
    bonnes normes.

    EXAMPLES::

        sage: from slabbe.mult_cont_frac_cocycle import distorsion
        sage: M = matrix(3, (1,2,3,4,5,6,7,8,9))
        sage: M
        [1 2 3]
        [4 5 6]
        [7 8 9]
        sage: distorsion(M)
        3/2
        sage: (3+6+9) / (1+4+7)
        3/2
        sage: distorsion(M, p=oo)
        9/7
    """
    norms = [c.norm(p=p) for c in M.columns()]
    return max(norms) / min(norms)

def is_pisot(m):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_cocycle import is_pisot
        sage: M = matrix(3, (1,2,3,4,5,6,7,8,9))
        sage: is_pisot(M)
        False
    """
    S = sorted((abs(e) for e in m.eigenvalues()), reverse=True)
    return S[0] > 1 and S[1] < 1

def perron_right_eigenvector(M):
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac_cocycle import perron_right_eigenvector
        sage: m = matrix(2,[-11,14,-26,29])
        sage: perron_right_eigenvector(m)           # tolerance 0.00001
        (15.0000000000000, (0.35, 0.6499999999999999))
    """
    from sage.modules.free_module_element import vector
    from sage.rings.real_mpfr import RR
    from sage.rings.all import CC
    import numpy
    eig, vec = numpy.linalg.eig(M)
    index = abs(eig).argmax()
    rightv = vec.transpose()[index]
    if eig[index].imag == 0:
        eig_sage = RR(eig[index].real)
        vec_sage = vector(a.real for a in rightv)
    else:
        eig_sage = CC(eig[index])
        vec_sage = vector(CC, rightv)
    return eig_sage, vec_sage/sum(vec_sage)

def semi_norm_v(M, v,  p=2, verbose=False):
    r"""
    Return the semi norm on the hyperplane orthogonal to v.

    EXAMPLES::

        sage: from slabbe.mult_cont_frac_cocycle import semi_norm_v
        sage: A1 = matrix(3, [1,-1,-1, 0,1,0, 0,0,1]).inverse()
        sage: semi_norm_v(A1, vector( (1,1,1)))      # tolerance 0.0001
        0.9999999999890247
        sage: semi_norm_v(A1, vector( (1,1,1)), p=1)   # tolerance 0.0001
        0.9999394820959548
        sage: semi_norm_v(A1, vector( (1,1,1)), p=oo)   # tolerance 0.0001
        1.0

    """
    from sage.modules.free_module_element import vector
    from sage.numerical.optimize import minimize_constrained
    def func(z):
        vz = vector(z)
        return - (M*vz).norm(p) / vz.norm(p)
    cons = [lambda z: v * vector(z),
            lambda z: - v * vector(z)]
    x0 = range(len(v))
    x0[0] = v[1]
    x0[1] = -v[0]
    rep = minimize_constrained(func, cons, x0)
    if verbose:
        print rep, rep.norm(), rep*v
    return -func(rep)

def semi_norm_cone(M, cone,  p=2, verbose=False):
    r"""
    Return the semi norm on the hyperplane orthogonal to v where v lives in
    the cone.

    EXAMPLES:

    For Arnoux-Rauzy, only the 1-norm works::

        sage: from slabbe.mult_cont_frac_cocycle import semi_norm_cone
        sage: A1 = matrix(3, [1,1,1, 0,1,0, 0,0,1])
        sage: cone = A1
        sage: semi_norm_cone(A1.transpose(), cone, p=1)    # tolerance 0.00001
        0.9999999999999998
        sage: semi_norm_cone(A1.transpose(), cone, p=oo)   # tolerance 0.0001
        1.9999757223144654
        sage: semi_norm_cone(A1.transpose(), cone, p=2)   # tolerance 0.00001
        1.3065629648763757

    For Poincaré, all norms work::

        sage: P21 = matrix(3, [1,1,1, 0,1,1, 0,0,1])
        sage: H21 = matrix(3, [1,0,0, 0,1,0, 1,0,1])
        sage: cone = P21 * H21
        sage: semi_norm_cone(P21.transpose(), cone, p=1)   # tolerance 0.00001
        0.9999957276014074
        sage: semi_norm_cone(P21.transpose(), cone, p=oo)   # tolerance 0.00001
        1.0
        sage: semi_norm_cone(P21.transpose(), cone, p=2)   # tolerance 0.00001
        0.9999999999670175

    For Poincaré on the whole cone, it works for some norms::

        sage: P21 = matrix(3, [1,1,1, 0,1,1, 0,0,1])
        sage: cone = P21
        sage: semi_norm_cone(P21.transpose(), cone, p=1)   # tolerance 0.0001
        1.9999675644077723
        sage: semi_norm_cone(P21.transpose(), cone, p=2)   # tolerance 0.00001
        1.6180339887021953
        sage: semi_norm_cone(P21.transpose(), cone, p=oo)   # tolerance 0.00001
        1.0

    For a product, all norms work::

        sage: A1 = matrix(3, [1,1,1, 0,1,0, 0,0,1])
        sage: P21 = matrix(3, [1,1,1, 0,1,1, 0,0,1])
        sage: H21 = matrix(3, [1,0,0, 0,1,0, 1,0,1])
        sage: M = A1 * P21
        sage: cone = A1 * P21 * H21
        sage: semi_norm_cone(M.transpose(), cone, p=1)   # tolerance 0.00001
        0.999993244882415
        sage: semi_norm_cone(M.transpose(), cone, p=oo)   # tolerance 0.00001
        0.9999935206958908
        sage: semi_norm_cone(M.transpose(), cone, p=2)   # tolerance 0.00001
        0.7529377601317161
    """
    from sage.modules.free_module_element import vector
    from sage.numerical.optimize import minimize_constrained
    a,b,c = cone.columns()
    ab = vector(matrix((a,b)).right_kernel_matrix().row(0))
    ac = vector(matrix((a,c)).right_kernel_matrix().row(0))
    bc = vector(matrix((b,c)).right_kernel_matrix().row(0))
    middle = a+b+c
    cons = []
    if ab * middle < 0:
        cons.append(lambda z: -ab * vector(z))
    else:
        cons.append(lambda z: ab * vector(z))
    if ac * middle < 0:
        cons.append(lambda z: -ac * vector(z))
    else:
        cons.append(lambda z: ac * vector(z))
    if bc * middle < 0:
        cons.append(lambda z: -bc * vector(z))
    else:
        cons.append(lambda z: bc * vector(z))
    if not all(con(middle) > 0 for con in cons):
        raise ValueError("the middle should be in the cone")
    func = lambda v : - semi_norm_v(M,vector(v),p)
    x0 = middle
    rep = minimize_constrained(func, cons, x0)
    if not all((con(rep) >= 0 or abs(con(rep)) < 1e-7) for con in cons):
        raise ValueError("the answer (={}) should be in the cone".format(rep))
    if not all(r >= 0 or abs(r) < 1e-7 for r in rep):
        raise ValueError("the answer (={}) should be positive".format(rep))
    if verbose:
        print "optimal found at ", rep / rep.norm(p)
    return -func(rep)

####################
# Polyhedron Partition
####################
def arp_polyhedron(d=3):
    r"""
    Return the d-dimensional 1-cylinders of the ARP algorithm.

    EXAMPLES::

        sage: from slabbe.mult_cont_frac_cocycle import arp_polyhedron
        sage: A,P,L = arp_polyhedron(3)
        sage: A.vertices_list()
        [[0, 0, 0], [1/2, 1/2, 0], [1/2, 1/4, 1/4], [1, 0, 0]]
        sage: P.vertices_list()
        [[0, 0, 0], [1/2, 1/2, 0], [1/2, 1/4, 1/4], [1/3, 1/3, 1/3]]

    ::

        sage: A,P,L = arp_polyhedron(4)
        sage: A.vertices_list()
        [[0, 0, 0, 0],
         [1/2, 1/2, 0, 0],
         [1/2, 1/6, 1/6, 1/6],
         [1/2, 1/4, 1/4, 0],
         [1, 0, 0, 0]]
        sage: P.vertices_list()
        [[0, 0, 0, 0],
         [1/2, 1/2, 0, 0],
         [1/2, 1/4, 1/4, 0],
         [1/2, 1/6, 1/6, 1/6],
         [1/4, 1/4, 1/4, 1/4],
         [1/3, 1/3, 1/3, 0]]

    ::

        sage: A,P,L = arp_polyhedron(5)
        sage: A.vertices_list()
        [[0, 0, 0, 0, 0],
         [1/2, 1/2, 0, 0, 0],
         [1/2, 1/8, 1/8, 1/8, 1/8],
         [1/2, 1/6, 1/6, 1/6, 0],
         [1/2, 1/4, 1/4, 0, 0],
         [1, 0, 0, 0, 0]]
        sage: P.vertices_list()
        [[0, 0, 0, 0, 0],
         [1/2, 1/2, 0, 0, 0],
         [1/2, 1/6, 1/6, 1/6, 0],
         [1/2, 1/8, 1/8, 1/8, 1/8],
         [1/2, 1/4, 1/4, 0, 0],
         [1/3, 1/3, 1/3, 0, 0],
         [1/5, 1/5, 1/5, 1/5, 1/5],
         [1/4, 1/4, 1/4, 1/4, 0]]
    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    positive = [ [0]*i + [1] + [0]*(d-i) for i in range(1, d+1)]
    atmostone = [[1] + [-1]*d]
    ieq_arnoux = [[0]+[1]+[-1]*(d-1)]
    ieq_arnoux_not = [[0]+[-1]+[1]*(d-1)]
    ieq_sorted = [ [0]*i + [1,-1] + [0]*(d-i-1) for i in range(1,d)]
    A = Polyhedron(ieqs=positive + atmostone + ieq_sorted + ieq_arnoux)
    P = Polyhedron(ieqs=positive + atmostone + ieq_sorted + ieq_arnoux_not)
    L = Polyhedron(ieqs=positive + atmostone + ieq_sorted)
    return A,P,L

def cassaigne_polyhedron(d=3):
    r"""
    Return the d-dimensional 1-cylinders of the Cassaigne algorithm.

    (of the dual!)

    EXAMPLES::

        sage: from slabbe.mult_cont_frac_cocycle import cassaigne_polyhedron
        sage: L,La,Lb = cassaigne_polyhedron(3)
        sage: L.vertices_list()
        [[0, 0, 0], [0, 1/2, 1/2], [1/3, 1/3, 1/3], [1/2, 1/2, 0]]
        sage: La.vertices_list()
        [[0, 0, 0], [0, 1/2, 1/2], [1/3, 1/3, 1/3], [1/4, 1/2, 1/4]]
        sage: Lb.vertices_list()
        [[0, 0, 0], [1/3, 1/3, 1/3], [1/2, 1/2, 0], [1/4, 1/2, 1/4]]

    ::

        sage: L,La,Lb = cassaigne_polyhedron(4)
        sage: L.vertices_list()
        [[0, 0, 0, 0],
         [0, 1/3, 1/3, 1/3],
         [1/3, 1/3, 1/3, 0],
         [1/4, 1/4, 1/4, 1/4],
         [1/5, 2/5, 1/5, 1/5],
         [1/5, 1/5, 2/5, 1/5]]

    ::

        sage: L,La,Lb = cassaigne_polyhedron(5)
        sage: L.vertices_list()
        [[0, 0, 0, 0, 0],
         [0, 1/4, 1/4, 1/4, 1/4],
         [1/4, 1/4, 1/4, 1/4, 0],
         [1/6, 1/6, 1/3, 1/6, 1/6],
         [1/5, 1/5, 1/5, 1/5, 1/5],
         [1/6, 1/3, 1/6, 1/6, 1/6],
         [1/7, 2/7, 2/7, 1/7, 1/7],
         [1/7, 2/7, 1/7, 2/7, 1/7],
         [1/7, 1/7, 2/7, 2/7, 1/7],
         [1/6, 1/6, 1/6, 1/3, 1/6]]
    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    # [-1,7,3,4] represents the inequality 7x_1+3x_2+4x_3>= 1.
    positive = [ [0]*i + [1] + [0]*(d-i) for i in range(1, d+1)]
    atmostone = [[1] + [-1]*d]
    ai_lt_a1d = [[0]+[1]+[0]*(i-2)+[-1]+[0]*(d-i-1)+[1] for i in range(2,d)]
    ai_gt_a1 = [[0]+[-1]+[0]*(i-2)+[1]+[0]*(d-i-1)+[0] for i in range(2,d)]
    ai_gt_ad = [[0]+[0]+[0]*(i-2)+[1]+[0]*(d-i-1)+[-1] for i in range(2,d)]
    a1_gt_ad = [[0]+[1]+[0]*(d-2)+[-1]]
    a1_lt_ad = [[0]+[-1]+[0]*(d-2)+[1]]
    L = Polyhedron(ieqs=positive + atmostone + ai_lt_a1d + ai_gt_a1 + ai_gt_ad)
    La = Polyhedron(ieqs=positive+atmostone+ai_lt_a1d+ai_gt_a1+ai_gt_ad+a1_lt_ad)
    Lb = Polyhedron(ieqs=positive+atmostone+ai_lt_a1d+ai_gt_a1+ai_gt_ad+a1_gt_ad)
    return L, La, Lb

