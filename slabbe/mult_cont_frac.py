# -*- coding: utf-8 -*-
r"""
Multidimensional Continued Fraction Algorithms (Python code)

EXAMPLES::

    sage: from slabbe.mult_cont_frac import Brun
    sage: algo = Brun()

Drawing the natural extension::

    sage: fig = algo.natural_extension_plot(3000, norm_xyz=1, axis_off=True)
    sage: fig
    <matplotlib.figure.Figure object at ...>
    sage: fig.savefig('a.png')  # not tested

Drawing the invariant measure::

    sage: fig = algo.invariant_measure_wireframe_plot(10^6, 50)
    sage: fig
    <matplotlib.figure.Figure object at ...>
    sage: fig.savefig('a.png')  # not tested

Word with given frequencies::

    sage: algo.s_adic_word((1,e,pi))
    word: 1232323123233231232332312323123232312323...

Construction of the same s-adic word from the substitutions and the coding
iterator::

    sage: from itertools import repeat
    sage: D = algo.substitutions()
    sage: it = algo.coding_iterator((1,e,pi))
    sage: words.s_adic(it, repeat(1), D)
    word: 1232323123233231232332312323123232312323...

AUTHORS:

 - Sébastien Labbé, Externalize Python only functions (pip install takes now 33s
   instead of 51s), August 2016
"""
#*****************************************************************************
#       Copyright (C) 2013-2016 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function

from sage.structure.dynamic_class import dynamic_class
from slabbe import mult_cont_frac_pyx

PGF_COLORS = ["red", "green", "blue", "cyan", "brown", "gray", "orange", "pink",
"yellow", "black", "white", "darkgray", "lightgray",
"lime", "olive", "magenta", "purple", "teal", "violet"]

####################
# MCF Python methods
####################
class _MCFAlgorithm_methods(object):
    ######################
    # COMBINATORICS METHODS
    ######################
    def matrix_cocycle(self):
        r"""
        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ARP
            sage: ARP().matrix_cocycle()
            Cocycle with 9 gens over Regular language over [1, 2, 3, 
            123, 132, 213, 231, 312, 321]
            defined by: Automaton with 7 states

        ::

            sage: from slabbe.mult_cont_frac import Sorted_Brun
            sage: Sorted_Brun().matrix_cocycle()
            Cocycle with 3 gens over Language of finite words over alphabet
            [1, 2, 3]

        ::

            sage: from slabbe.mult_cont_frac import Brun
            sage: Brun().matrix_cocycle()
            Cocycle with 6 gens over Regular language over 
            [123, 132, 213, 231, 312, 321]
            defined by: Automaton with 6 states
        """
        from .matrix_cocycle import cocycles
        try:
            f = getattr(cocycles, self.class_name())
        except AttributeError:
            msg = "Matrix cocyle not implemented for {}"
            msg = msg.format(self.class_name())
            raise NotImplementedError(msg)
        return f()

    def n_matrix(self, start, n_iterations):
        r"""
        Return the n-matrix associated to the direction v.

        INPUT:

        - ``start`` -- iterable of three real numbers
        - ``n_iterations`` - integer, number of iterations

        OUTPUT:

            matrix

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ARP
            sage: ARP().n_matrix((1,e,pi), 10)
            [ 31  40   7]
            [ 84 109  19]
            [ 97 126  22]
        """
        it = self.coding_iterator(start)
        L = [next(it) for _ in range(n_iterations)]
        cocycle = self.matrix_cocycle()
        return cocycle.word_to_matrix(L)

    def s_adic_word(self, v=None, n_iterations=100, nth_letter=1):
        r"""
        Return the s-adic word obtained from application of the MCF
        algorithm on the vector v.

        INPUT:

        - ``v`` - initial vector (default: ``None``), if None, then
          initial point is random
        - ``n_iterations`` - integer (default: ``100``), number of
          iterations
        - ``nth_letter`` - letter (default: ``1``)

        OUTPUT:

            word

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun, ARP, Cassaigne
            sage: Brun().s_adic_word((.414578,.571324,.65513))
            word: 1232312312323123123312323123123312323123...
            sage: Brun().s_adic_word((1,e,pi))
            word: 1232323123233231232332312323123232312323...
            sage: ARP().s_adic_word((1,e,pi))
            word: 1232323123233231232332312323123232312323...
            sage: Cassaigne().s_adic_word((1,e,pi))
            word: 2323213232323132323213232321323231323232...

        On integer entries::

            sage: w = Brun().s_adic_word((3,4,5))
            sage: w
            word: 123233123123
            sage: w.abelian_vector()
            [3, 4, 5]

        When the gcd is not one::

            sage: Brun().s_adic_word((2,4,10))
            word: 1233323312333233

        ::

            sage: from slabbe.mult_cont_frac import FullySubtractive
            sage: FullySubtractive().s_adic_word((1,2,5))
            Traceback (most recent call last):
            ...
            ValueError: On input=(1, 2, 5), algorithm Fully Subtractive
            loops on (1.0, 0.0, 3.0)

        ::

            sage: from slabbe.mult_cont_frac import Reverse
            sage: algo = Reverse()
            sage: algo.s_adic_word((18,1,1))
            word: 31111111111211111111
            sage: _.abelian_vector()
            [18, 1, 1]

        ::

            sage: Reverse().s_adic_word((3,1,1))
            Traceback (most recent call last):
            ...
            ValueError: On input=(3, 1, 1), algorithm Reverse reaches non
            integer entries (0.5, 0.5, 0.5)

        TESTS::

            sage: v = ARP().s_adic_word((1,e,pi))
            sage: w = Brun().s_adic_word((1,e,pi))
            sage: v.longest_common_prefix(w, 'finite').length()
            212
        """
        from sage.combinat.words.word_generators import words
        from sage.rings.integer_ring import ZZ
        if v is None:
            v = (random(), random(), random())
        D = self.substitutions()
        if all(a in ZZ for a in v):
            S = []
            it = self.cone_orbit_iterator(v)
            previousA = None
            for P,b in it:
                A,B = P.to_tuple()
                if not all(a in ZZ for a in A):
                    raise ValueError("On input={}, algorithm {} reaches"
                            " non integer entries {}".format(v, self.name(), A))
                if A == previousA:
                    break
                S.append(b)
                previousA = A
            (x,y,z) = A
            if x == 0 == y: 
                letter = 3
                the_gcd = ZZ(z)
            elif x == 0 == z:
                letter = 2
                the_gcd = ZZ(y)
            elif y == 0 == z:
                letter = 1
                the_gcd = ZZ(x)
            else:
                raise ValueError("On input={}, algorithm {} loops"
                                 " on {}".format(v, self.name(), A))
            return words.s_adic(S, [letter], D)**the_gcd

        else:
            it = self.coding_iterator(v)
            S = [next(it) for _ in range(n_iterations)]
            letter = nth_letter
            L = [letter]
            for key in reversed(S):
                letter = D[key](letter)[0]
                L.append(letter)
            L.pop()
            L.reverse()
            return words.s_adic(S, L, D)

    def discrepancy_statistics(self, length):
        r"""
        Return the discrepancy of words of given length.

        INPUT:

        - ``length`` -- integer

        OUTPUT:

            dict

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
            sage: Brun().discrepancy_statistics(5)
            {[1, 1, 3]: 6/5,
             [1, 2, 2]: 4/5,
             [1, 3, 1]: 4/5,
             [2, 1, 2]: 4/5,
             [2, 2, 1]: 4/5,
             [3, 1, 1]: 4/5}
        """
        from .finite_word import discrepancy
        from sage.combinat.composition import Compositions
        D = {}
        for c in Compositions(length, length=3, min_part=1):
            w = self.s_adic_word(c)
            if c != w.abelian_vector(): 
                raise ValueError("c={} but vector is"
                      " {}".format(c,w.abelian_vector()))
            D[c] = discrepancy(w)
        return D

    def e_one_star_patch(self, v, n):
        r"""
        Return the n-th iterated patch of normal vector v.

        INPUT:

        - ``v`` -- vector, the normal vector
        - ``n`` -- integer

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import FullySubtractive
            sage: FullySubtractive().e_one_star_patch((1,e,pi), 4)
            Patch of 21 faces
        """
        from sage.combinat.e_one_star import E1Star, Patch, Face
        from sage.misc.misc_c import prod
        if v is None:
            v = (random(), random(), random())
        it = self.coding_iterator(v)
        keys = [next(it) for _ in range(n)]
        D = self.dual_substitutions()
        L = prod(D[key] for key in reversed(keys))
        dual_sub = E1Star(L)
        cube = Patch([Face((1,0,0),1), Face((0,1,0),2), Face((0,0,1),3)])
        return dual_sub(cube)

    def _base_translation_vectors(self):
        r"""
        Return the base translation vectors associated to each elementary matrix

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import FullySubtractive
            sage: algoFS = FullySubtractive()
            sage: algoFS._base_translation_vectors()
            {1: (1, 0, 0), 2: (0, 1, 0), 3: (0, 0, 1)}

        ::

            sage: from slabbe.mult_cont_frac import ARP
            sage: algoARP = ARP()
            sage: algoARP._base_translation_vectors()
            {1: (0, 1, 1),
             2: (1, 0, 1),
             3: (1, 1, 0),
             123: (1, 1, 0),
             132: (1, 0, 1),
             213: (1, 1, 0),
             231: (0, 1, 1),
             312: (1, 0, 1),
             321: (0, 1, 1)}
        """
        from sage.modules.free_module_element import vector
        def e(i):
            A = [0,0,0]
            A[i-1] = 1
            return vector(A)
        D = {}
        gens = self.matrix_cocycle().gens()
        for key,M in gens.iteritems():
            Minv = M.inverse()
            v = sum(e(list(c).index(1)+1) for c in Minv.columns() if -1 in c)
            D[key] = v
        return D

    def translation_vectors(self, start, n_iterations):
        r"""
        Return the base translation vectors associated to each elementary matrix

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import FullySubtractive, ARP, Brun
            sage: ARP().translation_vectors((1,e,pi), 5)
            [(1, 1, 0), (1, -1, 1), (-2, 1, 0), (1, 1, -1), (-3, 0, 1)]
            sage: Brun().translation_vectors((1,e,pi), 5)
            [(0, 1, 0), (1, 0, 0), (1, 0, 0), (-2, 1, 0), (0, -1, 1)]
            sage: FullySubtractive().translation_vectors((1,e,pi), 5)
            [(1, 0, 0), (1, 0, 0), (-2, 1, 0), (3, -1, 0), (-3, 0, 1)]
        """
        D = self._base_translation_vectors()
        it = self.coding_iterator(start)
        L = [next(it) for _ in range(n_iterations)]
        cocycle = self.matrix_cocycle()
        gens_inv = cocycle.gens_inverses()
        Minv = cocycle.identity_matrix()
        T = []
        for a in L:
            v = D[a] * Minv
            v.set_immutable()
            T.append(v)
            Minv = gens_inv[a] * Minv
        return T

    def discrete_plane_patches(self, start, n_iterations):
        r"""
        Return the patches subsets of a discrete plane using method from
        [JLP2016]_.

        INPUT:

        - ``start`` -- normal vector
        - ``n_iterations`` -- integer, number of iterations

        EXAMPLES:

        Partial sums for Brun forms a connected subset of Z^3::

            sage: from slabbe.mult_cont_frac import FullySubtractive, ARP, Brun
            sage: A = Brun().discrete_plane_patches((1,e,pi), 10)
            sage: _ = A.tikz().pdf()     # long time

        ::

            sage: A = FullySubtractive().discrete_plane_patches((1,e,pi), 10)
            sage: _ = A.tikz().pdf()     # not tested

        Result is not connected::

            sage: A = ARP().discrete_plane_patches((1,e,pi), 10)
            sage: _ = A.tikz().pdf()     # not tested

        REFERENCES:

        .. [JLP2016] Jamet, Damien, Nadia Lafrenière, et Xavier Provençal. « Generation
           of Digital Planes Using Generalized Continued-Fractions Algorithms
           ». In Discrete Geometry for Computer Imagery, édité par Nicolas
           Normand, Jeanpierre Guédon, et Florent Autrusseau, 45-56. Lecture
           Notes in Computer Science 9647. Springer International Publishing,
           2016. http://link.springer.com/chapter/10.1007/978-3-319-32360-2_4.
        """
        T = self.translation_vectors(start, n_iterations)
        zero = T[0].parent().zero()
        from sage.combinat.subset import Subsets
        S = Subsets(range(len(T)))
        V = [sum((T[s] for s in subset), zero) for subset in S]
        from slabbe import DiscreteSubset
        return DiscreteSubset.from_subset(V)

    ######################
    # DRAWINGS METHODS (python):
    ######################
    def invariant_measure_wireframe_plot(self, n_iterations, ndivs, norm=1):
        r"""
        Return a matplotlib graph of the invariant measure.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``ndvis`` - integer, number of divisions per dimension
        - ``norm`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the orbit points

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Reverse, Brun
            sage: Reverse().invariant_measure_wireframe_plot(1000000, 80)
            <matplotlib.figure.Figure object at ...>
            sage: Brun().invariant_measure_wireframe_plot(1000000, 40, norm=1)
            <matplotlib.figure.Figure object at ...>

        """
        D = self._invariant_measure_dict(n_iterations, ndivs, norm=norm)
        the_mean = n_iterations / float(len(D))

        X = [[i for i in range(ndivs+1)] for j in range(ndivs+1)]
        Y = [[j for i in range(ndivs+1)] for j in range(ndivs+1)]
        Z = [[D.get((i,j),0)/the_mean for i in range(ndivs+1)] for j in range(ndivs+1)]
        
        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)

        ax.text(ndivs, 0, 0, "$(1,0,0)$", color='black', va='top', size=20)
        ax.text(0, ndivs, 0, "$(0,1,0)$", color='black', ha='left', size=20)
        ax.text(0, 0, 0, "$(0,0,1)$", color='black', ha='right', size=20)

        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zlabel('Frequency')

        title = "Density of an orbit for\n{} algorithm".format(self.name())
        ax.set_title(title)

        return fig

    def invariant_measure_contour_plot(self, n_iterations, ndivs, norm=1):
        r"""
        Return a matplotlib graph of the invariant measure.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``ndvis`` - integer, number of divisions per dimension
        - ``norm`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the orbit points

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Reverse, Brun
            sage: Reverse().invariant_measure_contour_plot(1000000, 80)
            <matplotlib.figure.Figure object at ...>
            sage: Brun().invariant_measure_contour_plot(1000000, 40, norm=1)
            <matplotlib.figure.Figure object at ...>

        """
        D = self._invariant_measure_dict(n_iterations, ndivs, norm=norm)
        the_mean = n_iterations / float(len(D))
        S = sorted(D.values())
        V = [S[k]/the_mean for k in range(0, len(D), len(D)/10)]

        X = [[i for i in range(ndivs+1)] for j in range(ndivs+1)]
        Y = [[j for i in range(ndivs+1)] for j in range(ndivs+1)]
        Z = [[D.get((i,j),0)/the_mean for i in range(ndivs+1)] for j in range(ndivs+1)]
        
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        CS = plt.contour(X, Y, Z, V)
        plt.clabel(CS, inline=1, fontsize=10)

        title = "Density of an orbit for\n{} algorithm".format(self.name())
        ax.set_title(title)

        return fig

    def natural_extension_plot(self, n_iterations, norm_xyz=1,
            norm_uvw=1, axis_off=False):
        r"""
        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``norm_xyz`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the orbit points
        - ``norm_uvw`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the dual orbit points
        - ``axis_off`` - boolean (default: ``False``), 

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
            sage: Brun().natural_extension_plot(3000, norm_xyz=1, axis_off=True)
            <matplotlib.figure.Figure object at ...>

        ::

            sage: from slabbe.mult_cont_frac import Sorted_ARP
            sage: Sorted_ARP().natural_extension_plot(1000)
            <matplotlib.figure.Figure object at ...>
        """
        import matplotlib
        import matplotlib.pyplot as plt
        t = self._natural_extension_dict(n_iterations, norm_xyz=norm_xyz,
                norm_uvw=norm_uvw)
        domain_left, image_left, domain_right, image_right = t
        c = dict(zip(domain_left.keys(), ['b','r','g','c','m','y','k']))

        def create_sub_plot(axx, D, title, norm):
            for key, value in D.iteritems():
                value = D[key]
                X,Y = zip(*value)
                A = axx.plot(X, Y, 'o', markersize=2, color=c[key], label=key)
            axx.legend(markerscale=3,fontsize="xx-small")
            axx.set_title(title)

        fig, (ax,bx,cx,dx) = plt.subplots(1,4,sharex=True,sharey=True)
        fig.set_figheight(2)
        fig.set_figwidth(12)
        if axis_off:
            ax.set_axis_off()
            bx.set_axis_off()
            cx.set_axis_off()
            dx.set_axis_off()

        create_sub_plot(ax, domain_left,  "Algo IN",    norm=norm_xyz)
        create_sub_plot(bx, domain_right,    "NatExt IN",  norm=norm_uvw)
        create_sub_plot(cx, image_left, "Algo OUT",   norm=norm_xyz)
        create_sub_plot(dx, image_right,   "NatExt OUT", norm=norm_uvw)

        title = "Algo=%s, nbiter=%s" % (self.name(), n_iterations)
        fig.suptitle(title)

        return fig

    def natural_extension_tikz(self, n_iterations, norm_xyz=1,
            norm_uvw=1, marksize=0.2, legend_marksize=2, 
            group_size="4 by 1"):
        r"""

        INPUT:

        - ``n_iterations`` -- integer, number of iterations
        - ``norm_xyz`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the orbit points
        - ``norm_uvw`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the dual orbit points
        - ``marksize`` -- tikz marksize (default:``0.2``)
        - ``legend_marksize`` -- tikz legend marksize (default:``2``)
        - ``group_size`` -- string (default:``"4 by 1"``)

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import Brun
            sage: s = Brun().natural_extension_tikz(1000)
            sage: s
            \documentclass[tikz]{standalone}
            \usepackage{pgfplots}
            \usetikzlibrary{pgfplots.groupplots}
            \begin{document}
            \begin{tikzpicture}[scale=.7]
            \begin{groupplot}
            [group style={group size=4 by 1},
            height=7cm,width=8cm,
            xmin=-1.1,xmax=1.1,ymin=-.6,ymax=1.20,
            ...
            ... 4184 lines not printed (17... characters in total) ...
            ...
            \draw[draw=none] (group c2r1.center) --
            node {$\to$}     (group c3r1.center);
            \draw[draw=none] (group c3r1.center) --
            node {$\times$}  (group c4r1.center);
            \end{tikzpicture}
            \end{document}

        ::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.pdf')
            sage: _ = s.pdf(filename)
        """
        t = self._natural_extension_dict(n_iterations, norm_xyz=norm_xyz,
                norm_uvw=norm_uvw)
        domain_left, image_left, domain_right, image_right = t
        sqrt3 = 1.73205080756888
        r = 1.08
        color_dict = dict(zip(domain_left.keys(), PGF_COLORS))

        lines = []
        lines.append(r"\begin{tikzpicture}[scale=.7]")
        lines.append(r"\begin{groupplot}")
        lines.append(r"[group style={group size=%s}," % group_size)
        lines.append(r"height=7cm,width=8cm,")
        lines.append(r"xmin=-1.1,xmax=1.1,ymin=-.6,ymax=1.20,")
        lines.append(r"hide axis]")
        datas = [domain_left, domain_right, image_left, image_right]
        labels = [r"$\mathbf{x}_n$",r"$\mathbf{a}_n$",
                  r"$\mathbf{x}_{n+1}$",r"$\mathbf{a}_{n+1}$"]
        for data,label in zip(datas, labels):
            lines.append(r"\nextgroupplot")
            lines.append(r"\draw[dashed] ")
            lines.append(r"(axis cs:%s, %s)" % (-r*sqrt3/2,r*-.5))
            lines.append(r" node[left] {$\mathbf{e}_1$} -- ")
            lines.append("(axis cs:%s, %s)" % (r*sqrt3/2,r*-.5))
            lines.append(r" node[right] {$\mathbf{e}_2$} -- ")
            lines.append("(axis cs:%s, %s)" % (0, r))
            lines.append(r" node[above] {$\mathbf{e}_3$} -- cycle;")
            lines.append(r"\node at (axis cs:%s, %s) {%s};" % (-.5,.8,label))
            for key,value in data.iteritems():
                lines.append(r"\addplot+[")
                lines.append(r"legend image post style={mark size=%s}," % legend_marksize)
                lines.append(r"only marks,mark=*,")
                lines.append(r"mark size=%s," % marksize)
                lines.append(r"mark options={color=%s}] " % color_dict[key])
                lines.append(r"coordinates {%s};" % '\n'.join(map(str, value)))
                lines.append(r"\addlegendentry{%s}" % key)
        lines.append(r"\end{groupplot}")
        if group_size == "4 by 1":
            lines.append(r"\draw[draw=none] (group c1r1.center) -- ")
            lines.append(r"node {$\times$}  (group c2r1.center);")
            lines.append(r"\draw[draw=none] (group c2r1.center) -- ")
            lines.append(r"node {$\to$}     (group c3r1.center);")
            lines.append(r"\draw[draw=none] (group c3r1.center) -- ")
            lines.append(r"node {$\times$}  (group c4r1.center);")
        lines.append(r"\end{tikzpicture}")
        from slabbe import TikzPicture
        return TikzPicture('\n'.join(lines), usepackage=['pgfplots'],
                usetikzlibrary=['pgfplots.groupplots'])

    def natural_extension_part_tikz(self, n_iterations, part=3, 
                                    norm_xyz=1, norm_uvw=1,
                                    marksize='1pt', limit_nb_points=None,
                                    verbose=False):
        r"""
        Return a pgfplots or some part of an orbit in the natural
        extension.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``part`` - integer, taking value 0, 1, 2 or 3
        - ``norm_xyz`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the orbit points
        - ``norm_uvw`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the dual orbit points
        - ``marksize`` -- string (default: ``'1pt'``), pgfplots mark size value
        - ``limit_nb_points`` -- None or integer (default: ``None``), limit
          number of points per patch
        - ``verbose`` -- string (default: ``False``)

        EXAMPLES::

            sage: from slabbe.mult_cont_frac import ARP
            sage: s = ARP().natural_extension_part_tikz(1000, part=3)
            sage: view(s, tightpage=True)    # not tested

        also ::

            sage: s = ARP().natural_extension_part_tikz(1000, part=3, marksize='.2pt', limit_nb_points=1200, verbose=True)
            Taking ... points for key ...
            Taking ... points for key ...
            Taking ... points for key ...
            Taking ... points for key ...
            Taking ... points for key ...
            Taking ... points for key ...
            Taking ... points for key ...
            Taking ... points for key ...
            Taking ... points for key ...
        """
        t = self._natural_extension_dict(n_iterations, norm_xyz=norm_xyz,
                norm_uvw=norm_uvw)
        domain_left, image_left, domain_right, image_right = t
        sqrt3 = 1.73205080756888
        r = 1.08
        color_dict = dict(zip(domain_left.keys(), PGF_COLORS))
        data = t[part]

        s = ''
        s += "\\begin{tikzpicture}[scale=.7]\n"
        s += ("\\begin{axis}[height=7cm,width=8cm,\n"
               "xmin=-1.1,xmax=1.1,ymin=-.6,ymax=1.20,\n"
               "hide axis]\n")
        s += ("\\draw[dashed] \n"
              "(axis cs:%s, %s)" % (-r*sqrt3/2,r*-.5) +
              " node[left] {$\\mathbf{e}_1$} -- \n"
              "(axis cs:%s, %s)" % (r*sqrt3/2,r*-.5)  +
              " node[right] {$\\mathbf{e}_2$} -- \n"
              "(axis cs:%s, %s)" % (0, r)             +
              " node[above] {$\\mathbf{e}_3$} -- cycle;\n")
        for key,value in data.iteritems():
            if limit_nb_points and len(value) > limit_nb_points:
                if verbose:
                    print("Taking only {} points instead of {} for key {}".format(
                            limit_nb_points, len(value), key))
                value = value[:limit_nb_points]
            elif verbose:
                print("Taking {} points for key {}".format(len(value),
                    key))

            s += "\\addplot+[only marks,mark=*,mark options={color=%s}," % color_dict[key]
            s += "mark size=%s]\n" % marksize
            s += "coordinates {\n"
            s += '\n'.join(map(str, value))
            s += "};\n" 
            s += "\\addlegendentry{%s}\n " % key
        s += "\\end{axis}\n"
        s += "\\end{tikzpicture}\n"
        from slabbe import TikzPicture
        return TikzPicture(s, usepackage=['pgfplots'])

    def natural_extension_part_png(self, n_iterations, draw,
                                    norm_xyz=1, norm_uvw=1,
                                    xrange=(-.866, .866),
                                    yrange=(-.5, 1.),
                                    urange=(-.866, .866),
                                    vrange=(-.5, 1.),
                                    color_dict=None,
                                    branch_order=None,
                                    ndivs=1024,
                                    verbose=False):
        r"""
        Return a png or some part of an orbit in the natural extension.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``draw`` -- string (default: ``'image_right'``), possible values
          are:

          - ``'domain_left'`` - use x and y ranges
          - ``'domain_right'`` - use u and v ranges
          - ``'image_left'`` - use x and y ranges
          - ``'image_right'`` - use u and v ranges

        - ``norm_xyz`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the orbit points
        - ``norm_uvw`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the dual orbit points
        - ``xrange`` -- tuple (default: ``(-.866, .866)``), interval of
          values for x
        - ``yrange`` -- tuple (default: ``(-.5, 1.)``), interval of
          values for y
        - ``urange`` -- tuple (default: ``(-.866, .866)``), interval of
          values for u
        - ``vrange`` -- tuple (default: ``(-.5, 1.)``), interval of
          values for v
        - ``color_dict`` -- dict (default: ``None``), dict from branches
          int to color (as RGB tuples)
        - ``branch_order`` -- list (default: ``None``), list of branches int
        - ``ndivs`` -- int (default: ``1024``), number of pixels
        - ``verbose`` -- string (default: ``False``)

        BENCHMARK:

        A minute (1min 13s) for a picture with 10^7 points::

            sage: from slabbe.mult_cont_frac import ARP
            sage: c = {}
            sage: c[1] = c[2] = c[3] = [0,0,0]
            sage: c[12] = c[13] = c[23] = c[21] = c[31] = c[32] = [255,0,0]
            sage: b = [1,2,3,12,13,21,23,31,32]
            sage: P = ARP().natural_extension_part_png(10^7, draw='image_right',  # not tested
            ....:   branch_order=b, color_dict=c, urange=(-.6,.6), vrange=(-.6,.6))  # not tested

        Half a minute (27s) for a picture zoomed in the orbit of 10^8
        points::

            sage: P = ARP().natural_extension_part_png(10^8, draw='image_right', # not tested
            ....:   branch_order=b, color_dict=c, urange=(.2,.3), vrange=(.2,.3))   # not tested

        EXAMPLES::

            sage: c = {}
            sage: c[1] = c[2] = c[3] = [0,0,0]
            sage: c[123] = c[132] = c[231] = c[213] = c[312] = c[321] = [255,0,0]
            sage: b = [1,2,3,123,132,213,231,312,321]
            sage: opt = dict(urange=(-.6,.6), vrange=(-.6,.6), color_dict=c, branch_order=b)
            sage: P = ARP().natural_extension_part_png(10^5, draw='domain_left', **opt)
            sage: P = ARP().natural_extension_part_png(10^5, draw='domain_right', **opt)
            sage: P = ARP().natural_extension_part_png(10^5, draw='image_left', **opt)
            sage: P = ARP().natural_extension_part_png(10^5, draw='image_right', **opt)
            sage: P.show() # not tested

        """
        L = self.simplex_orbit_filtered_list(n_iterations=n_iterations,
            norm_xyz=norm_xyz, norm_uvw=norm_uvw,
            xmin=xrange[0], xmax=xrange[1],
            ymin=yrange[0], ymax=yrange[1],
            umin=urange[0], umax=urange[1],
            vmin=vrange[0], vmax=vrange[1],
            ndivs=ndivs)
        if branch_order is None:
            branch_order = []
            raise NotImplementedError
        if color_dict is None:
            from random import randint
            color_dict = {}
            for key in branch_order:
                color_dict[key] = [randint(0,255),randint(0,255),randint(0,255)]

        #http://stackoverflow.com/questions/434583/what-is-the-fastest-way-to-draw-an-image-from-discrete-pixel-values-in-python
        import numpy as np
        import scipy.misc as smp

        # Create a 1024x1024x3 array of 8 bit unsigned integers
        data = np.zeros( (ndivs,ndivs,3), dtype=np.uint8 )
        data += 255   # white as default color

        if draw.startswith('domain'):
            L.sort(key=lambda a:branch_order.index(a[5]))
        elif draw.startswith('image'):
            L.sort(key=lambda a:branch_order.index(a[4]))
        else:
            raise ValueError("Unkown value for draw(={})".format(draw))

        if draw == 'domain_left':
            for x,y,u,v,prev_br,next_br in L:
                data[y,x] = color_dict[next_br]
        elif draw == 'domain_right':
            for x,y,u,v,prev_br,next_br in L:
                data[v,u] = color_dict[next_br]
        elif draw == 'image_left':
            for x,y,u,v,prev_br,next_br in L:
                data[y,x] = color_dict[prev_br]
        elif draw == 'image_right':
            for x,y,u,v,prev_br,next_br in L:
                data[v,u] = color_dict[prev_br]
        else:
            raise ValueError("Unkown value for draw(={})".format(draw))

        img = smp.toimage( data )       # Create a PIL image
        #img.show()                      # View in default viewer
        return img

    def measure_evaluation(self, n_iterations, draw,
                                norm_xyz=1, norm_uvw=1,
                                xrange=(-.866, .866),
                                yrange=(-.5, 1.),
                                urange=(-.866, .866),
                                vrange=(-.5, 1.),
                                ndivs=1024,
                                verbose=False):
        r"""
        Return the measure of a box according to an orbit.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``draw`` -- string (default: ``'image_right'``), possible values
          are:

          - ``'domain_left'`` - use x and y ranges
          - ``'domain_right'`` - use u and v ranges
          - ``'image_left'`` - use x and y ranges
          - ``'image_right'`` - use u and v ranges

        - ``norm_xyz`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the orbit points
        - ``norm_uvw`` -- string (default: ``1``), either ``0`` or
          ``1``, the p-norm used for the dual orbit points
        - ``xrange`` -- tuple (default: ``(-.866, .866)``), interval of
          values for x
        - ``yrange`` -- tuple (default: ``(-.5, 1.)``), interval of
          values for y
        - ``urange`` -- tuple (default: ``(-.866, .866)``), interval of
          values for u
        - ``vrange`` -- tuple (default: ``(-.5, 1.)``), interval of
          values for v
        - ``ndivs`` -- int (default: ``1024``), number of pixels
        - ``verbose`` -- string (default: ``False``)

        BENCHMARK::

            sage: from slabbe.mult_cont_frac import ARP
            sage: opt = dict(urange=(-.15,.25), vrange=(-.05,.05))
            sage: ARP().measure_evaluation(10^8, draw='right', ndivs=100, **opt) # optional long
            0.435...
            sage: ARP().measure_evaluation(10^8, draw='right', ndivs=1000, **opt) # optional long
            0.357...
            sage: ARP().measure_evaluation(10^8, draw='right', ndivs=2000, **opt) # optional long
            0.293...
            sage: ARP().measure_evaluation(10^8, draw='right', ndivs=4000, **opt) # optional long
            0.177...

        """
        L = self.simplex_orbit_filtered_list(n_iterations=n_iterations,
            norm_xyz=norm_xyz, norm_uvw=norm_uvw,
            xmin=xrange[0], xmax=xrange[1],
            ymin=yrange[0], ymax=yrange[1],
            umin=urange[0], umax=urange[1],
            vmin=vrange[0], vmax=vrange[1],
            ndivs=ndivs)

        if draw.endswith('left'):
            S = set((p[0],p[1]) for p in L)
        elif draw.endswith('right'):
            S = set((p[2],p[3]) for p in L)
        else:
            raise ValueError("Unkown value for draw(={})".format(draw))

        if verbose:
            print("nombre diterations dans la fenetre : ", len(L))
            print("{} pixels touchés parmi limage {}^2 ".format(len(S),
                ndivs))

        print(float(len(S) / ndivs**2))

########################
# Dynamic class creation
########################
def _dynamic_MCF_class(base):
    r"""
    Return the MCF algorithm associated to some Cython base class including all
    the Python methods.

    This is where the fusion of cython methods with python methods is done with
    the use of ``dynamic_class``.

    INPUT:

    - ``base`` -- MCF algorithm base class

    OUTPUT:

        class

    EXAMPLES::

        sage: from slabbe.mult_cont_frac_pyx import Brun
        sage: from slabbe.mult_cont_frac import _dynamic_MCF_class
        sage: cls = _dynamic_MCF_class(Brun)
        sage: algo = cls()
        sage: algo
        Brun 3-dimensional continued fraction algorithm
        sage: algo.matrix_cocycle()
        Cocycle with 6 gens over Regular language over [123, 132, 213, 231, 312, 321]
        defined by: Automaton with 6 states
    """
    class_name = base().class_name()
    return dynamic_class(class_name, (base, _MCFAlgorithm_methods))

################
# MCF Algorithms
################
def Brun():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Brun
        sage: Brun()
        Brun 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Brun)()
def Reverse():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Reverse
        sage: Reverse()
        Reverse 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Reverse)()
def ARP():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import ARP
        sage: ARP()
        Arnoux-Rauzy-Poincar\'e 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.ARP)()
def ArnouxRauzy():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import ArnouxRauzy
        sage: ArnouxRauzy()
        ArnouxRauzy 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.ArnouxRauzy)()
def Poincare():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Poincare
        sage: Poincare()
        Poincar\'e 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Poincare)()
def Selmer():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Selmer
        sage: Selmer()
        Selmer 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Selmer)()
def FullySubtractive():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import FullySubtractive
        sage: FullySubtractive()
        Fully Subtractive 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.FullySubtractive)()
def Cassaigne():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Cassaigne
        sage: Cassaigne()
        Cassaigne 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Cassaigne)()
def Sorted_Brun():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_Brun
        sage: Sorted_Brun()
        Sorted_Brun 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_Brun)()
def Sorted_BrunMulti():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_BrunMulti
        sage: Sorted_BrunMulti()
        Sorted_BrunMulti 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_BrunMulti)()
def Sorted_Selmer():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_Selmer
        sage: Sorted_Selmer()
        Sorted_Selmer 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_Selmer)()
def Sorted_FullySubtractive():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_FullySubtractive
        sage: Sorted_FullySubtractive()
        Sorted_FullySubtractive 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_FullySubtractive)()
def Sorted_ArnouxRauzy():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_ArnouxRauzy
        sage: Sorted_ArnouxRauzy()
        Sorted_ArnouxRauzy 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_ArnouxRauzy)()
def Sorted_ArnouxRauzyMulti():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_ArnouxRauzyMulti
        sage: Sorted_ArnouxRauzyMulti()
        Sorted_ArnouxRauzyMulti 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_ArnouxRauzyMulti)()
def Sorted_ARP():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_ARP
        sage: Sorted_ARP()
        Sorted_ARP 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_ARP)()
def Sorted_ARPMulti():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_ARPMulti
        sage: Sorted_ARPMulti()
        Sorted_ARPMulti 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_ARPMulti)()
def Sorted_Poincare():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_Poincare
        sage: Sorted_Poincare()
        Sorted_Poincare 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_Poincare)()
def Sorted_ARrevert():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_ARrevert
        sage: Sorted_ARrevert()
        Sorted_ARrevert 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_ARrevert)()
def Sorted_ARrevertMulti():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_ARrevertMulti
        sage: Sorted_ARrevertMulti()
        Sorted_ARrevertMulti 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_ARrevertMulti)()
def Sorted_ARMonteil():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_ARMonteil
        sage: Sorted_ARMonteil()
        Sorted_ARMonteil 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_ARMonteil)()
def Sorted_Delaunay():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import Sorted_Delaunay
        sage: Sorted_Delaunay()
        Sorted_Delaunay 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.Sorted_Delaunay)()
def JacobiPerron():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import JacobiPerron
        sage: JacobiPerron()
        JacobiPerron 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.JacobiPerron)()
def JacobiPerronAdditif():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import JacobiPerronAdditif
        sage: JacobiPerronAdditif()
        JacobiPerronAdditif 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.JacobiPerronAdditif)()
def JacobiPerronAdditifv2():
    r"""
    EXAMPLES::

        sage: from slabbe.mult_cont_frac import JacobiPerronAdditifv2
        sage: JacobiPerronAdditifv2()
        JacobiPerronAdditifv2 3-dimensional continued fraction algorithm
    """
    return _dynamic_MCF_class(mult_cont_frac_pyx.JacobiPerronAdditifv2)()
