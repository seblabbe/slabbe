# -*- coding: utf-8 -*-
r"""
Multidimensional Continued Fraction Algorithms

EXAMPLES::

    sage: from slabbe.mcf import algo
    sage: algo.brun.plot_natural_extension(3000, norm_algo='1', axis_off=True)
    Creation du fichier nat_ext_brun_iter3000.png

"""
#*****************************************************************************
#       Copyright (C) 2014 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from slabbe.mcf_pyx import MCFAlgorithm_pyx

from sage.misc.table import table
from sage.misc.latex import latex, LatexExpr
latex.add_to_preamble('\\usepackage{tikz}')
latex.add_to_preamble('\\usepackage{pgfplots}')
latex.add_to_preamble('\\usetikzlibrary{pgfplots.groupplots}')


PGF_COLORS = ["red", "green", "blue", "cyan", "brown", "gray", "orange", "pink",
"yellow", "black", "white", "darkgray", "lightgray",
"lime", "olive", "magenta", "purple", "teal", "violet"]

class MCFAlgorithm(MCFAlgorithm_pyx):
    def plot_invariant_measure(self, n_iterations, ndivs, norm='sup'):
        r"""
        Return a matplotlib graph of the invariant measure.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``ndvis`` - integer, number of divisions per dimension
        - ``norm`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points

        EXAMPLES::

            sage: algo.arrevert.plot_invariant_measure(1000000, 80) # not tested
            Creation du fichier mesure_arrevert_iter1000000_div80.png
            sage: algo.brun.plot_invariant_measure(1000000, 40, norm='sup')
            Creation du fichier mesure_brun_iter1000000_div40.png

        """
        D = self.invariant_measure_dict(n_iterations, ndivs, norm=norm)
        mx,my = map(max, zip(*D.keys()))
        len_D = float(len(D))
        the_mean = n_iterations / len_D

        X = [[i for i in range(mx+1)] for j in range(my+1)]
        Y = [[j for i in range(mx+1)] for j in range(my+1)]
        Z = [[D.get((i,j),0)/the_mean for i in range(mx+1)] for j in range(my+1)]

        title = "Algo=%s, nbiter=%s, ndivs=%s" % (self.name(), n_iterations, ndivs)
        filename = 'mesure_%s_iter%s_div%s.png' % (self.name(), n_iterations, ndivs)

        self._plot_wideframe(X,Y,Z,title,filename)


    def plot_invariant_measure_inverse(self, n_iterations, ndivs, norm='sup'):
        r"""
        Return a matplotlib graph of the inverse of the invariant measure.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``ndvis`` - integer, number of divisions per dimension
        - ``norm`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points

        EXAMPLES::

            sage: algo.arrevert.plot_invariant_measure_inverse(1000000, 80) # not tested
            Creation du fichier mesure_arrevert_iter1000000_div80.png

        """
        D = self.invariant_measure_dict(n_iterations, ndivs, norm=norm)
        mx,my = map(max, zip(*D.keys()))
        len_D = float(len(D))
        the_mean = n_iterations / len_D

        E = dict((key,the_mean/value) for key,value in D.iteritems())

        X = [[i for i in range(mx+1)] for j in range(my+1)]
        Y = [[j for i in range(mx+1)] for j in range(my+1)]
        Z = [[E.get((i,j),0) for i in range(mx+1)] for j in range(my+1)]

        title = "Algo=%s, nbiter=%s, ndivs=%s" % (self.name(), n_iterations, ndivs)
        filename = 'mesure_inversed_%s_iter%s_div%s.png' % (self.name(), n_iterations, ndivs)

        self._plot_wideframe(X,Y,Z,title,filename)

    def _plot_wideframe(self, X,Y,Z, title, filename):
        r"""
        EXAMPLES::
        """
        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        ax.set_title(title)
        plt.savefig(filename)
        print "Creation du fichier %s" % filename

    def plot_natural_extension(self, n_iterations, norm_algo='1',
            norm_ext='1', axis_off=False):
        r"""
        EXAMPLES::

            sage: algo.arp.plot_natural_extension(1000)
        """
        import matplotlib
        import matplotlib.pyplot as plt
        t = self.natural_extension(n_iterations, norm_algo=norm_algo, norm_ext=norm_ext)
        domain_in, domain_out, dual_in, dual_out = t
        c = dict(zip(domain_in.keys(), ['b','r','g','c','m','y','k']))

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

        create_sub_plot(ax, domain_in,  "Algo IN",    norm=norm_algo)
        create_sub_plot(bx, dual_in,    "NatExt IN",  norm=norm_ext)
        create_sub_plot(cx, domain_out, "Algo OUT",   norm=norm_algo)
        create_sub_plot(dx, dual_out,   "NatExt OUT", norm=norm_ext)

        title = "Algo=%s, nbiter=%s" % (self.name(), n_iterations)
        fig.suptitle(title)
        filename = 'nat_ext_%s_iter%s.png' % (self.name(), n_iterations)
        plt.savefig(filename)
        #plt.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)
        print "Creation du fichier %s" % filename

    def tikz_natural_extension(self, n_iterations, norm_algo='1', norm_ext='1'):
        r"""
        EXAMPLES::

            sage: from slabbe.mcf import algo
            sage: s = algo.brun.tikz_natural_extension(1000)
            sage: view(s, tightpage=True)
        """

        t = self.natural_extension(n_iterations, norm_algo=norm_algo, norm_ext=norm_ext)
        domain_in, domain_out, dual_in, dual_out = t
        sqrt3 = 1.73205080756888
        r = 1.08
        color_dict = dict(zip(domain_in.keys(), PGF_COLORS))

        s = ''
        s += "\\begin{tikzpicture}[scale=.7]\n"
        s += ("\\begin{groupplot}[group style={group size=4 by 1},"
               "height=7cm,width=8cm,"
               "xmin=-1.1,xmax=1.1,ymin=-.6,ymax=1.20,"
               "hide axis]\n")
        for P in [domain_in, dual_in, domain_out, dual_out]:
            s += "\\nextgroupplot\n"
            for key,value in P.iteritems():
                s += ("\\draw[dashed] "
                 "(axis cs:%s, %s)" % (-r*sqrt3/2,r*-.5) + " node[left]  {$e_1$} -- "
                 "(axis cs:%s, %s)" % (r*sqrt3/2,r*-.5)  + " node[right] {$e_2$} -- "
                 "(axis cs:%s, %s)" % (0, r)             + " node[above] {$e_3$} -- cycle;\n")
                s += "\\addplot+[only marks,mark=*,mark options={color=%s}] " % color_dict[key]
                s += "coordinates {%s};\n" % '\n'.join(map(str, value))
                s += "\\addlegendentry{%s}\n " % key
        s += "\\end{groupplot}\n"
        s += "\\draw[draw=none] (group c1r1.center) -- node {$\\times$} (group c2r1.center);\n"
        s += "\\draw[draw=none] (group c2r1.center) -- node {$\\to$} (group c3r1.center);\n"
        s += "\\draw[draw=none] (group c3r1.center) -- node {$\\times$} (group c4r1.center);\n"
        s += "\\end{tikzpicture}\n"
        return LatexExpr(s)

    def sample_lyapounov_exponent(self, ntimes, n_iterations=1000):
        r"""
        Return two lists of values for theta1 and theta2

        INPUT:

        - ``ntimes`` -- integer, number of orbits
        - ``n_iterations`` -- integer, length of each orbit

        OUTPUT:

        EXAMPLES::

            sage: algo.brun.sample_lyapounov_exponent(10, 100000)

        """
        Theta1 = []
        Theta2 = []
        Uniform = []
        for _ in range(ntimes):
            l1,l2,l2surl1 = self.lyapounov_exponents(n_iterations)
            Theta1.append(l1)
            Theta2.append(l2)
            Uniform.append(1-l2surl1)
        return Theta1, Theta2, Uniform

    def table_lyapounov_exponent(self, ntimes, n_iterations=1000):
        r"""
        Return a table of values of 1 - theta2/theta1.

        INPUT:

        - ``ntimes`` -- integer, number of orbits
        - ``n_iterations`` -- integer, length of each orbit

        OUTPUT:

        liapounov exponents

        EXAMPLES::

            sage: algo.brun.table_lyapounov_exponent(10, 100000)
            n          | 10
            iterations | 100000
            min        | -0.369748622574
            mean       | -0.367620100939
            max        | -0.364250082762
            std        | 0.00155230985651

        """
        L = self.sample_lyapounov_exponent(ntimes, n_iterations)
        import numpy as np
        a = np.array(L)
        header = ['n', 'iterations', 'min','mean','max','std']
        cols = (ntimes, n_iterations, a.min(), a.mean(), a.max(), a.std()),
        t = table(columns=cols,header_column=header)
        return t

    def plot_dual_domain(self):
        r"""
        Return a plot of the dual domain
        
        EXAMPLES::

            sage: algo.armonteil.plot_dual_domain()
            Creation du fichier transposed_mat_armonteil.png
        """
        from sage.plot.point import point2d
        from sage.plot.colors import hue
        from sage.plot.graphics import Graphics
        from sage.plot.polygon import polygon
        from sage.plot.line import line
        from sage.plot.text import text
        from sage.misc.latex import latex
        from sage.matrix.constructor import matrix 
        from sage.modules.free_module_element import vector
        from sage.rings.rational_field import QQ
        D = self.dual_domain()
        len_D = len(D) * 1.
        G = Graphics()
        for i, (key, value) in enumerate(D.iteritems()):
            color = hue(i/len_D)
            value2d = [pt[:2] for pt in value]
            G += polygon(value2d, color=color, legend_label=str(key))
            G += point2d(value2d, color=color, legend_label=str(key), size=200)
            center = sum(map(vector, value2d)) / 3
            M = matrix(QQ, value).transpose()
            G += text("%s:\n%s" %(key,M), center, color='black')
        title = "Transposed matrices of algo=%s" % self.name()
        filename = 'transposed_mat_%s.png' % self.name()
        G.save(filename, title=title)
        print "Creation du fichier %s" % filename
        #return G


# Les algos
class Algo(object):
    sorted_brun = MCFAlgorithm('sorted_brun')
    sorted_arp = MCFAlgorithm('sorted_arp')
    sorted_arrevert = MCFAlgorithm('sorted_arrevert')
    sorted_selmer = MCFAlgorithm('sorted_selmer')
    sorted_meester = MCFAlgorithm('sorted_meester')
    sorted_poincare = MCFAlgorithm('sorted_poincare')
    sorted_arnouxrauzy = MCFAlgorithm('sorted_arnouxrauzy')
    brun = MCFAlgorithm('brun')
    brunmulti = MCFAlgorithm('brunmulti')
    selmer = MCFAlgorithm('selmer')
    meester = MCFAlgorithm('meester')
    poincare = MCFAlgorithm('poincare')
    arp = MCFAlgorithm('arp')
    arpmulti = MCFAlgorithm('arpmulti')
    arrevert = MCFAlgorithm('arrevert')
    arrevertmulti = MCFAlgorithm('arrevertmulti')
    armonteil = MCFAlgorithm('armonteil')
    delaunay = MCFAlgorithm('delaunay')
    jacobi = MCFAlgorithm('jacobi')
    jacobiadditif = MCFAlgorithm('jacobiadditif')
    jacobiadditifv2 = MCFAlgorithm('jacobiadditifv2')
    def __iter__(self):
        L = [   self.brun,
                self.brunmulti,     self.selmer,    self.meester,
                self.poincare,      self.arp,       self.arpmulti, self.arrevert,
                self.arrevertmulti, self.armonteil, self.jacobi,
                self.jacobiadditif, self.jacobiadditifv2]
        return iter(L)
algo = Algo()
