# -*- coding: utf-8 -*-
r"""
Multidimensional Continued Fraction Algorithms

EXAMPLES::

    sage: from slabbe.mcf import algo
    sage: algo.brun.plot_natural_extension(3000, norm_left='1', axis_off=True)
    Creation du fichier nat_ext_brun_iter3000.png

::

    sage: from slabbe.mcf import algo
    sage: from itertools import repeat
    sage: D = algo.arp.substitutions()
    sage: it = algo.arp.coding_iterator((1,e,pi))
    sage: words.s_adic(it, repeat(1), D)
    word: 1232323123233231232332312323123232312323...

TODO:

 - Move some code to the cython part of this module (grep TODO)

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

            sage: from slabbe.mcf import algo
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

            sage: from slabbe.mcf import algo
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

    def plot_natural_extension(self, n_iterations, norm_left='1',
            norm_right='1', axis_off=False):
        r"""
        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``norm_left`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``norm_right`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``axis_off`` - boolean, 

        EXAMPLES::

            sage: from slabbe.mcf import algo
            sage: algo.sorted_arp.plot_natural_extension(1000)
            Creation du fichier nat_ext_sorted_arp_iter1000.png
        """
        import matplotlib
        import matplotlib.pyplot as plt
        t = self.natural_extension(n_iterations, norm_left=norm_left,
                norm_right=norm_right)
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

        create_sub_plot(ax, domain_left,  "Algo IN",    norm=norm_left)
        create_sub_plot(bx, domain_right,    "NatExt IN",  norm=norm_right)
        create_sub_plot(cx, image_left, "Algo OUT",   norm=norm_left)
        create_sub_plot(dx, image_right,   "NatExt OUT", norm=norm_right)

        title = "Algo=%s, nbiter=%s" % (self.name(), n_iterations)
        fig.suptitle(title)
        filename = 'nat_ext_%s_iter%s.png' % (self.name(), n_iterations)
        plt.savefig(filename)
        #plt.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)
        print "Creation du fichier %s" % filename

    def tikz_natural_extension(self, n_iterations, norm_left='1',
            norm_right='1', marksize=0.2, legend_marksize=2):
        r"""

        INPUT:

        - ``n_iterations`` -- integer, number of iterations
        - ``norm_left`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``norm_right`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``marksize`` -- tikz marksize (default:``0.2``)
        - ``legend_marksize`` -- tikz legend marksize (default:``2``)

        EXAMPLES::

            sage: from slabbe.mcf import algo
            sage: s = algo.brun.tikz_natural_extension(1000)
            sage: view(s, tightpage=True)   # not tested
        """

        t = self.natural_extension(n_iterations, norm_left=norm_left,
                norm_right=norm_right)
        domain_left, image_left, domain_right, image_right = t
        sqrt3 = 1.73205080756888
        r = 1.08
        color_dict = dict(zip(domain_left.keys(), PGF_COLORS))

        s = ''
        s += "\\begin{tikzpicture}[scale=.7]\n"
        s += ("\\begin{groupplot}[group style={group size=4 by 1},"
               "height=7cm,width=8cm,"
               "xmin=-1.1,xmax=1.1,ymin=-.6,ymax=1.20,"
               "hide axis]\n")
        for data in [domain_left, domain_right, image_left, image_right]:
            s += "\\nextgroupplot\n"
            s += ("\\draw[dashed] "
                  "(axis cs:%s, %s)" % (-r*sqrt3/2,r*-.5) +
                  " node[left] {$\\mathbf{e}_1$} -- \n"
                  "(axis cs:%s, %s)" % (r*sqrt3/2,r*-.5)  +
                  " node[right] {$\\mathbf{e}_2$} -- \n"
                  "(axis cs:%s, %s)" % (0, r)             +
                  " node[above] {$\\mathbf{e}_3$} -- cycle;\n")
            for key,value in data.iteritems():
                s += "\\addplot+["
                s += "legend image post style={mark size=%s}," % legend_marksize
                s += "only marks,mark=*,"
                s += "mark size=%s," % marksize
                s += "mark options={color=%s}] " % color_dict[key]
                s += "coordinates {%s};\n" % '\n'.join(map(str, value))
                s += "\\addlegendentry{%s}\n " % key
        s += "\\end{groupplot}\n"
        s += "\\draw[draw=none] (group c1r1.center) -- node {$\\times$} (group c2r1.center);\n"
        s += "\\draw[draw=none] (group c2r1.center) -- node {$\\to$} (group c3r1.center);\n"
        s += "\\draw[draw=none] (group c3r1.center) -- node {$\\times$} (group c4r1.center);\n"
        s += "\\end{tikzpicture}\n"
        return LatexExpr(s)

    def tikz_natural_extension_part(self, n_iterations, part=3, 
                                    norm_left='1', norm_right='1',
                                    marksize='1pt', limit_nb_points=None,
                                    verbose=False):
        r"""
        Return a pgfplots or some part of an orbit in the natural
        extension.

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``part`` - integer, taking value 0, 1, 2 or 3
        - ``norm_left`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``norm_right`` -- string (default: ``'sup'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``marksize`` -- string (default: ``'1pt'``), pgfplots mark size value
        - ``limit_nb_points`` -- None or integer (default: ``None``), limit
          number of points per patch
        - ``verbose`` -- string (default: ``False``)

        EXAMPLES::

            sage: from slabbe.mcf import algo
            sage: s = algo.arp.tikz_natural_extension_part(1000, part=3)
            sage: view(s, tightpage=True)    # not tested

        also ::

            sage: s = algo.arp.tikz_natural_extension_part(1000, part=3, marksize='.2pt', limit_nb_points=1200, verbose=True)
            Taking ... points for key 32
            Taking ... points for key 1
            Taking ... points for key 2
            Taking ... points for key 3
            Taking ... points for key 12
            Taking ... points for key 13
            Taking ... points for key 21
            Taking ... points for key 23
            Taking ... points for key 31
        """
        t = self.natural_extension(n_iterations, norm_left=norm_left,
                norm_right=norm_right)
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
                    print "Taking only {} points instead of {} for key {}".format(
                            limit_nb_points, len(value), key)
                value = value[:limit_nb_points]
            elif verbose:
                print "Taking {} points for key {}".format(len(value), key)

            s += "\\addplot+[only marks,mark=*,mark options={color=%s}," % color_dict[key]
            s += "mark size=%s]\n" % marksize
            s += "coordinates {\n"
            s += '\n'.join(map(str, value))
            s += "};\n" 
            s += "\\addlegendentry{%s}\n " % key
        s += "\\end{axis}\n"
        s += "\\end{tikzpicture}\n"
        return LatexExpr(s)

    def png_natural_extension_part(self, n_iterations, draw,
                                    norm_left='1', norm_right='1',
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

        - ``norm_left`` -- string (default: ``'1'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``norm_right`` -- string (default: ``'1'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
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

            sage: from slabbe.mcf import algo
            sage: c = {}
            sage: c[1] = c[2] = c[3] = [0,0,0]
            sage: c[12] = c[13] = c[23] = c[21] = c[31] = c[32] = [255,0,0]
            sage: b = [1,2,3,12,13,21,23,31,32]
            sage: P = algo.arp.png_natural_extension_part(10^7, draw='image_right',  # not tested
            ....:   branch_order=b, color_dict=c, urange=(-.6,.6), vrange=(-.6,.6))  # not tested

        Half a minute (27s) for a picture zoomed in the orbit of 10^8
        points::

            sage: P = algo.arp.png_natural_extension_part(10^8, draw='image_right', # not tested
            ....:   branch_order=b, color_dict=c, urange=(.2,.3), vrange=(.2,.3))   # not tested

        EXAMPLES::

            sage: from slabbe.mcf import algo
            sage: c = {}
            sage: c[1] = c[2] = c[3] = [0,0,0]
            sage: c[12] = c[13] = c[23] = c[21] = c[31] = c[32] = [255,0,0]
            sage: b = [1,2,3,12,13,21,23,31,32]
            sage: opt = dict(urange=(-.6,.6), vrange=(-.6,.6), color_dict=c, branch_order=b)
            sage: P = algo.arp.png_natural_extension_part(10^5, draw='domain_left', **opt)
            sage: P = algo.arp.png_natural_extension_part(10^5, draw='domain_right', **opt)
            sage: P = algo.arp.png_natural_extension_part(10^5, draw='image_left', **opt)
            sage: P = algo.arp.png_natural_extension_part(10^5, draw='image_right', **opt)
            sage: P.show() # not tested

        """
        L = self.orbit_filtered_list(n_iterations,
            norm_left=norm_left, norm_right=norm_right,
            xmin=xrange[0], xmax=xrange[1],
            ymin=yrange[0], ymax=yrange[1],
            umin=urange[0], umax=urange[1],
            vmin=vrange[0], vmax=vrange[1],
            ndivs=ndivs)
        if branch_order is None:
            raise NotImplementedError
            branch_order = []
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

    def measure_evalution(self, n_iterations, draw,
                                norm_left='1', norm_right='1',
                                xrange=(-.866, .866),
                                yrange=(-.5, 1.),
                                urange=(-.866, .866),
                                vrange=(-.5, 1.),
                                ndivs=1024,
                                verbose=False):
        r"""

        INPUT:

        - ``n_iterations`` - integer, number of iterations
        - ``draw`` -- string (default: ``'image_right'``), possible values
          are:

          - ``'domain_left'`` - use x and y ranges
          - ``'domain_right'`` - use u and v ranges
          - ``'image_left'`` - use x and y ranges
          - ``'image_right'`` - use u and v ranges

        - ``norm_left`` -- string (default: ``'1'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
        - ``norm_right`` -- string (default: ``'1'``), either ``'sup'`` or
          ``'1'``, the norm used for the orbit points
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

            sage: from slabbe.mcf import algo
            sage: opt = dict(urange=(-.15,.25), vrange=(-.05,.05))
            sage: algo.arp.measure_evalution(10^8, draw='right', ndivs=100, **opt) # optional long
            0.435...
            sage: algo.arp.measure_evalution(10^8, draw='right', ndivs=1000, **opt) # optional long
            0.357...
            sage: algo.arp.measure_evalution(10^8, draw='right', ndivs=2000, **opt) # optional long
            0.293...
            sage: algo.arp.measure_evalution(10^8, draw='right', ndivs=4000, **opt) # optional long
            0.177...

        """
        L = self.orbit_filtered_list(n_iterations,
            norm_left=norm_left, norm_right=norm_right,
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
            print "nombre diterations dans la fenetre : ", len(L)
            print "{} pixels touchés parmi limage {}^2 ".format(len(S), ndivs)

        print float(len(S) / ndivs**2)

    def sample_lyapounov_exponent(self, ntimes, n_iterations=1000):
        r"""
        Return two lists of values for theta1 and theta2

        INPUT:

        - ``ntimes`` -- integer, number of orbits
        - ``n_iterations`` -- integer, length of each orbit

        OUTPUT:

        EXAMPLES::

            sage: from slabbe.mcf import algo
            sage: T1, T2, U = algo.brun.sample_lyapounov_exponent(10, 100000)
            sage: T1[0]
            0.303680940345907
            sage: T2[0]
            -0.1119022164698561 
            sage: U[0]
            1.3684861366090153

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

            sage: from slabbe.mcf import algo
            sage: algo.sorted_brun.table_lyapounov_exponent(10, 100000)
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

    def plot_dual_domain(self, savefile=False):
        r"""
        Return a plot of the dual domain
        
        EXAMPLES::

            sage: from slabbe.mcf import algo
            sage: algo.sorted_armonteil.plot_dual_domain()
            Graphics object consisting of 18 graphics primitives
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
        if savefile:
            title = "Transposed matrices of algo=%s" % self.name()
            filename = 'transposed_mat_%s.png' % self.name()
            G.save(filename, title=title)
            print "Creation of the file %s" % filename
        return G

# Les algos
class Algo(object):
    sorted_brun = MCFAlgorithm('sorted_brun')
    sorted_arp = MCFAlgorithm('sorted_arp')
    sorted_arrevert = MCFAlgorithm('sorted_arrevert')
    sorted_selmer = MCFAlgorithm('sorted_selmer')
    sorted_meester = MCFAlgorithm('sorted_meester')
    sorted_poincare = MCFAlgorithm('sorted_poincare')
    sorted_arnouxrauzy = MCFAlgorithm('sorted_arnouxrauzy')
    sorted_armonteil = MCFAlgorithm('sorted_armonteil')
    brun = MCFAlgorithm('brun')
    brunmulti = MCFAlgorithm('brunmulti')
    selmer = MCFAlgorithm('selmer')
    meester = MCFAlgorithm('meester')
    poincare = MCFAlgorithm('poincare')
    arp = MCFAlgorithm('arp')
    arpmulti = MCFAlgorithm('arpmulti')
    arrevert = MCFAlgorithm('arrevert')
    arrevertmulti = MCFAlgorithm('arrevertmulti')
    delaunay = MCFAlgorithm('delaunay')
    jacobi = MCFAlgorithm('jacobi')
    jacobiadditif = MCFAlgorithm('jacobiadditif')
    jacobiadditifv2 = MCFAlgorithm('jacobiadditifv2')
    def __iter__(self):
        L = [   self.brun,
                self.brunmulti,     self.selmer,    self.meester,
                self.poincare,      self.arp,       self.arpmulti, self.arrevert,
                self.arrevertmulti, self.sorted_armonteil, self.jacobi,
                self.jacobiadditif, self.jacobiadditifv2]
        return iter(L)
algo = Algo()
