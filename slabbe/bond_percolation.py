# coding=utf-8
r"""
Bond Percolation

This is an implementation of bond percolation. See Chapter 3 of [POG]_. See
also my blog post `Percolation and self-avoiding walks`__ related to this
code.

__ http://www.liafa.univ-paris-diderot.fr/~labbe/blogue/2012/12/percolation-and-self-avoiding-walks/

AUTHORS:

- Sebastien Labbe (2012-12-17): initial version, for `pog`__ lecture group

__ http://www.liafa.univ-paris-diderot.fr/~jolivet/pog/

REFERENCES:

.. [POG] Geoffrey Grimmett, Probability on Graphs,
   http://www.statslab.cam.ac.uk/~grg/books/pgs.html

EXAMPLES:

Bond percolation sample
-----------------------

We construct a bond percolation sample in dimension d=2 with probability of
open edges p=0.3::

    sage: from slabbe import BondPercolationSample
    sage: S = BondPercolationSample(p=0.3, d=2)
    sage: S
    Bond percolation sample d=2 p=0.300

An edge is defined uniquely as a starting point in Z^d and an axis
direction given by an integer i such that 1 <= i <= d. One may ask if a
given edge is in the sample S::

    sage: ((34,56), 2) in S     # random
    False

The result is cached so the same answer is returned again::

    sage: ((34,56), 2) in S     # random
    False
    sage: ((34,56), 2) in S     # random
    False

The cluster containing the point zero is returned as an iterator::

    sage: S.cluster()
    <generator object at ...>

It may be finite of infinite. If you believe it is finite, you may compute
its cardinality. If the cluster is infinite, it will not halt::

    sage: S.cluster_cardinality() # not tested, might not halt
    13

For larger values of p, the cluster might be larger if not infinite. In
this case you may want to stop the computation at a certain point
determined in advance. The following method does this. And it returns the
cardinality if it is smaller than the stop value::

    sage: S = BondPercolationSample(p=0.45, d=2)
    sage: S.cluster_cardinality_stop_at(stop=10)          # random
    '>=10'
    sage: S.cluster_cardinality_stop_at(stop=100)          # random
    '>=100'
    sage: S.cluster_cardinality_stop_at(stop=1000)          # random
    625

Bond percolation samples
------------------------

Construction of 20 bond percolation samples. For each of them, compute the
cardinality of the open cluster containing zero::

    sage: from slabbe import BondPercolationSamples
    sage: S20 = BondPercolationSamples(p=0.4, d=2, n=20)
    sage: S20.cluster_cardinality(stop=100)               # random
    [4, 2, 1, 4, 1, 10, 62, 71, 1, 25, 19, 2, 2, 42, '>=100', 1, 18, 2, '>=100', 20]

By considering "larger than 100" to be an infinite cluster, this gives a
value of 2/20 = 0.10 for the percolation probability::

    sage: S20.percolation_probability(stop=100)        # random
    0.100

By increasing the stop value, the computations can be redone again *on the
same samples*. In this case, by using a stop value of 1000, we get a value
of 0 for the percolation probability of p=0.4::

    sage: S20.cluster_cardinality(stop=1000)             # random
    [4, 2, 1, 4, 1, 10, 62, 71, 1, 25, 19, 2, 2, 42, 176, 1, 18, 2, 186, 20]
    sage: S20.percolation_probability(stop=1000)         # random
    0.000

Percolation probability
-----------------------

One can define the percolation probability function for a given dimension
d. It will generate n samples and consider the cluster to be infinite if
its cardinality is larger than the given stop value::

    sage: from slabbe import PercolationProbability
    sage: T = PercolationProbability(d=2, n=10, stop=100)
    sage: T
    Percolation Probability $\theta(p)$
    d = dimension = 2
    n = # samples = 10
    stop counting at = 100

Compute the value for a certain probability p::

    sage: T(0.4534)             # random
    0.300

Of course, this value will change for another equal percolation
probability since it depends on the samples::

    sage: T = PercolationProbability(d=2, n=10, stop=100)
    sage: T(0.4534)              # random
    0.600

Anyway, it is usefull to draw the plot of the percolation probability::

    sage: T = PercolationProbability(d=2, n=10, stop=100)
    sage: T.return_plot((0,1), adaptive_recursion=4, plot_points=4)     # optional long
    Graphics object consisting of 2 graphics primitives

Here we use Sage adaptative recursion algorithm for drawing plots which
finds the particular important intervals to ask for more values of the
function. See help section of plot function for details. Because T might be
long to compute we start with only 4 points

TODO
----

 - Make it 100% doctested (presently 21/24 = 87%)
 - Base it on DiscreteSubset code
 - Fix tikz2pdf use

Do we want to use?::

    sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
    sage: from itertools import count
    sage: S = EnumeratedSetFromIterator(count)
    sage: S
    {0, 1, 2, 3, 4, ...}

and ?::

    sage: M = FiniteSetMaps(["a", "b"], [3, 4, 5]); M
    Maps from {'a', 'b'} to {3, 4, 5}

Methods and classes
-------------------

"""
#*****************************************************************************
#       Copyright (C) 2012 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
import os
import itertools
from sage.structure.sage_object import SageObject
from sage.combinat.backtrack import TransitiveIdeal
from sage.rings.infinity import Infinity
from sage.graphs.graph import Graph
from sage.functions.generalized import sgn
from sage.functions.other import abs
from sage.misc.functional import numerical_approx
from sage.misc.cachefunc import cached_method
from sage.misc.latex import LatexExpr
from sage.plot.graphics import Graphics
from sage.plot.point import point
from sage.plot.circle import circle
from sage.plot.line import line
from sage.plot.text import text
from sage.plot.plot import graphics_array, plot


class BondPercolationSample(SageObject):
    r"""
    Let $L^d = (Z^d,E^d)$ be the hypercubic lattice.

    A sample contained in the set {0,1}^{E^d}.

    An edge e in E is open (=1) in the sample with probability p.

    Cached __contains__ method below does the job of memory.
    """
    def __init__(self, p, d=2):
        r"""
        INPUT:

        - ``p`` - real number in [0,1]
        """
        self._p = p
        self._dimension = d

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import BondPercolationSample
            sage: S = BondPercolationSample(0.4,2)
            sage: S
            Bond percolation sample d=2 p=0.400
        """
        s = "Bond percolation sample d=%s p=%.3f" % (self._dimension, self._p)
        return s

    @cached_method
    def __contains__(self, arg):
        r"""
        Return True with probability p.

        INPUT:

        - ``arg`` - a 2-tuple containing:

          - ``pt`` - point in Z^d
          - ``direction`` - integer, possible values are 1, 2, ..., d and -1,
            -2, ..., -d.

        EXAMPLES::

            sage: from slabbe import BondPercolationSample
            sage: S = BondPercolationSample(0.5)
            sage: (S.zero(), 1) in S        # random
            True
            sage: (S.zero(), -1) in S        # random
            True
        """
        pt, direction = arg
        # here because creates docbuild error when the import is global
        from random import random 
        return random() < self._p

    def neighbor(self, pt, d):
        r"""
        Return the neighbors of the point pt in direction d.

        INPUT:

        - ``pt`` - tuple, point in Z^d
        - ``direction`` - integer, possible values are 1, 2, ..., d and -1,
          -2, ..., -d.

        EXAMPLES:

            sage: from slabbe import BondPercolationSample
            sage: S = BondPercolationSample(0.5,2)
            sage: S.neighbor((2,3),1)
            (3, 3)
            sage: S.neighbor((2,3),2)
            (2, 4)
            sage: S.neighbor((2,3),-1)
            (1, 3)
            sage: S.neighbor((2,3),-2)
            (2, 2)

        """
        R = xrange(self._dimension)
        a = sgn(d)
        d = abs(d)
        return tuple(pt[k]+a if k==d-1 else pt[k] for k in R)

    def children(self, pt):
        r"""
        Return an iterator over open neighbors of the point pt.

        INPUT:

        - ``pt`` - tuple, point in Z^d
        - ``m`` - integer, limit

        EXAMPLES:

        The result is consistent::

            sage: from slabbe import BondPercolationSample
            sage: S = BondPercolationSample(0.5)
            sage: list(S.children((0,0)))        # random
            [(1, 0), (-1, 0), (0, -1)]

        Might be different for another sample::

            sage: S = BondPercolationSample(0.5)
            sage: list(S.children((0,0)))        # random
            [(-1, 0)]

        In dimension 3::

            sage: S = BondPercolationSample(0.5,3)
            sage: list(S.children((0,0,0)))        # random
            [(1, 0, 0), (0, -1, 0), (0, 0, -1)]
            sage: S = BondPercolationSample(0.5,3)
            sage: list(S.children((0,0,0)))        # random
            [(1, 0, 0), (-1, 0, 0)]
            sage: list(S.children((0,0,0)))        # random
            [(1, 0, 0), (-1, 0, 0)]

        ::

            sage: S = BondPercolationSample(1,2)
            sage: list(S.children((0,0)))        # random
            [(1, 0), (-1, 0), (0, 1), (0, -1)]
            sage: S = BondPercolationSample(1,3)
            sage: list(S.children((0,0,0)))        # random
            [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
        """
        for d in xrange(1, self._dimension+1):
            if (pt, d) in self:
                yield self.neighbor(pt, d)
            opposite_pt = self.neighbor(pt, -d)
            if (opposite_pt, d) in self:
                yield opposite_pt

    def zero(self):
        r"""
        EXAMPLES::

            sage: from slabbe import BondPercolationSample
            sage: S = BondPercolationSample(0.5,3)
            sage: S.zero()
            (0, 0, 0)

        ::

            sage: S = BondPercolationSample(0.5,5)
            sage: S.zero()
            (0, 0, 0, 0, 0)
        """
        return (0,) * self._dimension

    def cluster(self, pt=None):
        r"""
        Return an iterator over the open cluster containing the point pt.

        INPUT:

        - ``pt`` - tuple, point in Z^d. If None, pt=zero is considered.

        EXAMPLES:

            sage: from slabbe import BondPercolationSample
            sage: S = BondPercolationSample(0.5)
            sage: it = S.cluster()
            sage: next(it)
            (0, 0)

        """
        if pt is None:
            pt = self.zero()
        generators = [pt]
        return iter(TransitiveIdeal(self.children, generators))

    def cluster_cardinality(self, pt=None):
        r"""
        INPUT:

        - ``pt`` - tuple, point in Z^d. If None, pt=zero is considered.

        EXAMPLES::

            sage: from slabbe import BondPercolationSample
            sage: BondPercolationSample(0.01).cluster_cardinality() # random
            1
            sage: BondPercolationSample(0.4).cluster_cardinality()  # random
            28
        """
        return sum(1 for _ in self.cluster(pt=pt))

    def cluster_cardinality_stop_at(self, stop, pt=None):
        r"""
        Return the cardinality of the cluster or the strin ">=STOP" if the size is
        larger than stop value.

        INPUT:

        - ``stop`` - integer
        - ``pt`` - tuple, point in Z^d. If None, pt=zero is considered.

        EXAMPLES::

            sage: from slabbe import BondPercolationSample
            sage: S = BondPercolationSample(0.3,2)
            sage: S.cluster_cardinality()                # random
            13
            sage: S.cluster_cardinality_stop_at(1000)    # random
            13
            sage: S.cluster_cardinality_stop_at(100)     # random
            13
            sage: S.cluster_cardinality_stop_at(10)      # random
            '>=10'
        """
        it = itertools.islice(self.cluster(pt=pt), stop)
        size = sum(1 for _ in it)
        if size == stop:
            return ">=%s" % stop
        else:
            return size

    def edges_in_box(self, m):
        r"""
        Return an iterator over all edges in the primal box [-m,m]^d.

        INPUT:

        - ``m`` - integer

        EXAMPLES::

            sage: from slabbe import BondPercolationSample
            sage: S = BondPercolationSample(1,2)
            sage: for a in S.edges_in_box(1): print(a)
            ((-1, -1), (0, -1))
            ((-1, -1), (-1, 0))
            ((-1, 0), (0, 0))
            ((-1, 0), (-1, 1))
            ((0, -1), (1, -1))
            ((0, -1), (0, 0))
            ((0, 0), (1, 0))
            ((0, 0), (0, 1))
        """
        R = xrange(-m, m)
        L = [R]*self._dimension
        positive_directions = xrange(1, self._dimension+1)
        for pt in itertools.product(*L):
            for d in positive_directions:
                if (pt,d) in self:
                    neighbor = self.neighbor(pt,d)
                    yield (pt, neighbor)

    def cluster_in_box(self, m, pt=None):
        r"""
        Return the cluster (as a list) in the primal box [-m,m]^d
        containing the point pt.

        INPUT:

        - ``m`` - integer
        - ``pt`` - tuple, point in Z^d. If None, pt=zero is considered.

        EXAMPLES::

            sage: from slabbe import BondPercolationSample
            sage: S = BondPercolationSample(0.3,2)
            sage: S.cluster_in_box(2)         # random
            [(-2, -2), (-2, -1), (-1, -2), (-1, -1), (-1, 0), (0, 0)]
        """
        G = Graph()
        G.add_edges(self.edges_in_box(m))
        if pt is None:
            pt = self.zero()
        if pt in G:
            return G.connected_component_containing_vertex(pt)
        else:
            return []

    def plot(self, m, pointsize=100, thickness=3, axes=False):
        r"""
        Return 2d graphics object contained in the primal box [-m,m]^d.

        INPUT:

        - ``pointsize``, integer (default:``100``),
        - ``thickness``, integer (default:``3``),
        - ``axes``, bool (default:``False``),

        EXAMPLES::

            sage: from slabbe import BondPercolationSample
            sage: S = BondPercolationSample(0.5,2)
            sage: S.plot(2)           # optional long

        It works in 3d!!::

            sage: S = BondPercolationSample(0.5,3)
            sage: S.plot(3, pointsize=10, thickness=1)     # optional long
            Graphics3d Object

        """
        s = ""
        s += "\\begin{tikzpicture}\n"
        s += "[inner sep=0pt,thick,\n"
        s += "reddot/.style={fill=red,draw=red,circle,minimum size=5pt}]\n"
        s += "\\clip %s rectangle %s;\n" % ((-m-.4,-m-.4), (m+.4,m+.4))
        G = Graphics()
        for u in self.cluster_in_box(m+1):
            G += point(u, color='blue', size=pointsize)
        for (u,v) in self.edges_in_box(m+1):
            G += line((u,v), thickness=thickness, alpha=0.8)
        G += text("p=%.3f" % self._p, (0.5,1.03), axis_coords=True, color='black')
        G += circle((0,0), 0.5, color='red', thickness=thickness)
        if self._dimension == 2:
            G.axes(axes)
        return G

    def tikz(self, m):
        r"""
        Return tikz code.

        EXAMPLES::

            sage: from slabbe import BondPercolationSample
            sage: S = BondPercolationSample(0.5,2)
            sage: S.tikz(2)
            \begin{tikzpicture}
            [inner sep=0pt,thick, 
             reddot/.style={fill=red,draw=red,circle,minimum size=5pt}]
            \clip (-2.4, -2.4) rectangle (2.4, 2.4);
            \draw (..., ...) -- (..., ...);
            ...
            \node[circle,fill=none,draw=red,minimum size=0.8cm,ultra thick,inner sep=0pt] at (0,0) {};
            \node[above right] at (0,0) {$(0, 0)$};
            \end{tikzpicture}

        """
        s = ""
        s += "\\begin{tikzpicture}\n"
        s += "[inner sep=0pt,thick,\n"
        s += "reddot/.style={fill=red,draw=red,circle,minimum size=5pt}]\n"
        s += "\\clip %s rectangle %s;\n" % ((-m-.4,-m-.4), (m+.4,m+.4))
        for (u,v) in self.edges_in_box(m+1):
            s += "\\draw %s -- %s;\n" % (u,v)
        for u in self.cluster_in_box(m+1):
            s += "\\node[reddot] at %s {};\n" % (u,)
        s += "\\node[circle,fill=none,draw=red,minimum size=0.8cm,ultra thick,inner sep=0pt] at (0,0) {};\n"
        s += "\\node[above right] at (0,0) {$%s$};\n" % (self.zero(),)
        s += "\\end{tikzpicture}\n"
        return LatexExpr(s)

    def save_pdf(self, m):
        r"""
        EXAMPLES::

            sage: from slabbe import BondPercolationSample
            sage: BondPercolationSample(0.3, d=2).save_pdf(20)    # optional long
            Creation du fichier tikz_sample_d2_p300_m20.tikz
            Using template '/Users/slabbe/.tikz2pdf.tex'.
            tikz2pdf: calling pdflatex...
            tikz2pdf: Output written to 'tikz_sample_d2_p300_m20.pdf'.
        """
        tikz = self.tikz(m)
        prefix = "tikz_sample_d%s_p%3d_m%s" % (self._dimension,1000*self._p, m)
        prefix = prefix.replace(' ','0')
        filename = prefix + '.tikz'
        with open(filename, 'w') as f:
            f.write(tikz)
            print("Creation du fichier %s" % filename)
        os.system("tikz2pdf %s" % filename)
        #os.system("convert %s.pdf %s.png" % (prefix, prefix))

class BondPercolationSamples(SageObject):
    r"""
    Return a list of n BondPercolationSample of given parameter p and
    dimension d.

    EXAMPLES::

        sage: from slabbe import BondPercolationSamples
        sage: BondPercolationSamples(0.2,2,3)
        <slabbe.bond_percolation.BondPercolationSamples object at ...>
    """
    def __init__(self, p, d, n):
        r"""
        """
        self._p = p
        self._dimension = d
        self._n = n
        self._list = [BondPercolationSample(p,d) for _ in range(n)]

    @cached_method
    def cluster_cardinality(self, stop):
        r"""
        Return the list of cardinality of the cluster for each sample
        or +Infinity if the size is larger than stop value.

        INPUT:

        - ``stop`` - integer

        EXAMPLES::

            sage: from slabbe import BondPercolationSamples
            sage: S20 = BondPercolationSamples(p=0.2, d=2, n=20)
            sage: S20.cluster_cardinality(100)       # random
            [1, 4, 1, 2, 2, 6, 5, 1, 5, 9, 2, 2, 2, 1, 2, 4, 4, 3, 2, 1]

        ::

            sage: d = 2
            sage: n = 5
            sage: for p in srange(0,1,0.1): print(p,BondPercolationSamples(p,d,n).cluster_cardinality(100)) # optional long
            0.000000000000000 [1, 1, 1, 1, 1]
            0.100000000000000 [5, 1, 2, 1, 2]
            0.200000000000000 [3, 1, 4, 1, 1]
            0.300000000000000 [11, 7, 2, 1, 3]
            0.400000000000000 [3, 3, 35, 10, '>=100']
            0.500000000000000 ['>=100', '>=100', '>=100', '>=100', 26]
            0.600000000000000 ['>=100', '>=100', '>=100', '>=100', '>=100']
            0.700000000000000 ['>=100', '>=100', '>=100', '>=100', '>=100']
            0.800000000000000 ['>=100', '>=100', '>=100', '>=100', '>=100']
            0.900000000000000 ['>=100', '>=100', '>=100', '>=100', '>=100']

        """
        return [S.cluster_cardinality_stop_at(stop) for S in self._list]

    def ntimes_over_size(self, stop):
        r"""
        EXAMPLES::

            sage: from slabbe import BondPercolationSamples
            sage: S = BondPercolationSamples(0.2,2,20)
            sage: S.ntimes_over_size(100)    # random
            0
            sage: S = BondPercolationSamples(0.4,2,20)
            sage: S.ntimes_over_size(100)    # random
            1
            sage: S = BondPercolationSamples(0.5,2,20)
            sage: S.ntimes_over_size(100)    # random
            17
        """
        L = self.cluster_cardinality(stop)
        return sum(1 for i in L if isinstance(i, str))

    def percolation_probability(self, stop):
        return numerical_approx(self.ntimes_over_size(stop) / self._n, digits=3)

def percolation_graphics_array(range_p, d, m, ncols=3):
    r"""
    EXAMPLES::

        sage: from slabbe.bond_percolation import percolation_graphics_array
        sage: percolation_graphics_array(srange(0.1,1,0.1), d=2, m=5)    # optional long
        sage: P = percolation_graphics_array(srange(0.45,0.55,0.01), d=2, m=5)  # optional long
        sage: P.save('array_p45_p55_m5.png')     # not tested
        sage: P = percolation_graphics_array(srange(0.45,0.55,0.01), d=2, m=10) # optional long
        sage: P.save('array_p45_p55_m10.png')    # not tested
    """
    pointsize=20
    thickness=1
    L = [BondPercolationSample(p,d).plot(m,pointsize=pointsize,thickness=thickness) for p in range_p]
    nrows = (len(range_p)-1) // ncols + 1
    return graphics_array(L, n=nrows, m=ncols)

def compute_percolation_probability(range_p, d, n, stop):
    r"""
    EXAMPLES::

        sage: from slabbe.bond_percolation import compute_percolation_probability
        sage: compute_percolation_probability(srange(0,0.8,0.1), d=2, n=5, stop=100) # random
        d = 2, n = number of samples = 5
        stop counting at = 100
        p=0.0000, Theta=0.000, if |C|< 100 then max|C|=1
        p=0.1000, Theta=0.000, if |C|< 100 then max|C|=1
        p=0.2000, Theta=0.000, if |C|< 100 then max|C|=5
        p=0.3000, Theta=0.000, if |C|< 100 then max|C|=6
        p=0.4000, Theta=0.000, if |C|< 100 then max|C|=31
        p=0.5000, Theta=1.00, if |C|< 100 then max|C|=-Infinity
        p=0.6000, Theta=1.00, if |C|< 100 then max|C|=-Infinity
        p=0.7000, Theta=1.00, if |C|< 100 then max|C|=-Infinity

    ::

        sage: range_p = srange(0,0.8,0.1)
        sage: compute_percolation_probability(range_p, d=2, n=5, stop=100) # not tested
        d = 2, n = number of samples = 5
        stop counting at = 100
        p=0.0000, Theta=0.000, if |C|< 100 then max|C|=1
        p=0.1000, Theta=0.000, if |C|< 100 then max|C|=1
        p=0.2000, Theta=0.000, if |C|< 100 then max|C|=5
        p=0.3000, Theta=0.000, if |C|< 100 then max|C|=6
        p=0.4000, Theta=0.000, if |C|< 100 then max|C|=31
        p=0.5000, Theta=1.00, if |C|< 100 then max|C|=-Infinity
        p=0.6000, Theta=1.00, if |C|< 100 then max|C|=-Infinity
        p=0.7000, Theta=1.00, if |C|< 100 then max|C|=-Infinity

    ::

        sage: range_p = srange(0.45,0.55,0.01)
        sage: compute_percolation_probability(range_p, d=2, n=10, stop=1000) # not tested
        d = 2, n = number of samples = 10
        stop counting at = 1000
        p=0.4500, Theta=0.000, if |C|< 1000 then max|C|=378
        p=0.4600, Theta=0.000, if |C|< 1000 then max|C|=475
        p=0.4700, Theta=0.000, if |C|< 1000 then max|C|=514
        p=0.4800, Theta=0.100, if |C|< 1000 then max|C|=655
        p=0.4900, Theta=0.700, if |C|< 1000 then max|C|=274
        p=0.5000, Theta=0.700, if |C|< 1000 then max|C|=975
        p=0.5100, Theta=0.700, if |C|< 1000 then max|C|=16
        p=0.5200, Theta=0.700, if |C|< 1000 then max|C|=125
        p=0.5300, Theta=0.900, if |C|< 1000 then max|C|=4
        p=0.5400, Theta=0.700, if |C|< 1000 then max|C|=6

    ::

        sage: range_p = srange(0.475,0.485,0.001)
        sage: compute_percolation_probability(range_p, d=2, n=10, stop=1000) # not tested
        d = 2, n = number of samples = 10
        stop counting at = 1000
        p=0.4750, Theta=0.200, if |C|< 1000 then max|C|=718
        p=0.4760, Theta=0.200, if |C|< 1000 then max|C|=844
        p=0.4770, Theta=0.200, if |C|< 1000 then max|C|=566
        p=0.4780, Theta=0.500, if |C|< 1000 then max|C|=257
        p=0.4790, Theta=0.200, if |C|< 1000 then max|C|=566
        p=0.4800, Theta=0.300, if |C|< 1000 then max|C|=544
        p=0.4810, Theta=0.300, if |C|< 1000 then max|C|=778
        p=0.4820, Theta=0.500, if |C|< 1000 then max|C|=983
        p=0.4830, Theta=0.300, if |C|< 1000 then max|C|=473
        p=0.4840, Theta=0.500, if |C|< 1000 then max|C|=411

    ::

        sage: range_p = srange(0.47,0.48,0.001)
        sage: compute_percolation_probability(range_p, d=2, n=20, stop=2000)  # not tested
        d = 2, n = number of samples = 20
        stop counting at = 2000
        p=0.4700, Theta=0.0500, if |C|< 2000 then max|C|=1666
        p=0.4710, Theta=0.100, if |C|< 2000 then max|C|=1665
        p=0.4720, Theta=0.000, if |C|< 2000 then max|C|=1798
        p=0.4730, Theta=0.0500, if |C|< 2000 then max|C|=1717
        p=0.4740, Theta=0.150, if |C|< 2000 then max|C|=1924
        p=0.4750, Theta=0.150, if |C|< 2000 then max|C|=1893
        p=0.4760, Theta=0.150, if |C|< 2000 then max|C|=1458
        p=0.4770, Theta=0.150, if |C|< 2000 then max|C|=1573
        p=0.4780, Theta=0.200, if |C|< 2000 then max|C|=1762
        p=0.4790, Theta=0.250, if |C|< 2000 then max|C|=951
    """
    print("d = %s, n = number of samples = %s" % (d, n))
    print("stop counting at = %s" % stop)
    for p in range_p:
        p = numerical_approx(p, digits=4)
        S = BondPercolationSamples(p,d,n)
        L = [a for a in S.cluster_cardinality(stop) if not isinstance(a, str)]
        Y = max(L) if L else -Infinity
        theta = S.percolation_probability(stop)
        print("p=%s, Theta=%s, if |C|< %s then max|C|=%s" % (p, theta, stop, Y))


class PercolationProbability(SageObject):
    def __init__(self, d, n, stop, verbose=False):
        r"""
        EXAMPLES::

            sage: from slabbe import PercolationProbability
            sage: f = PercolationProbability(d=2, n=10, stop=100)
            sage: f
            Percolation Probability $\theta(p)$
            d = dimension = 2
            n = # samples = 10
            stop counting at = 100
        """
        self._dimension = d
        self._n = n
        self._stop = stop
        self._verbose = verbose

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import PercolationProbability
            sage: f = PercolationProbability(d=2, n=10, stop=100)
            sage: f
            Percolation Probability $\theta(p)$
            d = dimension = 2
            n = # samples = 10
            stop counting at = 100
        """
        s = "Percolation Probability $\\theta(p)$\n"
        s += "d = dimension = %s\n" % self._dimension
        s += "n = # samples = %s\n" % self._n
        s += "stop counting at = %s" % self._stop
        return s

    @cached_method
    def __call__(self, p ):
        r"""
        EXAMPLES::

            sage: from slabbe import PercolationProbability
            sage: f = PercolationProbability(d=2, n=10, stop=100)
            sage: f
            Percolation Probability $\theta(p)$
            d = dimension = 2
            n = # samples = 10
            stop counting at = 100
            sage: f(0.4534)         # random
            0.300
        """
        if self._verbose:
            print(p)
        S = BondPercolationSamples(p, self._dimension, self._n)
        theta = S.percolation_probability(self._stop)
        return theta

    def return_plot(self, interval=(0,1), adaptive_recursion=4,
            plot_points=4,
            adaptive_tolerance=0.10):
        r"""
        Return a plot of percolation probability using basic sage plot settings.

        INPUT:

        - ``interval``, default=(0,1)
        - ``adaptive_recursion``, default=0
        - ``plot_points``, default=10
        - ``adaptive_tolerance`` default=0.10

        EXAMPLES::

            sage: from slabbe import PercolationProbability
            sage: T = PercolationProbability(d=2, n=10, stop=100)
            sage: T.return_plot()           # optional long
            Graphics object consisting of 1 graphics primitive
        """
        P = plot(self, interval, adaptive_recursion=adaptive_recursion,
                 plot_points=plot_points,
                 adaptive_tolerance=adaptive_tolerance)
        P += text(repr(self), (0.8,0.2))
        return P


