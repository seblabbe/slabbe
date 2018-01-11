# -*- coding: utf-8 -*-
r"""
Christoffel graph

This module was developped for the article on a d-dimensional extension of
Christoffel Words written with Christophe Reutenauer [LR2014]_.

.. [LR2014] Labbé, Sébastien, and Christophe Reutenauer. A d-dimensional Extension of
   Christoffel Words. arXiv:1404.4021__ (April 15, 2014).

__ http://arxiv.org/abs/1404.4021

EXAMPLES:

Christoffel graph in 2d (tikz code)::

    sage: from slabbe import ChristoffelGraph, DiscreteBox
    sage: C = ChristoffelGraph((2,5))
    sage: b = DiscreteBox([-5,5],[-5,5])
    sage: I = C & b
    sage: point_kwds = {'label':lambda p:C.level_value(p),'label_pos':'above right'}
    sage: tikz = I.tikz_noprojection(scale=0.8,point_kwds=point_kwds)

Christoffel graph in 3d (tikz code)::

    sage: C = ChristoffelGraph((2,3,5))
    sage: tikz = C.tikz_kernel()

TODO:

    - Clean kernel_vector method of ChristoffelGraph

"""
#*****************************************************************************
#       Copyright (C) 2013-2014 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
import itertools
from copy import copy
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from sage.functions.other import abs
from slabbe.discrete_subset import DiscreteSubset, DiscreteTube
from slabbe.discrete_plane import DiscretePlane 
from slabbe.matrices import M3to2, M4to3
################################################
# Christoffel Graph
################################################
class ChristoffelGraph(DiscreteSubset):
    r"""
    Subset of a discrete object such that its projection by a matrix is
    inside a certain box.

    INPUT:

    - ``v`` - vector, normal vector

    EXAMPLES::

        sage: from slabbe import ChristoffelGraph
        sage: ChristoffelGraph((2,5))
        Christoffel set of edges for normal vector v=(2, 5)

    ::

        sage: C = ChristoffelGraph((2,5))
        sage: it = C.edges_iterator()
        sage: it.next()
        ((0, 0), (1, 0))

    ::

        sage: C = ChristoffelGraph((2,5,8))
        sage: it = C.edges_iterator()
        sage: it.next()
        ((0, 0, 0), (1, 0, 0))

    ::

        sage: from slabbe import DiscreteBox
        sage: C = ChristoffelGraph((2,5))
        sage: b = DiscreteBox([-5,5],[-5,5])
        sage: I = C & b
        sage: point_kwds = {'label':lambda p:C.level_value(p),'label_pos':'above right'}
        sage: tikz = I.tikz_noprojection(scale=0.8,point_kwds=point_kwds)

    TEST:

    This was once a bug. We make sure it is fixed::

        sage: from slabbe import DiscreteSubset
        sage: C = ChristoffelGraph((2,3,5))
        sage: isinstance(C, DiscreteSubset)
        True

    """
    def __init__(self, v, mod=None):
        r"""
        Constructor.

        Return the Christoffel set of edges.

        EXAMPLES::

            sage: from slabbe import ChristoffelGraph
            sage: ChristoffelGraph((2,5))
            Christoffel set of edges for normal vector v=(2, 5)

        """
        self._v = vector(v)
        if mod:
            self._sum_v = mod
        else:
            self._sum_v = sum(abs(a) for a in v)
        DiscreteSubset.__init__(self, dimension=len(v), edge_predicate=self.has_edge)

    def level_value(self, p):
        r"""
        Return the level value of a point p.

        INPUT:

        - ``p`` - point in the space

        EXAMPLES::

            sage: from slabbe import ChristoffelGraph
            sage: C = ChristoffelGraph((2,5,8))
            sage: C.level_value(vector((2,3,4)))
            6
            sage: C.level_value(vector((1,1,1)))
            0
        """
        return p.dot_product(self._v) % self._sum_v

    def has_edge(self, p, s):
        r"""
        Returns whether it has the edge (p, s) where s-p is a canonical
        vector.

        INPUT:

        - ``p`` - point in the space
        - ``s`` - point in the space

        EXAMPLES::

            sage: from slabbe import ChristoffelGraph
            sage: C = ChristoffelGraph((2,5,8))
            sage: C.has_edge(vector((0,0,0)), vector((0,0,1)))
            True
            sage: C.has_edge(vector((0,0,0)), vector((0,0,2)))
            False
            sage: C.has_edge(vector((0,0,0)), vector((0,0,-1)))
            False

        ::

            sage: C = ChristoffelGraph((2,5))
            sage: C.has_edge(vector((0,0)),vector((1,0)))
            True
            sage: C.has_edge(vector((0,0)),vector((-1,0)))
            False
            sage: C.has_edge(vector((-1,1)),vector((1,0)))
            False
        """
        F = self.level_value
        sum_s_minus_p = sum(s-p)
        if sum_s_minus_p == 1:
            return F(p) < F(s)
        elif sum_s_minus_p == -1:
            return F(p) > F(s)
        else:
            return False
            #raise ValueError("vector p=(%s) and s(=%s) should difer by a canonical vector" % (p,s))

    def _repr_(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from slabbe import ChristoffelGraph
            sage: ChristoffelGraph((2,5))
            Christoffel set of edges for normal vector v=(2, 5)
            sage: ChristoffelGraph((2,5,8))
            Christoffel set of edges for normal vector v=(2, 5, 8)
        """
        s = "Christoffel set of edges for normal vector v=%s" % self._v
        return s

    def kernel_vector(self, way='LLL', verbose=False):
        r"""
        todo: clean this

        EXAMPLES::

            sage: from slabbe import ChristoffelGraph
            sage: C = ChristoffelGraph((2,5,7))
            sage: C.kernel_vector()
            [(-1, -1, 1), (3, -4, 0)]

        """
        from sage.arith.misc import gcd
        if way == 'vect_gcd':
            a,b,c = self._v
            gcd_ac = gcd(a,c)
            gcd_bc = gcd(b,c)
            U = ua,ub,uc = vector((c,0,-a)) / gcd(a,c)
            V = va,vb,vc = vector((0,c,-b)) / gcd(b,c)
            rows = U,V
        elif way == 'echelon':
            a,b,c = self._v
            m = matrix(ZZ, 4, [1,1,1,c,0,-a,0,c,-b,b,-a,0])
            me = m.echelon_form()
            if verbose:
                print(me)
            rows = me[1],me[2]
        elif way == 'LLL':
            dim = self.dimension()
            if dim == 3:
                a,b,c = self._v
                M = matrix(ZZ, 4, [1,1,1,c,0,-a,0,c,-b,b,-a,0])
            elif dim == 4:
                a,b,c,d = self._v
                M = matrix(ZZ, 7, (1,1,1,1,b,-a,0,0,c,0,-a,0,0,c,-b,0,d,0,0,-a,0,d,0,-b,0,0,d,-c))
            else:
                raise ValueError("dimension (=%s) must be 3 or 4" % dim)
            rows = M.LLL().rows()
            VS = rows[0].parent()
            zero = VS(0)
            un = VS((1,)*dim)
            assert zero in rows, "(0,0,0) not in LLL result"
            assert un in rows, "(1,1,1) not in LLL result"
            while zero in rows: rows.remove(zero)
            while un in rows: rows.remove(un)
        elif way == 'vect':
            a,b,c = self._v
            U = ua,ub,uc = vector((c,0,-a))
            V = va,vb,vc = vector((0,c,-b))
            rows = U,V
        else:
            raise ValueError("unknown way")
        R = matrix(rows)
        if sum(map(abs, R.minors(dim-1))) != sum(map(abs,self._v)):
            print(R)
            print(R.minors(dim-1))
            print(sum(map(abs, R.minors(dim))))
            print(sum(map(abs,self._v)))
            raise Exception("The result (=%s) is false " % rows)
        return rows

    def tikz_kernel(self, projmat=M3to2, scale=1, edges=True,
            points=True, label=False, point_kwds={}, edge_kwds={}, extra_code='',
            way='LLL', kernel_vector=None):
        r"""
        INPUT:

        - ``projmat`` -- (default: M3to2) 2 x dim projection
          matrix where dim is the dimensoin of self, the isometric
          projection is used by default
        - ``scale`` -- real number (default: 1), scaling constant for the
          whole figure
        - ``edges`` - bool (optional, default: ``True``), whether to draw
          edges
        - ``points`` - bool (optional, default: ``True``), whether to draw
          points
        - ``point_kwds`` - dict (default: ``{}``)
        - ``edge_kwds`` - dict (default: ``{}``)
        - ``extra_code`` -- string (default: ``''``), extra tikz code to add
        - ``way`` -- string (default: ``'LLL'``), the way the base of the
          kernel is computed
        - ``kernel_vector`` -- list (default: ``None``), the vectors, if
          None it uses ``kernel_vector()`` output.

        EXAMPLES::

            sage: from slabbe import ChristoffelGraph
            sage: C = ChristoffelGraph((2,3,5))
            sage: tikz = C.tikz_kernel()
            sage: tikz
            \documentclass[tikz]{standalone}
            \usepackage{amsmath}
            \begin{document}
            \begin{tikzpicture}
            [scale=1]
            \clip (-5.239453693, 2.025) -- (-5.239453693, -1.025) -- (-4.330127019, -1.55) -- (0.0, -1.05) -- (0.909326674, -0.525) -- (0.909326674, 2.525) -- (0.0, 3.05) -- (-4.330127019, 2.55) -- cycle;
             \draw[very thick, blue] (0.00000, 0.00000) -- (-0.86603, -0.50000);
            \draw[very thick, blue] (0.00000, 0.00000) -- (0.86603, -0.50000);
            ...
            ... 296 lines not printed (24077 characters in total) ...
            ...
            \node[circle,fill=black,draw=black,minimum size=0.8mm,inner sep=0pt,] at (-6.92820, -3.00000) {};
             \filldraw[fill=white,very thick,dotted,opacity=0.5,even odd rule]
              (-5.239453693, 2.025) -- (-5.239453693, -1.025) -- (-4.330127019, -1.55) -- (0.0, -1.05) -- (0.909326674, -0.525) -- (0.909326674, 2.525) -- (0.0, 3.05) -- (-4.330127019, 2.55) -- cycle
              (-4.330127019, 1.5) -- (-4.330127019, -0.5) -- (0.0, 0.0) -- (0.0, 2.0) -- cycle;
             \end{tikzpicture}
            \end{document}
        """
        # kernel to clip
        if kernel_vector:
            K = kernel_vector
        else:
            K = self.kernel_vector(way=way)

        if self.dimension() == 4:
            projmat_tube = M4to3
        elif self.dimension() == 3:
            projmat_tube = M3to2
        else:
            raise NotImplementedError("dim must be 3 or 4")

        zero = self.an_element().parent(0)
        para_base_edges = [(zero, vector(U)) for U in K]
        para = [sum(p) for p in itertools.product(*para_base_edges)]
        base_edges = self.base_edges()
        hex = base_edges + [-v for v in base_edges]
        clip_ambiant_dim = [p+1.05*h for p in para for h in hex]
        clip_proj_dim = [projmat_tube * c for c in clip_ambiant_dim]

        Z = zip(*clip_proj_dim)
        minZ = vector(RR, map(min, Z)) - vector(RR, (2,)*len(Z))
        maxZ = vector(RR, map(max, Z)) + vector(RR, (2,)*len(Z))
        tube = DiscreteTube(*zip(minZ,maxZ), projmat=projmat_tube)
        plane = DiscretePlane(self._v, sum(self._v))
        I = self & tube & plane
        #I = Intersection3d(I)

        if self.dimension() == 4:
            projmat = M3to2*M4to3
        elif self.dimension() == 3:
            projmat = M3to2
        else:
            raise NotImplementedError("dim must be 3 or 4")
        tikz = I.tikz(projmat=projmat, scale=scale, edges=edges,
                points=points, point_kwds=point_kwds, edge_kwds=edge_kwds,
                clip=clip_ambiant_dim, contour=para, extra_code=extra_code)

        return tikz


