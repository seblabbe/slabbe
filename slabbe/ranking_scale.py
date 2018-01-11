# coding=utf-8
r"""
Ranking Scale

EXEMPLES::

    sage: from slabbe.ranking_scale import *
    sage: R = RankingScale_CQU4_2011()
    sage: R = RankingScale_USAU_2013()
    sage: R = RankingScale_CQU4_2014()

AUTHOR : 

- Sébastien Labbé, Fall 2011, first version
"""
#*****************************************************************************
#       Copyright (C) 2010-2014 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function

import csv
import itertools
from sage.functions.other import sqrt
from sage.symbolic.constants import e

def RankingScale_CQU4_2011():
    r"""
    EXEMPLES::
        
        sage: from slabbe.ranking_scale import RankingScale_CQU4_2011
        sage: R = RankingScale_CQU4_2011()
        sage: R.table()
          Position   Grand Chelem   Mars Attaque   La Flotte   Petit Chelem
          1          1000           1000           800         400
          2          938            938            728         355
          3          884            884            668         317
          4          835            835            614         285
          5          791            791            566         256
          6          750            750            522         230
          7          711            711            482         206
          8          675            675            444         184
          9          641            641            409         164
          10         609            609            377         146
        ...
          95         0              11             0           0
          96         0              10             0           0
          97         0              9              0           0
          98         0              8              0           0
          99         0              7              0           0
          100        0              6              0           0
          101        0              4              0           0
          102        0              3              0           0
          103        0              2              0           0
          104        0              1              0           0
    """
    from sage.rings.integer_ring import ZZ
    serieA = [0] + discrete_curve(50, 1000, K=1, R=sqrt(2), base=e) #Movember, Bye Bye, Cdf, MA
    la_flotte = [0] + discrete_curve(32, 800, K=1, R=sqrt(2), base=e)  # la flotte
    serieB = [0] + discrete_curve(24, 400, K=1, R=sqrt(2), base=e) # october fest, funenuf, la viree

    nb_ma = 104
    pivot_x = 40
    pivot_y = serieA[pivot_x]
    slope = (1-pivot_y) / (nb_ma-pivot_x)
    L = [pivot_y + slope * (p-pivot_x) for p in range(pivot_x, nb_ma+1)]
    L = [ZZ(round(_)) for _ in L]
    mars_attaque = serieA[:pivot_x] + L

    scales = serieA, mars_attaque, la_flotte, serieB
    names = ["Grand Chelem", "Mars Attaque", "La Flotte", "Petit Chelem"]
    return RankingScale(names, scales)

def RankingScale_USAU_2013():
    r"""
    EXAMPLES::

        sage: from slabbe.ranking_scale import RankingScale_USAU_2013
        sage: R = RankingScale_USAU_2013()
        sage: R.table()
          Position   Serie 1500   Serie 1000   Serie 500   Serie 250
          1          1500         1000         500         250
          2          1366         911          455         228
          3          1252         835          417         209
          4          1152         768          384         192
          5          1061         708          354         177
          6          979          653          326         163
          7          903          602          301         151
          8          832          555          278         139
          9          767          511          256         128
          10         706          471          236         118
          11         648          432          217         109
          12         595          397          199         100
          13         544          363          182         91
          14         497          332          166         83
          15         452          302          151         76
          16         410          274          137         69
          17         371          248          124         62
          18         334          223          112         56
          19         299          199          100         50
          20         266          177          89          45
          21         235          157          79          40
          22         206          137          69          35
          23         178          119          60          30
          24         152          102          51          26
          25         128          86           43          22
          26         106          71           36          18
          27         85           57           29          15
          28         66           44           22          12
          29         47           32           16          9
          30         31           21           11          6
          31         15           10           6           3
          32         1            1            1           1
    """
    L1500 = [0] + discrete_curve(32, 1500, K=1, R=sqrt(2), base=e)
    L1000 = [0] + discrete_curve(32, 1000, K=1, R=sqrt(2), base=e)
    L500 = [0] + discrete_curve(32, 500, K=1, R=sqrt(2), base=e)
    L250 = [0] + discrete_curve(32, 250, K=1, R=sqrt(2), base=e)

    scales = L1500, L1000, L500, L250
    names = ['Serie 1500', 'Serie 1000', 'Serie 500', 'Serie 250']
    return RankingScale(names, scales)

def RankingScale_CQU4_2014(R=1, base=e):
    r"""
    EXAMPLES::

        sage: from slabbe.ranking_scale import RankingScale_CQU4_2014
        sage: R = RankingScale_CQU4_2014()
    """
    L1000 = [0] + discrete_curve(100, 1000, K=1, R=R, base=base)
    L666 = [0] + discrete_curve(50, 666, K=1, R=R, base=base)
    L333 = [0] + discrete_curve(24, 333, K=1, R=R, base=base)

    scales = L1000, L666, L333
    names = ['Serie 1000', 'Serie 666', 'Serie 333']
    return RankingScale(names, scales)


def RankingScale_CQU4_2015_v1():
    r"""
    EXAMPLES::

        sage: from slabbe.ranking_scale import RankingScale_CQU4_2015_v1
        sage: R = RankingScale_CQU4_2015_v1()
    """
    L1000 = [0] + discrete_curve_2(64, 1000)
    L500 = [0] + discrete_curve_2(32, 500)
    L250 = [0] + discrete_curve_2(16, 250)

    scales = L1000, L500, L250
    names = ['Serie 1000', 'Serie 500', 'Serie 250']
    return RankingScale(names, scales)
def RankingScale_CQU4_2015_v2():
    r"""
    EXAMPLES::

        sage: from slabbe.ranking_scale import RankingScale_CQU4_2015_v2
        sage: R = RankingScale_CQU4_2015_v2()
    """
    L1000 = [0] + discrete_curve_2(100, 1000)
    L500 = [0] + discrete_curve_2(40, 500)
    L250 = [0] + discrete_curve_2(12, 250)

    scales = L1000, L500, L250
    names = ['Serie 1000', 'Serie 500', 'Serie 250']
    return RankingScale(names, scales)
def RankingScale_CQU4_2015_v3():
    r"""
    EXAMPLES::

        sage: from slabbe.ranking_scale import RankingScale_CQU4_2015_v3
        sage: R = RankingScale_CQU4_2015_v3()
    """
    L1000 = [0] + discrete_curve_2(100, 1000)
    L666 = [0] + discrete_curve_2(40, 666)
    L333 = [0] + discrete_curve_2(12, 333)

    scales = L1000, L666, L333
    names = ['Serie 1000', 'Serie 666', 'Serie 333']
    return RankingScale(names, scales)


def RankingScale_CQU4_2015_v4():
    r"""
    EXAMPLES::

        sage: from slabbe.ranking_scale import RankingScale_CQU4_2015_v4
        sage: R = RankingScale_CQU4_2015_v4()
    """
    R=1
    base=e
    L1000 = [0] + discrete_curve(64, 1000, K=1, R=R, base=base)
    L666 = [0] + discrete_curve(32, 666, K=1, R=R, base=base)
    L333 = [0] + discrete_curve(16, 333, K=1, R=R, base=base)

    scales = L1000, L666, L333
    names = ['Serie 1000', 'Serie 666', 'Serie 333']
    return RankingScale(names, scales)


def RankingScale_CQU4_2015_v5():
    r"""
    EXAMPLES::

        sage: from slabbe.ranking_scale import RankingScale_CQU4_2015_v5
        sage: R = RankingScale_CQU4_2015_v5()
    """
    R=1
    base=e
    L1000 = [0] + discrete_curve(64, 1000, K=1, R=R, base=base)
    L500 = [0] + discrete_curve(32, 500, K=1, R=R, base=base)
    L250 = [0] + discrete_curve(16, 250, K=1, R=R, base=base)

    scales = L1000, L500, L250
    names = ['Serie 1000', 'Serie 500', 'Serie 250']
    return RankingScale(names, scales)


class RankingScale(object):
    def __init__(self, scale_names, scales):
        r"""
        INPUT:

        - ``scale_names`` -- iterable of str
        - ``scales`` -- iterable of list

        NOTE: the ordering of both input must match
        """
        self._scale_names = tuple(scale_names)
        self._scales = tuple(scales)

    def length(self):
        r"""
        EXEMPLES::

            sage: from slabbe.ranking_scale import RankingScale_CQU4_2011
            sage: R = RankingScale_CQU4_2011()
            sage: R.length()
            105
        """
        return max(map(len, self._scales))

    def pointage_csv(self, filename='pointage.csv', dialect='excel'):
        r"""
        EXEMPLES::

            sage: from slabbe.ranking_scale import RankingScale_CQU4_2011
            sage: R = RankingScale_CQU4_2011()
            sage: R.pointage_csv()          # not tested
            Creation of file pointage.csv
        """
        with open(filename, 'w') as f:
            csv_writer = csv.writer(f, dialect=dialect)
            M = self.length()
            names = ('Position',) + self._scale_names
            scales = (range(M),) + self._scales
            rows.append(names)
            Z = itertools.izip_longest(*scales, fillvalue=0)
            csv_writer.writerow(names)
            Z = itertools.izip_longest(*scales, fillvalue=0)
            Z.next() # enlever la premiere ligne de zeros
            csv_writer.writerows(Z)
            print("Creation of file %s" % filename)

    def table(self):
        r"""
        EXAMPLES::

            sage: from slabbe.ranking_scale import RankingScale_CQU4_2011
            sage: R = RankingScale_CQU4_2011()
            sage: R.table()
              Position   Grand Chelem   Mars Attaque   La Flotte   Petit Chelem
              1          1000           1000           800         400
              2          938            938            728         355
              3          884            884            668         317
              4          835            835            614         285
              5          791            791            566         256
              6          750            750            522         230
              7          711            711            482         206
              8          675            675            444         184
              9          641            641            409         164
              10         609            609            377         146
            ...
              95         0              11             0           0
              96         0              10             0           0
              97         0              9              0           0
              98         0              8              0           0
              99         0              7              0           0
              100        0              6              0           0
              101        0              4              0           0
              102        0              3              0           0
              103        0              2              0           0
              104        0              1              0           0
        """
        from sage.misc.table import table
        rows = []
        M = self.length()
        names = ('Position',) + self._scale_names
        scales = (range(M),) + self._scales
        rows.append(names)
        Z = itertools.izip_longest(*scales, fillvalue=0)
        Z.next() # enlever la premiere ligne de zeros
        rows.extend(Z)
        return table(rows)

    def plot(self, pointsize=20):
        from sage.plot.graphics import Graphics
        G = Graphics()
        m = len(self._scales)
        for i,(name,scale) in enumerate(zip(self._scale_names, self._scales)):
            G += list_plot(scale, color=hue(1.*i/m), legend_label=name,
                    pointsize=pointsize)
        return G

######################
# Fonctions de courbes
######################
def curve(nb_equipes, max_points=100, K=1, R=2, base=2, verbose=False):
    r"""
    INPUT:

    - ``nb_equipes`` -- integer
    - ``max_points`` -- integer
    - ``K`` -- the value at ``p = nb_equipes``
    - ``R`` -- real (default: ``2``), curve parameter
    - ``base`` -- 2
    - ``verbose`` - bool 

    EXEMPLES::

        sage: from slabbe.ranking_scale import curve
        sage: curve(20, 100)
        -99*(p*(log(40) + 1) - p*log(p) - 20*log(40) + 20*log(20) -
        20)/(19*log(40) - 20*log(20) + 19) + 1
        sage: curve(64, 100)
        -33*(p*(7*log(2) + 1) - p*log(p) - 64*log(2) - 64)/(19*log(2) + 21) + 1

    ::

        sage: curve(64, 100)(p=64)
        1
        sage: curve(64, 100)(p=1)
        100
        sage: curve(64, 100)(p=2)
        66*(26*log(2) + 31)/(19*log(2) + 21) + 1
        sage: n(curve(64, 100)(p=2))     # abs tol 1e-10
        95.6871477097753

    ::

        sage: curve(64, 100, verbose=True)
        fn = -(p*(7*log(2) + 1) - p*log(p) - 64*log(2) - 64)/log(2)
        aire = 147.889787576005
        fn normalise = -33*(p*(7*log(2) + 1) - p*log(p) - 64*log(2) - 64)/(19*log(2) + 21) + 1
        -33*(p*(7*log(2) + 1) - p*log(p) - 64*log(2) - 64)/(19*log(2) + 21) + 1

    The base argument seems to be useless (why?)::

        sage: curve(100,100,base=3)
        -99*(p*(log(200) + 1) - p*log(p) - 100*log(200) + 200*log(10) - 
        100)/(99*log(200) - 200*log(10) + 99) + 1
        sage: curve(100,100,base=2)
        -99*(p*(log(200) + 1) - p*log(p) - 100*log(200) + 200*log(10) -
        100)/(99*log(200) - 200*log(10) + 99) + 1
        
    """
    from sage.symbolic.assumptions import forget, assume
    from sage.misc.functional import integrate, n
    from sage.functions.log import log
    from sage.calculus.var import var
    x,p = var('x,p')
    forget()
    assume(p - 1 > 0)
    assume(p-nb_equipes < 0)
    fn = integrate(log(R*nb_equipes, base=base) - log(x, base=base), x, p, nb_equipes)
    if verbose: print("fn = %s" % fn)
    aire = fn(p=1)
    if verbose: print("aire = %s" % n(aire))
    fn_normalise = fn / aire * (max_points - K) + K
    if verbose: print("fn normalise = %s" % fn_normalise)
    return fn_normalise

def discrete_curve(nb_equipes, max_points=100, K=1, R=2, base=2, verbose=False):
    r"""
    INPUT:

    - ``nb_equipes`` -- integer
    - ``max_points`` -- integer
    - ``K`` -- the value at ``p = nb_equipes``
    - ``R`` -- real (default: ``2``), curve parameter
    - ``base`` -- 2
    - ``verbose`` - bool 

    EXEMPLES::
        
        sage: from slabbe.ranking_scale import discrete_curve
        sage: A = discrete_curve(64, 100,verbose=True)
        First difference sequence is
        [1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1,
        2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 2, 3, 3, 3, 4, 4, 4]

    ::

        sage: A = discrete_curve(64, 100)
        sage: B = discrete_curve(32, 50)
        sage: C = discrete_curve(16, 25)

    ::

        sage: A
        [100, 96, 92, 88, 85, 82, 79, 77, 74, 71, 69, 67, 64, 62, 60, 58,
        56, 54, 52, 50, 49, 47, 45, 44, 42, 41, 39, 38, 36, 35, 33, 32, 31,
        29, 28, 27, 26, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12,
        11, 10, 9, 8, 8, 7, 6, 5, 5, 4, 3, 2, 2, 1]
        sage: B
        [50, 46, 43, 40, 37, 35, 33, 31, 29, 27, 25, 23, 21, 20, 18, 17,
        16, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 2, 1]
        sage: C
        [25, 22, 19, 17, 15, 13, 11, 9, 8, 7, 6, 5, 4, 3, 2, 1]

    ::

        sage: A = discrete_curve(64+2, 100) #Movember, Bye Bye, Cdf, MA
        sage: B = discrete_curve(32+2, 70)  # la flotte
        sage: C = discrete_curve(16+2, 40) # october fest, funenuf, la viree
    """
    from sage.misc.functional import round
    from sage.rings.integer_ring import ZZ
    from sage.combinat.words.word import Word
    fn_normalise = curve(nb_equipes, max_points, K=K, R=R, base=base)
    L = [ZZ(round(fn_normalise(p=i))) for i in range(1,nb_equipes+1)]
    if verbose: 
        print("First difference sequence is")
        print(list(Word(L).reversal().finite_differences()))
    return L

def discrete_curve_2(nb_equipes, max_points=100):
    r"""
    """
    from sage.misc.functional import round
    from sage.functions.log import log
    from sage.rings.integer_ring import ZZ
    f = lambda p:(max_points-1)*(1-log(p)/log(nb_equipes))+1
    L = [ZZ(round(f(p=i))) for i in range(1,nb_equipes+1)]
    return L

######################
# Should be in Sage
######################
def table_to_csv(self, filename, dialect='excel'):
    r"""
    """
    with open(filename, 'w') as f:
        csv_writer = csv.writer(f, dialect=dialect)
        csv_writer.writerows(self._rows)
        print("Creation of file %s" % filename)

