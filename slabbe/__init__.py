from __future__ import absolute_import, print_function

# Eventually, toward support of Python 3
#from __future__ import division, absolute_import, print_function, unicode_literals

# For sphinx to work, we first need to import the sage library
from sage.all_cmdline import *

from sage.misc.latex import latex
latex.add_to_preamble('\\usepackage{tikz}')
latex.add_to_preamble('\\usepackage{pgfplots}')
latex.add_to_preamble('\\usetikzlibrary{pgfplots.groupplots}')

from .discrete_subset import DiscreteSubset, DiscreteBox, DiscreteTube, Intersection
from .billiard import BilliardCube
from .discrete_plane import DiscretePlane, DiscreteLine, DiscreteHyperplane
from .christoffel_graph import ChristoffelGraph
from .bispecial_extension_type import ExtensionType, ExtensionType1to1, ExtensionTypeLong
from .double_square_tile import DoubleSquare, christoffel_tile
from .fruit import Fruit, Banana, Strawberry
from .joyal_bijection import Endofunctions, Endofunction, DoubleRootedTree
from .bond_percolation import (BondPercolationSamples, 
                             BondPercolationSample, 
                             PercolationProbability)
from .tikz_picture import TikzPicture

from .matrices import M3to2, M2to3, M4to2, M4to3

from .substitution_2d import Substitution2d
from .wang_tiles import WangTiling, WangTileSolver, WangTileSet, HexagonalWangTileSet
from .polyhedron_partition import PolyhedronPartition, rotation_mod


from .partial_injection import (number_of_partial_injection,
                                random_partial_injection,
                                random_cyclically_reduced_stallings_graph)

# BUG (sometimes, cython code does not work properly)
try:
    from .kolakoski_word import KolakoskiWord
except ImportError:
    print("There was an error while importing KolakoskiWord cython module")
    print("Maybe that's because the pwd contains the slabbe folder?")
    print("We ignore the problem but this module won't work.")

# for doctext to work, we import convex_boundary
# from .discrete_subset import convex_boundary

# do not import module names just the above stuff
#__all__ = []
