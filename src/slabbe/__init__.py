from sage.misc.latex import latex
latex.add_to_preamble('\\usepackage{tikz}')
latex.add_to_preamble('\\usepackage{pgfplots}')
latex.add_to_preamble('\\usetikzlibrary{pgfplots.groupplots}')

from discrete_subset import DiscreteSubset, DiscreteBox, DiscreteTube, Intersection
from billiard import BilliardCube
from discrete_plane import DiscretePlane, DiscreteLine, DiscreteHyperplane
from christoffel_graph import ChristoffelGraph
from bispecial_extension_type import ExtensionType, ExtensionType1to1, ExtensionTypeLong
from double_square_tile import DoubleSquare, christoffel_tile
from fruit import Fruit, Banana, Strawberry
from joyal_bijection import Endofunctions, Endofunction, DoubleRootedTree
from bond_percolation import (BondPercolationSamples, 
                             BondPercolationSample, 
                             PercolationProbability)
from kolakoski_word import KolakoskiWord
from tikz_picture import TikzPicture

# for doctext to work, we import convex_boundary
# from discrete_subset import convex_boundary

# do not import module names just the above stuff
#__all__ = []
