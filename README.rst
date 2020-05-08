===============================
Sébastien Labbé's Research Code
===============================

README
======

This is an optional package for SageMath containing code I wrote for research.
It contains modules on the following topics.

**Discrete dynamical systems**
  diophantine approximation, Markov transformations, Wang tilings, Lyapunov
  exponents, matrix cocycles, multidimensional continued fraction algorithms,
  polyhedron exchange transformations, 

**Combinatorics**
  2d substitutions, bispecial factors, bond percolation, Dyck word in 3D,
  words, Joyal bijection, languages, Oldenburger sequence, ostrowski
  numeration, partial injections,

**Digital geometry**
  Christoffel graph, discrete subset, discrete plane, double square tiles,
  polyhedron partitions, 

**Vizualization**
  tikzpicture

**Miscellaneous**
  analyze Sage build time, fruit Python classes example, ranking scale

Links: 

 - documentation: http://www.slabbe.org/docs/
 - PyPI: http://pypi.python.org/pypi/slabbe
 - github: http://github.com/seblabbe/slabbe
 - www: http://www.slabbe.org/Sage/

Prerequisites
-------------

Installing slabbe requires a working SageMath installation (with Cython and
gcc). Depending on the usage, it might be necessary to install the optional
packages dot2tex__, glucose__, cryptominisat__ and latte_int__::

    sage -i dot2tex glucose cryptominisat latte_int

__ https://dot2tex.readthedocs.io/en/latest/
__ https://www.labri.fr/perso/lsimon/glucose/
__ https://www.msoos.org/cryptominisat5/
__ https://www.math.ucdavis.edu/~latte/

as well as the external packages ImageMagick__, Graphviz__ and Gurobi__. See
this `thematic tutorial`__ to use Gurobi in SageMath.

On Debian or Ubuntu, one may do::

    sudo apt install graphviz imagemagick

On OSX, one may do after installing Homebrew__::

    sudo brew install graphviz imagemagick

Note that graphviz must be installed *before* dot2tex.

__ https://imagemagick.org/
__ https://graphviz.org/
__ http://www.gurobi.com/
__ http://doc.sagemath.org/html/en/thematic_tutorials/linear_programming.html#using-cplex-or-gurobi-through-sage
__ https://brew.sh/

Installation
------------

The module is distributed on PyPI and is easily installed through the Python
package manager pip::

    sage -pip install slabbe

To install the module in your user space (which does not require administrator
rights)::

    sage -pip install slabbe --user

To install the most recent development version::

    sage -pip install --upgrade git+https://github.com/seblabbe/slabbe

Usage::

    sage: from slabbe import *

It builds on SageMath
---------------------

It depends heavily on the SageMath library as it uses the following modules:
combinat, functions, geometry, graphs, matrix, misc, modules, numerical,
parallel, plot, probability, rings, sat, sets, structure, symbolic.

SageMath__ is free open source math software that supports research and
teaching in algebra, geometry, number theory, cryptography, and related areas.  

__ http://www.sagemath.org/

Follows the Best practices for scientific computing
---------------------------------------------------

It follows as much as possible the `SageMath general conventions`__ and the
`Best Practices for Scientific Computing`__. Each module is fully documented
and doctested. Before each new release, we make sure that all examples are
working. As the `ReScience Journal`__ says: "*Reproducible science is good.
Replicated Science is better*".

__ http://doc.sagemath.org/html/en/developer/coding_basics.html
__ https://doi.org/10.1371/journal.pbio.1001745
__ http://rescience.github.io/

Future inclusion into Sage
--------------------------

Some modules may have a wider interest to the SageMath community
(``tikz_picture.py`` for example) and could be included in SageMath at some
point. Please contact the author if you want to act as a reviewer for some
module(s) and I will create a ticket on trac__ for its inclusion into SageMath.

__ https://trac.sagemath.org/

Release history
---------------

*Version 0.6.1 (May 8, 2020)*
  New modules to deal with the coding of `Z^d`-action by PETs, `d`-dimensional
  sturmian configurations. Improved the computation of induced polyhedron partition
  and induced polyhedron exchange transformation. New modules containing the
  code for the articles `arxiv:1903.06137`__ and `arXiv:1906.01104`__

__ https://arxiv.org/abs/1903.06137
__ https://arxiv.org/abs/1906.01104

*Version 0.6 (November 22, 2019)*
  Make the package work with Python 3. Most of the tests pass with Python 3 now.

*Version 0.5.1 (May 30, 2019)*
  Few fixes for the publication of "Induction of `Z^2`-actions on partitions of
  the 2-torus". Improved html documentation.

*Version 0.5 (April 10, 2019)*
  Few fixes for the version 2 of "Substitutive structure of Jeandel-Rao
  aperiodic tilings". New additions includes solving Wang tilings problem
  using SAT solvers and a class for Polyhedron exchange transformations.

*Version 0.4.4 (September 28, 2018)*
  Make ``import slabbe`` work in Sage with Python 3.

*Version 0.4.3 (August 22, 2018)*
  Few fixes for the publication of "Substitutive structure of Jeandel-Rao
  aperiodic tilings".

*Version 0.4.2 (July 20, 2018)*
  Few fixes for the version 2 of "A self-similar aperiodic set of 19 Wang
  tiles".

*Version 0.4.1 (February 9, 2018)*
  Few fixes for the publication of "A self-similar aperiodic set of 19 Wang
  tiles".  New module to solve the Magic hexagon problem.

*Version 0.4 (January 20, 2018)*
  Version ``0.4`` includes new modules for Wang tilings, 2d substitutions,
  polyhedron partitions, partial injections, ostrowski numeration and many
  improvements to other modules.

*Version 0.3b2 (December 11, 2016)*
  Version ``0.3b2`` includes a new module for diophantine approximations,
  random point generation inside polytopes, analyzing sage building time, and
  many improvements to previous modules.

*Version 0.3b1 (June 12, 2016)*
  Version ``0.3b1`` is now a Python package available in the Python Package
  Index (PyPI). It was migrated from the previous sage optional spkg old-style
  format. It also adds code to deal with bispecial factors, some new methods
  of graphs, substitutions and matrices.

*Version 0.2 (November 25, 2015)*
  slabbe-0.2.spkg__ (documentation__) provides modules on multidimensional
  continued fraction algorithms, matrix cocycles, languages and tikzpictures.  

__ http://www.slabbe.org/Sage/slabbe-0.2.spkg
__ http://www.slabbe.org/Sage/slabbe-0.2.pdf

*Version 0.1.1 (June 3, 2015)*
  slabbe-0.1.1.spkg__ fixes a bug with ``gcd`` import error.

__ http://www.slabbe.org/Sage/slabbe-0.1.1.spkg

*Version 0.1 (August 27, 2014)*
  slabbe-0.1.spkg__ (documentation__) contains modules on digital geometry,
  combinatorics on words and more. 

__ http://www.slabbe.org/Sage/slabbe-0.1.spkg
__ http://www.slabbe.org/Sage/slabbe-0.1.pdf

