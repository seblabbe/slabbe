===============================
Sébastien Labbé's Research Code
===============================

This is a optional package for SageMath containing part of my own research
code. It contains modules on 

**Digital geometry**
  Christoffel graph, discrete subset, discrete plane, double square tiles,
  polyhedron partitions, 

**Combinatorics**
  2d substitutions, bispecial factors, bond percolation, Dyck word in 3D,
  words, Joyal bijection, languages, Oldenburger sequence, ostrowski
  numeration, partial injections,

**Discrete dynamical systems**
  diophantine approximation, Markov transformations, Wang tilings, Lyapunov
  exponents, matrix cocycles, multidimensional continued fraction algorithms,
  polyhedron exchange transformations, 

**Vizualization**
  tikzpicture

**Miscellaneous**
  analyze Sage build time, fruit Python classes example, ranking scale

This package depends heavily on the SageMath__ library as it uses the following
modules: combinat, functions, geometry, graphs, matrix, misc, modules,
numerical, parallel, plot, probability, rings, sat, sets, structure, symbolic.

__ http://www.sagemath.org/

Links: slabbe.org__, `slabbe on Github`__, `slabbe on pypi`__

__ http://www.slabbe.org/Sage/
__ http://github.com/seblabbe/slabbe
__ http://pypi.python.org/pypi/slabbe

Future inclusion into Sage
--------------------------

Some modules may have a wider interest to the SageMath community
(``tikz_picture.py`` for example) and could be included in SageMath at some
point. Please contact the author if you want to act as a reviewer for some
module(s) and I will create a ticket on trac__ for its inclusion into SageMath.

__ https://trac.sagemath.org/

Prerequisites
-------------

Installing slabbe requires a working Sage installation (with Cython and gcc).
It is also recommanded to install the optional SageMath packages dot2tex,
glucose, cryptominisat and latte_int as well as the external packages
Gravphviz and Gurobi.

Installation
------------

The module is distributed on PyPI and is easily installed through the Python
package manager pip::

    sage -pip install slabbe

To install the module in your user space (and does not require administrator
rights)::

    sage -pip install slabbe --user

Usage::

    sage: from slabbe import *

Release history
---------------

**Version 0.5 (April 10, 2019)**
  Few fixes for the version 2 of "Substitutive structure of Jeandel-Rao
  aperiodic tilings". New additions includes solving Wang tilings problem
  using SAT solvers and a class for Polyhedron exchange transformations.

**Version 0.4.4 (September 28, 2018)**
  Make ``import slabbe`` work in Sage with Python 3.

**Version 0.4.3 (August 22, 2018)**
  Few fixes for the publication of "Substitutive structure of Jeandel-Rao
  aperiodic tilings".

**Version 0.4.2 (July 20, 2018)**
  Few fixes for the version 2 of "A self-similar aperiodic set of 19 Wang
  tiles".

**Version 0.4.1 (February 9, 2018)**
  Few fixes for the publication of "A self-similar aperiodic set of 19 Wang
  tiles".  New module to solve the Magic hexagon problem.

**Version 0.4 (January 20, 2018)**
  Version ``0.4`` includes new modules for Wang tilings, 2d substitutions,
  polyhedron partitions, partial injections, ostrowski numeration and many
  improvements to other modules.

**Version 0.3b2 (December 11, 2016)**
  Version ``0.3b2`` includes a new module for diophantine approximations,
  random point generation inside polytopes, analyzing sage building time, and
  many improvements to previous modules.

**Version 0.3b1 (June 12, 2016)**
  Version ``0.3b1`` is now a Python package available in the Python Package
  Index (PyPI). It was migrated from the previous sage optional spkg old-style
  format. It also adds code to deal with bispecial factors, some new methods
  of graphs, substitutions and matrices.

**Version 0.2 (November 25, 2015)**
  slabbe-0.2.spkg__ (documentation__) provides modules on multidimensional
  continued fraction algorithms, matrix cocycles, languages and tikzpictures.  

__ http://www.slabbe.org/Sage/slabbe-0.2.spkg
__ http://www.slabbe.org/Sage/slabbe-0.2.pdf

**Version 0.1.1 (June 3, 2015)**
  slabbe-0.1.1.spkg__ fixes a bug with ``gcd`` import error.

__ http://www.slabbe.org/Sage/slabbe-0.1.1.spkg

**Version 0.1 (August 27, 2014)**
  slabbe-0.1.spkg__ (documentation__) contains modules on digital geometry,
  combinatorics on words and more. 

__ http://www.slabbe.org/Sage/slabbe-0.1.spkg
__ http://www.slabbe.org/Sage/slabbe-0.1.pdf

