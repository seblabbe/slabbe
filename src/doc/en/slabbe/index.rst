==============================================
Sébastien Labbé Research Code Reference Manual
==============================================

This is the reference manual for the `Sébastien Labbé
<http://www.liafa.univ-paris-diderot.fr/~labbe>`_ Research Code extension to
the `Sage <http://sagemath.org/>`_ mathematical software system.  Sage is free
open source math software that supports research and teaching in algebra,
geometry, number theory, cryptography, and related areas.  
Sébastien Labbé Research Code implements digital geometry, combinatorics on
words and symbolic dynamical systems simulation code in Sage, via a set of new
Python classes. 
Many of the modules corresponds to research code written for published articles
(double square tiles, Christoffel graphs, factor complexity).
It is meant to be reused and reusable (full documentation including doctests)
but still has not be used for other purposes. Comments are welcome.

.. [BBL2012]_, [BBGL2011]_, [BGL2012]_

To use this module, you need to import it:: 

    from slabbe import *

.. WARNING:: 

    The above line is **mandatory** for any doctest of this reference
    manual to work.

This reference manual contains many examples that illustrate the usage of
slabbe spkg. The examples are all tested with each release of slabbe spkg, and
should produce exactly the same output as in this manual, except for line
breaks.

This work is licensed under a `Creative Commons Attribution-Share Alike
3.0 License`__.

__ https://creativecommons.org/licenses/by-sa/3.0/

Digital Geometry
================

.. toctree::
   :maxdepth: 1

   discrete_subset
   discrete_plane
   billiard
   christoffel_graph
   double_square_tile

Combinatorics on words
======================

.. toctree::
   :maxdepth: 1

   kolakoski_word
   bispecial_extension_type

Combinatorics
=============
.. toctree::
   :maxdepth: 1

   joyal_bijection
   bond_percolation

Python class inheritance
========================

.. toctree::
   :maxdepth: 1

   fruit

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
