=============================
Sébastien Labbé Research Code
=============================

This is the reference manual for the `Sébastien Labbé
<http://www.slabbe.org/Sage>`_'s Research Code extension to the `SageMath
<http://sagemath.org/>`_ mathematical software system. SageMath is free open
source math software that supports research and teaching in algebra, geometry,
number theory, cryptography, and related areas.  

This package implements digital geometry, combinatorics on words and symbolic
dynamical systems simulation code in Sage, via a set of new Python classes.
Many of the modules corresponds to research code written for published articles
(double square tiles, Christoffel graphs, factor complexity).
It is meant to be reused and reusable (full documentation including doctests).
Comments are welcome.

Installation::

    sage -pip install slabbe

To use this module, you need to import it:: 

    from slabbe import *

This reference manual contains many examples that illustrate the usage of
slabbe spkg. The examples are all tested with each release of slabbe spkg, and
should produce exactly the same output as in this manual.

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
   finite_word
   infinite_word
   language
   word_morphisms
   ostrowski

Combinatorics
=============
.. toctree::
   :maxdepth: 1

   joyal_bijection
   bond_percolation
   dyck_3d
   combinat
   graph
   partial_injection

Dynamical systems
=================

.. toctree::
   :maxdepth: 1

   matrix_cocycle
   mult_cont_frac
   diophantine_approx
   lyapunov
   markov_transformation
   matrices
   polyhedron_partition
   substitution_2d
   wang_tiles

Miscellaneous
=============

.. toctree::
   :maxdepth: 1

   analyze_sage_build
   tikz_picture
   ranking_scale
   fruit

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
