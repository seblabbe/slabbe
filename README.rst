Sébastien Labbé Research Code
=============================

Create the Sage package ``slabbe-0.1.spkg`` from this directory (named with the
version)::

    sage -pkg slabbe-0.1

Installation::

    sage -f slabbe-0.1.spkg

Testing the spkg (make sure the line ``from slabbe import *`` is in the file
``$DOT_SAGE/init.sage``)::

    sage -t --force-lib -i slabbe-0.1/src/slabbe/*

Check the doctest coverage::

    sage -coverage slabbe-0.1/src/slabbe/*

Build the documentation::

    sage -docbuild slabbe inventory
    sage -docbuild slabbe html
    sage -docbuild slabbe pdf

