
====
TODO
====

Transition à Python 3
---------------------

 - Ajouter des import print_chose pour prévenir le retour des vieilles méthodes

Transition spkg -> pip
----------------------

Make my spkg a Python package?

 - La documentation (?)
 - Est-ce que l'option -e marche bien? et pour les .pyx?
 - Uploader tout ça sur PyPI
 - Ajouter MANIFEST.in  (see metadata link below)
 - tests or sage tests ?
 - Move doc to docs (conf.py)

https://python-packaging.readthedocs.org/en/latest/
https://python-packaging.readthedocs.org/en/latest/minimal.html
https://python-packaging.readthedocs.org/en/latest/metadata.html
https://python-packaging.readthedocs.org/en/latest/testing.html
https://python-packaging.readthedocs.org/en/latest/non-code-files.html (like doc)

http://stackoverflow.com/questions/4505747/how-should-i-structure-a-python-package-that-contains-cython-code

Jeroen link:
https://packaging.python.org/en/latest/installing/
https://packaging.python.org/en/latest/distributing/

What to do with the documentation of an optional spkg?
https://groups.google.com/forum/#!topic/sage-devel/RKzN-fxIfHI

https://python-packaging-user-guide.readthedocs.io/en/latest/index.html

Other packages to compare with
------------------------------

https://github.com/pypa/sampleproject
https://github.com/williamstein/sage_modabvar
https://github.com/abelfunctions/abelfunctions 
https://github.com/nthiery/sage-semigroups/
https://github.com/sagemath/cysignals

Fichiers à ajouter
------------------

TODO ajouter les fichiers suivants::

    10 slabbe@pol ~/Blogue/Files $ find . -name *.sage
    ./2009/tamer.sage
    ./2012/bond_percolation.sage
    ./2012/joyal_bijection.sage
    ./2014/brigitte.sage
    10 slabbe@pol ~/Blogue/Files $ find . -name *.py
    ./2011/coverage_evolution.py
    ./2013/fruit.py
    ./2014/fichier.py

from my old hg patches::

    # Move to slabbe.spkg
    suffix_trie_gen-sl.patch #+experimental
    suffix_trie_new_classes-sl.patch #+conflict
    wordmorphism-is-pisot-methods-sl.patch #+conflict
    automatic_sequences.patch
    trac13346_kolakoski_cython-sl.patch
    trac_12697-string_process-fh.patch
    boyer_more-sl.patch
    tachyon-sl.patch
    word_infinite_to_finite-sl.patch #+conflict
    word_in_cython-sl.patch #+conflict #with 8418

NOTES
-----

:class:`~sage.combinat.words.morphism.WordMorphism` : sans le chemin au long
:class:`sage.combinat.words.morphism.WordMorphism`

Also works::

    tar -cjf slabbe-0.1.spkg slabbe-0.1

