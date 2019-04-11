
====
TODO
====

Use gitlab
----------

https://about.gitlab.com/

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

Semantic Versioning
-------------------

https://semver.org/

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

Conseils de Frédéric Chapoton pour le passage à Python 3
--------------------------------------------------------

---------- Forwarded message ---------
From: Frédéric Chapoton
Date: mer. 19 sept. 2018 à 11:24
Subject: Re: Passage à Python3
To: Samuel Lelièvre


Salut Samuel,

une réponse rapide. Je vois plusieurs façons de tester la compatibilité
d'un module :

* utiliser les version python3 (pip3) de pyflakes et de pycodestyle
(anciennement pep8), qui sont capables de trouver les erreurs de
syntaxe, dans une certaine mesure ;

* tout simplement essayer de charger le module dans sage-sous-python3 ;

* utiliser la technique simple proposée dans le ticket

https://trac.sagemath.org/ticket/15995#comment:2

puis éventuellement le script (read_deprecation-warnings-v2.py) en
attachement dans ce ticket pour extraire l'information si besoin.

* faire une branche git avec le module et faire tourner un patchbot sur
cette branche (avec juste les plugins --plugin-only, c'est pas trop
long). Le patchbot n'est hélas pas prévu pour tourner sur autre chose
que sage en entier...

Voila donc deja quelques idees, plus ou moins bonnes. Autrement, je me
sers souvent de "git grep" pour trouver les copies de problemes connus.
Mais je n'ai pas automatisé le processus. J'ai par contre essayé de
mettre des vérifications dans le patchbot.

Frédéric
