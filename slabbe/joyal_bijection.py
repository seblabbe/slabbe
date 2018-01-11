# -*- coding: utf-8 -*-
r"""
André Joyal's Bijection

Problem suggested by Doron Zeilberger during a talk done at CRM, Montreal,
May 11th, 2012 to compare code in different languages. This is a
implementation of the `Joyal's Bijection`__ using Sage. It will not win for
the most brief code, but it is object oriented, documented, reusable,
testable and allows introspection.

__ http://fr.wikipedia.org/wiki/Bijection_de_Joyal

AUTHOR:

    - Sébastien Labbé, May 12th 2012

TODO:

 - Base Endofunction class on sage's FiniteSetMap classes (for both element
   and parent)

EXAMPLES:

Creation of an endofunction
---------------------------

::

    sage: from slabbe import Endofunction
    sage: L = [7, 0, 6, 1, 4, 7, 2, 1, 5] 
    sage: f = Endofunction(L)
    sage: f
    Endofunction:
    [0..8] -> [7, 0, 6, 1, 4, 7, 2, 1, 5] 

Creation of a double rooted tree
--------------------------------

::

    sage: from slabbe import DoubleRootedTree
    sage: L = [(0,6),(2,1),(3,1),(4,2),(5,7),(6,4),(7,0),(8,5)]
    sage: D = DoubleRootedTree(L, 1, 7)
    sage: D
    Double rooted tree:
    Edges: [(0, 6), (2, 1), (3, 1), (4, 2), (5, 7), (6, 4), (7, 0), (8, 5)]
    RootA: 1
    RootB: 7

Joyal's bijection
-----------------

From the endofunction ``f``, we get a double rooted tree::

    sage: f.to_double_rooted_tree()
    Double rooted tree:
    Edges: [(0, 6), (2, 1), (3, 1), (4, 2), (5, 7), (6, 4), (7, 0), (8, 5)]
    RootA: 1
    RootB: 7

From the double rooted tree ``D``, we get an endofunction::

    sage: D.to_endofunction()
    Endofunction:
    [0..8] -> [7, 0, 6, 1, 4, 7, 2, 1, 5]

In fact, we got ``D`` from ``f`` and vice versa::

    sage: D == f.to_double_rooted_tree()
    True
    sage: f == D.to_endofunction()
    True

Endofunctions are defined on the set [0, 1, ..., n-1]
-----------------------------------------------------

As of now, the code supports only endofunctions defined on the set 
[0, 1, ..., n-1] ::

    sage: L = [1, 0, 3, 4, 5, 7, 1]
    sage: f = Endofunction(L)
    Traceback (most recent call last):
    ...
    ValueError: images of [0..6] must be 0 <= i < 7

Another example
---------------

From a list ``L``, we create an endofunction ``f`` ::

    sage: L = [12, 7, 8, 3, 3, 11, 11, 9, 5, 12, 0, 10, 9]
    sage: f = Endofunction(L)
    sage: f
    Endofunction:
    [0..12] -> [12, 7, 8, 3, 3, 11, 11, 9, 5, 12, 0, 10, 9]

From ``f``, we create a double rooted tree ``D``::

    sage: D = f.to_double_rooted_tree(); D
    Double rooted tree:
    Edges: [(0, 12), (1, 7), (2, 8), (3, 12), (4, 3), (5, 11), 
    (6, 11), (7, 9), (8, 5), (10, 0), (11, 10), (12, 9)]
    RootA: 9
    RootB: 3

And from ``D``, we create an endofunction::

    sage: D.to_endofunction()
    Endofunction:
    [0..12] -> [12, 7, 8, 3, 3, 11, 11, 9, 5, 12, 0, 10, 9]

We test that we recover the initial endofunction ``f``::

    sage: f == f.to_double_rooted_tree().to_endofunction()
    True

A random example
----------------

We define the set of all endofunctions on [0..7]::

    sage: from slabbe import Endofunctions
    sage: E = Endofunctions(8)
    sage: E
    Endofunctions of [0..7]

We choose a random endofunction on the set [0..7]::

    sage: f = E.random_element()
    sage: f                               # random
    Endofunction:
    [0..7] -> [5, 5, 0, 4, 5, 0, 1, 1]

We construct a double rooted tree from it::

    sage: f.to_double_rooted_tree()       # random
    Double rooted tree:
    Edges: [(1, 5), (2, 0), (3, 4), (4, 5), (5, 0), (6, 1), (7, 1)]
    RootA: 0
    RootB: 5

We recover an endofunction from the double rooted tree::

    sage: f.to_double_rooted_tree().to_endofunction()   # random
    Endofunction:
    [0..7] -> [5, 5, 0, 4, 5, 0, 1, 1]

Finally, we check the bijection::

    sage: f == f.to_double_rooted_tree().to_endofunction()
    True

Large random example
--------------------

::

    sage: E = Endofunctions(1000)
    sage: f = E.random_element()
    sage: f == f.to_double_rooted_tree().to_endofunction()
    True

TESTS:

We test the limit cases::

    sage: f = Endofunction([0])
    sage: f == f.to_double_rooted_tree().to_endofunction()
    True
    sage: f = Endofunction([0,1])
    sage: f == f.to_double_rooted_tree().to_endofunction()
    True
    sage: f = Endofunction([1,0])
    sage: f == f.to_double_rooted_tree().to_endofunction()
    True

More extensively::

    sage: E = Endofunctions(1)
    sage: all(f == f.to_double_rooted_tree().to_endofunction() for f in E)
    True
    sage: E = Endofunctions(2)
    sage: all(f == f.to_double_rooted_tree().to_endofunction() for f in E)
    True
    sage: E = Endofunctions(3)
    sage: all(f == f.to_double_rooted_tree().to_endofunction() for f in E)
    True
    sage: E = Endofunctions(4)
    sage: all(f == f.to_double_rooted_tree().to_endofunction() for f in E)
    True

TIMING TESTS:

When the extension of the file is .sage::

    sage: E = Endofunctions(3)
    sage: time all(f == f.to_double_rooted_tree().to_endofunction() for f in E) # not tested
    True
    Time: CPU 0.02 s, Wall: 0.02 s
    sage: E = Endofunctions(4)
    sage: time all(f == f.to_double_rooted_tree().to_endofunction() for f in E) # not tested
    True
    Time: CPU 0.22 s, Wall: 0.22 s
    sage: E = Endofunctions(5)
    sage: time all(f == f.to_double_rooted_tree().to_endofunction() for f in E) # not tested
    True
    Time: CPU 2.82 s, Wall: 2.82 s
    sage: E = Endofunctions(6)
    sage: time all(f == f.to_double_rooted_tree().to_endofunction() for f in E) # not tested
    True
    Time: CPU 45.66 s, Wall: 45.74 s

When the extension of the file is .spyx::

    sage: E = Endofunctions(3)
    sage: time all(f == f.to_double_rooted_tree().to_endofunction() for f in E) # not tested
    True
    Time: CPU 0.02 s, Wall: 0.02 s
    sage: E = Endofunctions(4)
    sage: time all(f == f.to_double_rooted_tree().to_endofunction() for f in E) # not tested
    True
    Time: CPU 0.21 s, Wall: 0.21 s
    sage: E = Endofunctions(5)
    sage: time all(f == f.to_double_rooted_tree().to_endofunction() for f in E) # not tested
    True
    Time: CPU 2.71 s, Wall: 2.72 s
    sage: E = Endofunctions(6)
    sage: time all(f == f.to_double_rooted_tree().to_endofunction() for f in E) # not tested
    True
    Time: CPU 44.08 s, Wall: 44.17 s

When the extension of the file is .sage::

    sage: E = Endofunctions(1000)
    sage: f = E.random_element()
    sage: time f == f.to_double_rooted_tree().to_endofunction() # not tested
    True
    Time: CPU 0.09 s, Wall: 0.09 s
    sage: E = Endofunctions(10000)
    sage: f = E.random_element()
    sage: time f == f.to_double_rooted_tree().to_endofunction()  # not tested
    True
    Time: CPU 2.23 s, Wall: 2.24 s

When the extension of the file is .spyx::

    sage: E = Endofunctions(1000)
    sage: f = E.random_element()
    sage: time f == f.to_double_rooted_tree().to_endofunction() # not tested
    True
    Time: CPU 0.11 s, Wall: 0.11 s
    sage: E = Endofunctions(10000)
    sage: f = E.random_element()
    sage: time f == f.to_double_rooted_tree().to_endofunction() # not tested
    True
    Time: CPU 2.91 s, Wall: 2.93 s

"""
#*****************************************************************************
#       Copyright (C) 2008 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
from copy import copy
from sage.misc.flatten import flatten
from sage.graphs.graph import Graph
from sage.misc.prandom import randint
from sage.sets.finite_set_maps import FiniteSetMaps
from sage.combinat.words.word import Word


class Endofunctions(object):
    r"""
    Returns the set of all endofunction on the set [0..n-1].

    INPUT:

    - ``n`` - positive integer

    EXAMPLES::

        sage: from slabbe import Endofunctions
        sage: Endofunctions(10)
        Endofunctions of [0..9]
    """
    def __init__(self, n):
        r"""
        EXAMPLES::

            sage: from slabbe import Endofunctions
            sage: Endofunctions(10)
            Endofunctions of [0..9]
        """
        self._n = n

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import Endofunctions
            sage: Endofunctions(1)
            Endofunctions of [0..0]
            sage: Endofunctions(10)
            Endofunctions of [0..9]
            sage: Endofunctions(0)
            Endofunctions of [0..-1]
        """
        return "Endofunctions of [0..%s]" % (self._n-1)

    def __iter__(self):
        for f in FiniteSetMaps(self._n):
            yield Endofunction(f)

    def random_element(self):
        r"""
        Return a random endofunction on [0..n-1].

        EXAMPLES::

            sage: from slabbe import Endofunctions
            sage: E = Endofunctions(10)
            sage: E.random_element()          # random
            Endofunction:
            [0..9] -> [2, 8, 7, 0, 0, 6, 2, 3, 5, 9]
            sage: E.random_element()          # random
            Endofunction:
            [0..9] -> [8, 7, 7, 5, 4, 1, 0, 3, 8, 6]
        """
        L = [randint(0,self._n-1) for a in range(self._n)]
        return Endofunction(L)


class Endofunction(object):
    r"""
    Returns an endofunction.

    INPUT:

    - ``L`` - list of length n containing images of the integers from 0 to n-1
      where the images belong to the integers from 0 to n-1.

    EXAMPLES::

        sage: from slabbe import Endofunction
        sage: L = [0, 5, 7, 2, 1, 1, 2, 6, 2, 4]
        sage: f = Endofunction(L)
        sage: f
        Endofunction:
        [0..9] -> [0, 5, 7, 2, 1, 1, 2, 6, 2, 4]
    """
    def __init__(self, L):
        r"""
        """
        if not isinstance(L, list):
            L = list(L)
        len_L = len(L)
        if not all(0 <= i < len_L for i in L):
            raise ValueError("images of [0..%s] must be 0 <= i < %s" % (len_L-1,len_L) )
        self._list = L

    def __repr__(self):
        r"""
        """
        s = "Endofunction:\n"
        s += "[0..%s] -> %s" % (len(self)-1, self._list)
        return s

    def __len__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import Endofunction
            sage: L = [1, 6, 5, 1, 6, 5, 0]
            sage: f = Endofunction(L)
            sage: len(f)
            7
        """
        return len(self._list)

    def __call__(self, a):
        r"""
        EXAMPLES::

            sage: from slabbe import Endofunction
            sage: L = [1, 6, 5, 1, 6, 5, 0]
            sage: f = Endofunction(L)
            sage: f(2)
            5
        """
        return self._list[a]

    def __eq__(self, other):
        r"""
        EXAMPLES::

            sage: from slabbe import Endofunction
            sage: L = [1, 6, 5, 1, 6, 5, 0]
            sage: f = Endofunction(L)
            sage: f == f
            True
            sage: f == 2
            False
        """
        return (isinstance(other, Endofunction) and
                self._list == other._list)

    def __ne__(self, other):
        return not self.__eq__(other)

    def two_cycle_elements(self):
        r"""
        Iterator over elements in a two-cycle.

        EXAMPLES::

            sage: from slabbe import Endofunction
            sage: L = [0, 5, 7, 2, 1, 1, 2, 6, 2, 4]
            sage: f = Endofunction(L)
            sage: list(f.two_cycle_elements())
            [1, 5]
        """
        for i,a in enumerate(self._list):
            if i != a and i == self._list[a] :
                yield i

    def cycle_elements(self):
        r"""
        Returns the list of all elements in a cycle for this endofunction.

        OUTPUT:

        list

        EXAMPLES::

            sage: from slabbe import Endofunction
            sage: L = [6, 5, 7, 2, 1, 1, 2, 6, 2, 4]
            sage: f = Endofunction(L)
            sage: f.cycle_elements()   # random order
            [0, 6, 7, 2, 1, 5]

        .. NOTE::

            ``G.cycle_basis()`` is not implemented for directed or
            multiedge graphs (in Networkx).  Hence, the ``cycle_basis``
            method is missing the 2-cycles.
        """
        G = Graph(loops=True)
        for i, a in enumerate(self._list):
            G.add_edge(i,a)
        C = G.cycle_basis()
        L = flatten(C)
        # adding the missing 2-cycle elements
        for i in self.two_cycle_elements():
            L.append(i)
        return L

    def skeleton(self):
        r"""
        Return the skeleton of the endofunction.

        OUTPUT:

        list

        EXAMPLES::

            sage: from slabbe import Endofunction
            sage: L = [0, 5, 7, 2, 1, 1, 2, 6, 2, 4]
            sage: f = Endofunction(L)
            sage: f.skeleton()
            [0, 5, 7, 1, 2, 6]
        """
        sorted_cycle_elements = sorted(self.cycle_elements())
        L = []
        for a in sorted_cycle_elements:
            L.append(self._list[a])
        return L

    def to_double_rooted_tree(self):
        r"""
        Return the double rooted tree following André Joyal Bijection.

        OUTPUT:

        Double rooted tree

        EXAMPLES::

            sage: from slabbe import Endofunction
            sage: L = [0, 5, 7, 2, 1, 1, 2, 6, 2, 4]
            sage: f = Endofunction(L)
            sage: f.to_double_rooted_tree()
            Double rooted tree:
            Edges: [(0, 5), (1, 2), (2, 6), (3, 2), (4, 1), (5, 7), (7, 1), (8, 2), (9, 4)]
            RootA: 6
            RootB: 0
        """
        tree_list = copy(self._list)
        skeleton = self.skeleton()
        rootA = skeleton[-1]
        rootB = skeleton[0]
        for i in range(len(skeleton)-1):
            tree_list[skeleton[i]] = skeleton[i+1]
        edges = []
        for i,a in enumerate(tree_list):
            if i != rootA:
                edges.append((i,a))
        return DoubleRootedTree(edges, rootA, rootB)

class DoubleRootedTree(object):
    r"""
    Returns a double rooted tree.

    INPUT:

    - ``edges`` - list of edges
    - ``rootA`` - root A
    - ``rootB`` - root B

    EXAMPLES::

        sage: from slabbe import DoubleRootedTree
        sage: edges = [(0,5),(1,2),(2,6),(3,2),(4,1),(5,7),(7,1),(8,2),(9,4)]
        sage: D = DoubleRootedTree(edges, 6, 0)
        sage: D
        Double rooted tree:
        Edges: [(0, 5), (1, 2), (2, 6), (3, 2), (4, 1), (5, 7), (7, 1), (8, 2), (9, 4)]
        RootA: 6
        RootB: 0
    """
    def __init__(self, edges, rootA, rootB):
        r"""
        EXAMPLES::

            sage: from slabbe import DoubleRootedTree
            sage: edges = [(0,5),(1,2),(2,6),(3,2),(4,1),(5,7),(7,1),(8,2),(9,4)]
            sage: D = DoubleRootedTree(edges, 6, 0)
            sage: D
            Double rooted tree:
            Edges: [(0, 5), (1, 2), (2, 6), (3, 2), (4, 1), (5, 7), (7, 1), (8, 2), (9, 4)]
            RootA: 6
            RootB: 0
        """
        self._edges = edges
        self._rootA = rootA
        self._rootB = rootB
        #
        self._graph = None

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from slabbe import DoubleRootedTree
            sage: edges = [(0,5),(1,2),(2,6),(3,2),(4,1),(5,7),(7,1),(8,2),(9,4)]
            sage: D = DoubleRootedTree(edges, 6, 0)
            sage: D
            Double rooted tree:
            Edges: [(0, 5), (1, 2), (2, 6), (3, 2), (4, 1), (5, 7), (7, 1), (8, 2), (9, 4)]
            RootA: 6
            RootB: 0
        """
        s = "Double rooted tree:\n"
        s += "Edges: %s\n" % self._edges
        s += "RootA: %s\n" % self._rootA
        s += "RootB: %s" % self._rootB
        return s

    def __eq__(self, other):
        r"""
        EXAMPLES::

            sage: from slabbe import DoubleRootedTree
            sage: edges = [(0,5),(1,2),(2,6),(3,2),(4,1),(5,7),(7,1),(8,2),(9,4)]
            sage: D = DoubleRootedTree(edges, 6, 0)
            sage: D == 3
            False
            sage: D == D
            True
        """
        return (isinstance(other, DoubleRootedTree) and
                self._rootA == other._rootA and
                self._rootB == other._rootB and
                set(self._edges) == set(other._edges))

    def __ne__(self, other):
        r"""
        EXAMPLES::

            sage: from slabbe import DoubleRootedTree
            sage: edges = [(0,5),(1,2),(2,6),(3,2),(4,1),(5,7),(7,1),(8,2),(9,4)]
            sage: D = DoubleRootedTree(edges, 6, 0)
            sage: D != 3
            True
            sage: D != D
            False
        """
        return not self.__eq__(other)

    def graph(self):
        r"""
        EXAMPLES::

            sage: from slabbe import DoubleRootedTree
            sage: edges = [(0,5),(1,2),(2,6),(3,2),(4,1),(5,7),(7,1),(8,2),(9,4)]
            sage: D = DoubleRootedTree(edges, 6, 0)
            sage: D.graph()
            Graph on 10 vertices
        """
        if self._graph is None:
            G = Graph()
            G.add_vertex(self._rootA)
            G.add_vertex(self._rootB)
            for a,b in self._edges:
                G.add_edge(a,b)
            self._graph = G
        return self._graph

    def skeleton(self):
        r"""
        EXAMPLES::

            sage: from slabbe import DoubleRootedTree
            sage: edges = [(0,5),(1,2),(2,6),(3,2),(4,1),(5,7),(7,1),(8,2),(9,4)]
            sage: D = DoubleRootedTree(edges, 6, 0)
            sage: D.skeleton()
            [0, 5, 7, 1, 2, 6]
        """
        G = self.graph()
        return G.shortest_path(self._rootB, self._rootA)

    def skeleton_cycles(self):
        r"""
        EXAMPLES::

            sage: from slabbe import DoubleRootedTree
            sage: edges = [(0,5),(1,2),(2,6),(3,2),(4,1),(5,7),(7,1),(8,2),(9,4)]
            sage: D = DoubleRootedTree(edges, 6, 0)
            sage: D.skeleton()
            [0, 5, 7, 1, 2, 6]
            sage: D.skeleton_cycles()
            [(0,), (1, 5), (2, 7, 6)]
        """
        skeleton = self.skeleton()
        sorted_skeleton = sorted(skeleton)
        w = Word(skeleton)
        cycles_standard = w.standard_permutation().to_cycles()
        cycles = []
        for cycle_standard in cycles_standard:
            t = tuple(sorted_skeleton[c-1] for c in cycle_standard)
            cycles.append(t)
        return cycles

    def to_endofunction(self):
        r"""

        EXAMPLES::

            sage: from slabbe import DoubleRootedTree
            sage: edges = [(0,5),(1,2),(2,6),(3,2),(4,1),(5,7),(7,1),(8,2),(9,4)]
            sage: D = DoubleRootedTree(edges, 6, 0)
            sage: D.to_endofunction()
            Endofunction:
            [0..9] -> [0, 5, 7, 2, 1, 1, 2, 6, 2, 4]

        TESTS::

            sage: D = DoubleRootedTree([], 0, 0)
            sage: D.to_endofunction()
            Endofunction:
            [0..0] -> [0]
        """
        function = {}
        # define the function for elements in the skeleton
        for cycle in self.skeleton_cycles():
            len_cycle = len(cycle)
            for i in range(len_cycle):
                a = cycle[i]
                b = cycle[(i+1) % len_cycle]
                function[a] = b
        # delete edges of the graph between elements of the skeleton
        skeleton = self.skeleton()
        G = self.graph()
        for i in range(len(skeleton)-1):
            a = skeleton[i]
            b = skeleton[i+1]
            G.delete_edge(a, b)
        # define the function for all other elements
        todo = skeleton
        while todo:
            a = todo.pop()
            for b in G[a]:
                function[b] = a
                G.delete_edge(a, b)
                todo.append(b)
        # reset graph attribute since we destroyed it
        self._graph = None
        # turn the function (dict) into a list
        function_list = [function[i] for i in range(len(function))]
        return Endofunction(function_list)
