# -*- coding: utf-8 -*-
r"""
Some fruits

This file is an example of Sage conventions for syntax and documentation.
It also serves as an example to explain class inheritance: methods defined
for fruits are automatically defined for bananas and strawberries.

AUTHORS:

- Sébastien Labbé, 2010-2013

EXAMPLES::

    sage: from slabbe import Fruit
    sage: f = Fruit(5)
    sage: f
    A fruit of 5 kilos.
    sage: f.weight()
    5
    sage: f.is_a_fruit()
    True

Because of class inheritance which says that a banana is a fruit and a
strawberry is a fruit, the methods written for fruits are defined for
bananas and strawberries automatically::

    sage: from slabbe import Strawberry
    sage: s = Strawberry(32)
    sage: s
    A strawberry of 32 kilos.
    sage: s.weight()
    32
    sage: s.is_a_fruit()
    True

::

    sage: from slabbe import Banana
    sage: b = Banana(13)
    sage: b
    A banana of 13 kilos.
    sage: b.weight()
    13
    sage: b.is_a_fruit()
    True

::

    sage: s = Strawberry(32)
    sage: t = Strawberry(7)
    sage: s
    A strawberry of 32 kilos.
    sage: t
    A strawberry of 7 kilos.
    sage: s + t
    A strawberry of 1073 kilos.

"""
#*****************************************************************************
#       Copyright (C) 2010-2013 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
from sage.structure.sage_object import SageObject

class Fruit(SageObject):
    r"""
    Creates a fruit.

    INPUT:

    - ``weight`` - number, in kilos

    OUTPUT:

        Fruit

    EXAMPLES::

        sage: from slabbe import Fruit
        sage: f = Fruit(5)
        sage: f
        A fruit of 5 kilos.

    """
    def __init__(self, weight=1):
        r"""
        Constructor.

        EXAMPLES::

            sage: from slabbe import Fruit
            sage: f = Fruit(3)

        TESTS:

        Testing that pickle works::

            sage: loads(dumps(f))
            A fruit of 3 kilos.
            sage: f == loads(dumps(f))
            True

        Running the test suite::

            sage: TestSuite(f).run()
        """
        self._weight = weight

    def _repr_(self):
        r"""
        Returns the string representation.

        OUTPUT:

            String

        EXAMPLES::

            sage: from slabbe import Fruit
            sage: f = Fruit(3)
            sage: f
            A fruit of 3 kilos.
        """
        return "A fruit of %s kilos."%self._weight

    def __eq__(self, other):
        r"""
        Return if self is equal to other.

        INPUT:

        - ``other`` - an object

        OUTPUT:

            Boolean

        EXAMPLES::

            sage: from slabbe import Fruit
            sage: s = Fruit(32)
            sage: t = Fruit(7)
            sage: u = Fruit(32)
            sage: s == u
            True
            sage: s == t
            False
        """
        return (isinstance(other, Fruit) and
                self.weight() == other.weight())
    def weight(self):
        r"""
        Return the weight.

        OUTPUT:

            Number

        EXAMPLES::

            sage: from slabbe import Fruit
            sage: f = Fruit(3)
            sage: f.weight()
            3
        """
        return self._weight

    def is_a_fruit(self):
        r"""
        Returns True if it is a fruit.

        OUTPUT:

            Boolean

        EXAMPLES::

            sage: from slabbe import Fruit
            sage: f = Fruit(3)
            sage: f.is_a_fruit()
            True
        """
        return True

class Banana(Fruit):
    r"""
    Creates a banana.

    INPUT:

    - ``weight`` - number, in kilos

    OUTPUT:

        Banana

    EXAMPLES::

        sage: from slabbe import Banana
        sage: b = Banana(9)
        sage: b
        A banana of 9 kilos.

    TESTS:

    Testing that pickle works::

        sage: loads(dumps(b))
        A banana of 9 kilos.
        sage: b == loads(dumps(b))
        True

    Running the test suite::

        sage: TestSuite(b).run()
    """
    def _repr_(self):
        r"""
        Returns the string representation.

        OUTPUT:

            String

        EXAMPLES::

            sage: from slabbe import Banana
            sage: b = Banana(9)
            sage: b
            A banana of 9 kilos.
        """
        return "A banana of %s kilos."%self._weight

    def __eq__(self, other):
        r"""
        Return if self is equal to other.

        INPUT:

        - ``other`` - an object

        OUTPUT:

            Boolean

        EXAMPLES::

            sage: from slabbe import Banana
            sage: s = Banana(32)
            sage: t = Banana(7)
            sage: u = Banana(32)
            sage: s == u
            True
            sage: s == t
            False
        """
        return (isinstance(other, Banana) and
                self.weight() == other.weight())
class Strawberry(Fruit):
    r"""
    Creates a strawberry.

    INPUT:

    - ``weight`` - number, in kilos

    OUTPUT:

        Strawberry

    EXAMPLES::

        sage: from slabbe import Strawberry
        sage: s = Strawberry(34)
        sage: s
        A strawberry of 34 kilos.

    TESTS:

    Testing that pickle works::

        sage: loads(dumps(s))
        A strawberry of 34 kilos.
        sage: s == loads(dumps(s))
        True

    Running the test suite::

        sage: TestSuite(s).run()
    """
    def _repr_(self):
        r"""
        Returns the string representation.

        OUTPUT:

            String

        EXAMPLES::

            sage: from slabbe import Strawberry
            sage: f = Strawberry(34)
            sage: f
            A strawberry of 34 kilos.
        """
        return "A strawberry of %s kilos."%self._weight

    def __eq__(self, other):
        r"""
        Return if self is equal to other.

        INPUT:

        - ``other`` - an object

        OUTPUT:

            Boolean

        EXAMPLES::

            sage: from slabbe import Strawberry
            sage: s = Strawberry(32)
            sage: t = Strawberry(7)
            sage: u = Strawberry(32)
            sage: s == u
            True
            sage: s == t
            False
        """
        return (isinstance(other, Strawberry) and
                self.weight() == other.weight())
    def __add__(self, other):
        r"""
        Fusion of a strawberry with another strawberry.

        INPUT:

        - ``other`` - a strawberry

        OUTPUT:

            Strawberry

        EXAMPLES::

            sage: from slabbe import Strawberry
            sage: f = Strawberry(34)
            sage: g = Strawberry(20)
            sage: f + g
            A strawberry of 1556 kilos.
        """
        return Strawberry(self.weight()**2 + other.weight()**2)

