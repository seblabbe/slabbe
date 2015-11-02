# -*- coding: utf-8 -*-
r"""
Should be in Sage methods

AUTHORS:

 - Sebastien Labbe, 2015

"""
#*****************************************************************************
#       Copyright (C) 2015 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

####################
# WordMorphism class
####################
def image_first_letter_to_letter(self):
    r"""
    EXAMPLES::

        sage: from slabbe.should_be_in_sage import image_first_letter_to_letter
        sage: m = WordMorphism({1: [1, 2], 2: [2], 3: [3]})
        sage: image_first_letter_to_letter(m)
        defaultdict(<type 'list'>, {1: [1], 2: [2], 3: [3]})

    ::

        sage: m = WordMorphism({1: [1], 2: [1, 3], 3: [2]})
        sage: image_first_letter_to_letter(m)
        defaultdict(<type 'list'>, {1: [1, 2], 2: [3]})

    ::

        sage: m = WordMorphism({1: [2], 2: [1, 3], 3: [3]})
        sage: image_first_letter_to_letter(m)
        defaultdict(<type 'list'>, {1: [2], 2: [1], 3: [3]})
    """
    from collections import defaultdict
    D = defaultdict(list)
    for a in self.domain().alphabet():
        image = self.image(a)
        D[image[0]].append(a)
    return D
