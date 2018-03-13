# -*- coding: utf-8 -*-
r"""
Writing to files

EXAMPLES::

    sage: from slabbe.write_to_file import write_str_to_file
    sage: filename = tmp_filename(ext='.tex')
    sage: write_str_to_file('aaa', filename)
    Creation of ...tex

AUTHOR:

    - Sébastien Labbé, March 13, 2018
"""
#*****************************************************************************
#       Copyright (C) 2018 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function

def write_str_to_file(s, filename):
    r"""
    INPUT:

    - ``s`` -- string
    - ``filename`` -- string, filename

    EXAMPLES::

        sage: from slabbe.write_to_file import write_str_to_file
        sage: filename = tmp_filename(ext='.tex')
        sage: write_str_to_file('aaa', filename)
        Creation of ...tex
    """
    with open(filename, 'w') as f:
        print("Creation of {}".format(filename))
        f.write(s)
