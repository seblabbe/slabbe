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

def print_set_to_file(s, filename, ncolumns=4, start=r'\{', end=r'\}',
        sep=r',', eol=r',\\&\hspace{10mm}'):
    r"""
    INPUT:

    - ``s`` -- string
    - ``filename`` -- string, filename

    EXAMPLES::

        sage: from slabbe.write_to_file import print_set_to_file
        sage: filename = tmp_filename(ext='.tex')
        sage: M = [random_matrix(ZZ, 3, 3) for _ in range(20)]
        sage: print_set_to_file(M, filename, ncolumns=8, 
        ....:          start=r'\left\{\begin{array}{cccccccc}',
        ....:          end=r'\end{array}\right\}',
        ....:          sep=r',&', eol=r',\\[4mm]')
        Creation of ...tex

    """
    from sage.misc.latex import latex
    L = [start]
    for i,a in enumerate(s):
        L.append(latex(a))
        if (i+1) == len(s):
            L.append(end)
        elif (i+1) % ncolumns == 0:
            L.append(eol)
        else:
            L.append(sep)
    output = '\n'.join(L)
    write_str_to_file(output, filename)

