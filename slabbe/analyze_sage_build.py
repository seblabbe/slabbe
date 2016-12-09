# -*- coding: utf-8 -*-
r"""
A time evolution picture of Sage building its packages

EXAMPLES::

    sage: from slabbe.analyze_sage_build import draw_sage_build
    sage: t = draw_sage_build() # not tested
    sage: t.pdf()               # not tested

AUTHOR:

    - Sébastien Labbé, December 9, 2016
"""
#*****************************************************************************
#       Copyright (C) 2016 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import
import os, datetime, re
from .tikz_picture import TikzPicture

def last_modified_datetime(path_to_file):
    r"""
    EXAMPLES::

        sage: from slabbe.analyze_sage_build import last_modified_datetime
        sage: from sage.env import SAGE_ROOT
        sage: last_modified_datetime(SAGE_ROOT+'/README.md')
        datetime.datetime(2016, 8, 25, 17, 13, 49)
        sage: last_modified_datetime(SAGE_ROOT+'/VERSION.txt')
        datetime.datetime(2016, 12, 2, 0, 46, 46)
    """
    t_stamp = os.path.getmtime(path_to_file)
    return datetime.datetime.fromtimestamp(t_stamp)

def build_duration(path_to_file):
    r"""
    EXAMPLES::

        sage: from slabbe.analyze_sage_build import build_duration
        sage: from sage.env import SAGE_LOGS
        sage: build_duration(SAGE_LOGS+'/sagelib-7.5.beta5.log')
        datetime.timedelta(0, 9, 747000)
        sage: build_duration(SAGE_LOGS+'/sqlite.log')
        Traceback (most recent call last):
        ...
        ValueError: No duration found in file /Users/slabbe/Applications/sage-git/logs/pkgs/sqlite.log
    """
    with open(path_to_file, 'r') as f:
        s = f.read()

    result = re.findall('real\t.*', s)
    if not result:
        raise ValueError("No duration found in file {}".format(path_to_file))

    line = result[-1]
    pos_m = line.find('m')
    pos_s = line.find('s')
    minutes = int(line[5:pos_m])
    seconds = float(line[pos_m+1:pos_s])
    return datetime.timedelta(minutes=minutes, seconds=seconds)

def sage_logs_datetime_list(verbose=False):
    r"""
    Return a dictionnary of duration and last modified information from the
    sage log files.

    EXAMPLES::

        sage: import datetime
        sage: from slabbe.analyze_sage_build import sage_logs_datetime_list
        sage: L = sage_logs_datetime_list()
        sage: L = sage_logs_datetime_list()
    """
    L = []
    from sage.env import SAGE_LOGS
    for file in os.listdir(SAGE_LOGS):
        path_to_file = os.path.join(SAGE_LOGS, file)
        B = last_modified_datetime(path_to_file)
        try:
            delta  = build_duration(path_to_file)
        except ValueError:
            if verbose:
                print("Warning: no duration found in file {} ...".format(file))
            delta = datetime.timedelta(0)
        entry = (B-delta,B,delta,file)
        L.append(entry)
    return L

def draw_sage_build(start=None, stop=None):
    r"""
    INPUT:

    - ``start`` -- datetime object (default is January 1st, 2000)
    - ``stop`` -- datetime object (default is now)

    OUTPUT:

    TikzPicture object

    EXAMPLES::

        sage: from slabbe.analyze_sage_build import draw_sage_build
        sage: t = draw_sage_build()

    During the previous 7 day::

        sage: import datetime
        sage: stop = datetime.datetime.now()
        sage: start = stop - datetime.timedelta(7r)
        sage: t = draw_sage_build(start, stop)       # not tested

    ::

        sage: t = draw_sage_build(stop=datetime.datetime(2016,12,9)) # not tested

    TESTS:

    It may fail if nothing is found during the period::

        sage: t = draw_sage_build()  # not tested
        Traceback (most recent call last):
        ...
        ValueError: no package found built between start date (=2016-12-08
        16:58:35.189127) and stop date (=2016-12-09 16:58:35.189127). Oldest is
        2016-02-16 11:38:00.673327. Newest is 2016-12-08 11:27:03.312340.
    """
    if stop is None:
        stop = datetime.datetime.now()
    if start is None:
        #start = stop - datetime.timedelta(1) # minus 1 day
        start = datetime.datetime(2000, 1, 1)
    assert start < stop, ("start date (={}) must precede "
                          "stop date (={})".format(start, stop))

    # Obtain datetime and timedelta from logs
    L_all = sage_logs_datetime_list()
    if not L_all:
        raise ValueError("no package found in the logs directory")
    L_all.sort()

    # Filter data according to desired period
    L = [(A,B,delta,file) for (A,B,delta,file) in L_all
                          if start <= A and B <= stop]
    if not L:
        raise ValueError("no package found built between start date (={}) "
                         "and stop date (={}). Oldest is {}. "
                         "Newest is {}.".format(start, stop, 
                                        L_all[0][0], L_all[-1][0]))
    # time to x axis float
    effective_start = L[0][0]
    effective_stop = L[-1][1]
    duration = (effective_stop-effective_start)
    TOTAL = duration.total_seconds()
    def timedelta_to_float(t):
        return t.total_seconds() / TOTAL * 40
    def datetime_to_float(t):
        return timedelta_to_float(t-effective_start)

    # Tikz code
    lines = []
    lines.append('\\begin{tikzpicture}[scale=0.5]')
    str_rect = '\\draw[fill=orange] ({}, {}) rectangle ({}, {}) ;'
    str_node = '\\node[{}] at ({}, {}) {{{}}};'

    # Titles
    lines.append(str_node.format('right',0,5,
                                 effective_start.isoformat()+' (start)'))
    lines.append(str_node.format('right',0,4,
                                 effective_stop.isoformat()+' (end)'))
    lines.append(str_node.format('right',0,3,
                                 '{} (duration)'.format(duration)))

    # Axis
    t0 = datetime_to_float(effective_start)
    tn = datetime_to_float(effective_stop)
    axis_y = 1.5
    lines.append('\\draw[thick,->] ({},{axis_y}) -- '
                                  '({},{axis_y});'.format(t0, tn, axis_y=axis_y))
    duration_minutes = int(duration.total_seconds()/60)
    step = duration_minutes // 4
    step = 10 * (step // 10 + 1)
    for n_minutes in range(0, duration_minutes, step):
        delta = datetime.timedelta(0,n_minutes*60)
        tj = timedelta_to_float(delta)
        lines.append('\\draw[thick] ({tj},{}) -- ({tj},{}) node[below] '
                '{{{delta}}};'.format(axis_y+.3, axis_y-.3, delta=delta, tj=tj))

    # Boxes
    for i,entry in enumerate(L):
        A,B,delta,file_ext = entry
        Af = datetime_to_float(A)
        Bf = datetime_to_float(B)
        lines.append(str_rect.format(Af,-i,Bf,-i-1))
        file,ext = os.path.splitext(file_ext)
        file = file.replace('_','\\_')
        file = '\\texttt{{{}}}'.format(file)
        lines.append(str_node.format('right',Bf,-i-.5,file))

    lines.append('\\end{tikzpicture}')
    return TikzPicture('\n'.join(lines))

