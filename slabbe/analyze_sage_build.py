# -*- coding: utf-8 -*-
r"""
A time evolution picture of packages built by Sage

EXAMPLES::

    sage: from slabbe.analyze_sage_build import draw_sage_build
    sage: t = draw_sage_build()
    sage: _ = t.pdf(view=False)

AUTHOR:

    - Sébastien Labbé, December 9-11, 2016
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
from __future__ import absolute_import, print_function
import os, datetime, re
from slabbe.tikz_picture import TikzPicture
from pytimeparse.timeparse import timeparse

def last_modified_datetime(path_to_file):
    r"""
    EXAMPLES::

        sage: from slabbe.analyze_sage_build import last_modified_datetime
        sage: from sage.env import SAGE_ROOT
        sage: last_modified_datetime(SAGE_ROOT+'/README.md')
        datetime.datetime(..., ..., ..., ..., ..., ...)
        sage: last_modified_datetime(SAGE_ROOT+'/VERSION.txt')
        datetime.datetime(..., ..., ..., ..., ..., ...)
    """
    t_stamp = os.path.getmtime(path_to_file)
    return datetime.datetime.fromtimestamp(t_stamp)

def build_duration_logs(path_to_file, pattern=None):
    r"""
    INPUT:

    - ``path_to_file`` -- string, file to parse
    - ``pattern`` -- string or ``None``, string to parse. If ``None``, it
      returns 0 seconds.

    OUTPUT:

        list of timedelta objects

    EXAMPLES::

        sage: L = [ ('/../ptestlong.log', 'Total time for all tests: '),
        ....:       ('/../dochtml.log', 'Elapsed time: '),
        ....:       ('/../start.log', None),
        ....:       ('/sagelib-7.5.beta6.log', 'real\t'),
        ....:       ('/sqlite.log', 'real\t')]

    ::

        sage: from slabbe.analyze_sage_build import build_duration_logs
        sage: from sage.env import SAGE_LOGS
        sage: import os
        sage: for (file, pattern) in L:                  # random
        ....:     print(file)
        ....:     if not os.path.exists(SAGE_LOGS+file):
        ....:         continue
        ....:     build_duration_logs(SAGE_LOGS+file, pattern)
        /../ptestlong.log
        /../dochtml.log
        [datetime.timedelta(0, 471, 700000)]
        /../start.log
        []
        /sagelib-7.5.beta6.log
        /sqlite.log
        []

    """
    if pattern is None:
        return []

    with open(path_to_file, 'r') as f:
        s = f.read()

    result = re.findall('^'+pattern+'.*', s, flags=re.MULTILINE)

    # line = "real\t0m41.027s"            # in logs/*.log
    # line = "Total time for all tests: 14922.7 seconds" # in ptestlong.log
    # line = "Elapsed time: 4156.8 seconds." # in dochtml.log
    L = []
    for line in result:
        time_string = line[len(pattern):].strip('.')
        seconds = timeparse(time_string)
        if seconds is None:
            raise ValueError('Unable to parse time delta in string '
                    '"{}" of file {}'.format(time_string, path_to_file))
        L.append(datetime.timedelta(seconds=seconds))

    return L

def sage_logs_datetime_list(consider='last', verbose=False):
    r"""
    Return a dictionnary of duration and last modified information from the
    sage log files.

    INPUT:

    - ``consider`` - string (default: ``'last'``), ``'last'`` or ``'all'``, in
      any log file, consider the last duration found or the sum of all of them
    - ``verbose`` - boolean (default: ``False``)

    EXAMPLES::

        sage: import datetime
        sage: from slabbe.analyze_sage_build import sage_logs_datetime_list
        sage: L = sage_logs_datetime_list()
        sage: L = sage_logs_datetime_list()
        sage: stop = datetime.datetime.now()
        sage: start = stop - datetime.timedelta(14r)
        sage: L_filtered = [(A,B,delta,file) for (A,B,delta,file) in L
        ....:                  if start <= A and B <= stop]
    """
    from sage.env import SAGE_LOGS
    K = []
    pattern = 'real\t'
    for file in os.listdir(SAGE_LOGS):
        path_to_file = os.path.join(SAGE_LOGS, file)
        K.append((path_to_file, pattern))

    path_to_file = os.path.join(SAGE_LOGS, '..', 'dochtml.log')
    if os.path.exists(path_to_file):
        pattern = 'Elapsed time: '
        K.append((path_to_file, pattern))

    path_to_file = os.path.join(SAGE_LOGS, '..', 'ptestlong.log')
    if os.path.exists(path_to_file):
        pattern = 'Total time for all tests: '
        K.append((path_to_file, pattern))

    path_to_file = os.path.join(SAGE_LOGS, '..', 'start.log')
    if os.path.exists(path_to_file):
        K.append((path_to_file, None))

    L = []
    for path_to_file, pattern in K:
        file = os.path.split(path_to_file)[-1]
        B = last_modified_datetime(path_to_file)
        deltas = build_duration_logs(path_to_file, pattern)
        if pattern is None:
            delta = datetime.timedelta(0)
        elif not deltas:
            if verbose:
                print('Warning: no duration found in file {} '
                      'with pattern "{}"'.format(file, pattern))
            delta = datetime.timedelta(0)
        elif consider == 'last':
            delta = deltas[-1]
        elif consider == 'all':
            delta = sum(deltas, datetime.timedelta(0))
        else:
            raise ValueError('Unknown value for consider(={})'.format(consider))
        entry = (B-delta,B,delta,file)
        L.append(entry)

    return L

def draw_sage_build(start=None, stop=None, consider='last', verbose=False):
    r"""
    Return a time evolution picture of packages built by Sage

    INPUT:

    - ``start`` -- datetime object (default is January 1st, 2000)
    - ``stop`` -- datetime object (default is now)
    - ``consider`` - string (default: ``'last'``), ``'last'`` or ``'all'``, in
      any log file, consider the last duration found or the sum of all of them
    - ``verbose`` - boolean (default: ``False``)

    OUTPUT:

    TikzPicture object

    EXAMPLES::

        sage: from slabbe.analyze_sage_build import draw_sage_build
        sage: t = draw_sage_build()

    During the previous 7 day::

        sage: import datetime
        sage: stop = datetime.datetime.now()
        sage: start = stop - datetime.timedelta(7r)
        sage: t = draw_sage_build(start, stop)        # not tested

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
    L_all = sage_logs_datetime_list(consider=consider, verbose=verbose)
    if not L_all:
        raise ValueError("no package found in the logs directory")
    L_all.sort()

    # Filter data according to desired period
    L = [(A,B,delta,file) for (A,B,delta,file) in L_all
                          if start <= A and B <= stop]
    if not L:
        raise ValueError("no package found built between start date (={}) "
                         "and stop date (={}). Oldest started at {}. "
                         "Newest started at {}.".format(start, stop,
                                        L_all[0][0], L_all[-1][0]))
    # time to x axis float
    effective_start = min(A for (A,B,delta,file) in L)
    effective_stop = max(B for (A,B,delta,file) in L)
    duration = (effective_stop-effective_start)
    TOTAL = duration.total_seconds()
    if TOTAL == 0:
        TOTAL += 1 # avoid division by zero

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
    L.sort()
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

