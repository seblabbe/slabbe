# -*- coding: utf-8 -*-
r"""
Sage build analyze

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
import os, datetime, re

def last_modified_datetime(path_to_file):
    r"""
    EXAMPLES::

        sage: from sage.env import SAGE_ROOT
        sage: last_modified_datetime(SAGE_ROOT+'/README.md')
        datetime.datetime(2016, 8, 25, 17, 13, 49)
        sage: last_modified_datetime(SAGE_ROOT+'/VERSION.txt')
        datetime.datetime(2016, 12, 2, 0, 46, 46)
    """
    t_stamp = os.path.getmtime(path_to_file)
    return datetime.datetime.fromtimestamp(t_stamp)

def duration_last_build(path_to_file):
    r"""
    EXAMPLES::

        sage: from sage.env import SAGE_LOGS
        sage: duration_last_build(SAGE_LOGS+'/sagelib-7.5.beta5.log')
        datetime.timedelta(0, 9, 747000)
        sage: duration_last_build(SAGE_LOGS+'/sqlite.log')
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

def two_dict():
    r"""
    EXAMPLES::

        sage: modified, duration = two_dict()
        Error: no duration found in file sqlite.log
        sage: modified
        {'sagelib-7.5.beta5.log': datetime.datetime(2016, 12, 2, 2, 10, 14),
         'sqlite.log': datetime.datetime(2016, 12, 2, 4, 36, 51)}
        sage: duration
        {'sagelib-7.5.beta5.log': datetime.timedelta(0, 9, 747000)}
    """
    modified = {}
    duration = {}
    from sage.env import SAGE_LOGS
    for file in os.listdir(SAGE_LOGS):
        path_to_file = os.path.join(SAGE_LOGS, file)
        modified[file] = last_modified_datetime(path_to_file)
        try:
            duration[file] = duration_last_build(path_to_file)
        except ValueError:
            print "Error: no duration found in file {}".format(file)
            duration[file] = 0
    return modified, duration

    






