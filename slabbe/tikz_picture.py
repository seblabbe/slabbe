# -*- coding: utf-8 -*-
r"""
TikzPicture

A Python Module for tikz pictures. A TikzPicture object is created from a string
starting with ``r'\begin{tikzpicture}'`` and ending with
``r'\end{tikzpicture}'``.

The module allows easy creation of tikz pictures from Sage objects like graphs
and posets. Conversion of tikz pictures to pdf and png format based on
standalone LaTeX document class.

EXAMPLES::

    sage: from slabbe import TikzPicture
    sage: V = [[1,0,1],[1,0,0],[1,1,0],[0,0,-1],[0,1,0],[-1,0,0],[0,1,1],[0,0,1],[0,-1,0]]
    sage: P = Polyhedron(vertices=V).polar()
    sage: s = P.projection().tikz([674,108,-731],112)
    sage: t = TikzPicture(s)

Creation of a pdf in a temporary directory. The returned value is a string
giving the file path::

    sage: path_to_file = t.pdf(view=False)   # long time (2s)

Setting ``view=True``, which is the default, opens the pdf in a viewer.

::

    sage: t
    \documentclass[tikz]{standalone}
    \usepackage{amsmath}
    \begin{document}
    \begin{tikzpicture}%
            [x={(0.249656cm, -0.577639cm)},
            y={(0.777700cm, -0.358578cm)},
            z={(-0.576936cm, -0.733318cm)},
            scale=2.000000,
    ...
    ... 80 lines not printed (4889 characters in total) ...
    ...
    \node[vertex] at (1.00000, 1.00000, -1.00000)     {};
    \node[vertex] at (1.00000, 1.00000, 1.00000)     {};
    %%
    %%
    \end{tikzpicture}
    \end{document}

Use ``print t`` to see the complete content of the file.

Adding a border avoids croping the vertices of a graph::

    sage: g = graphs.PetersenGraph()
    sage: s = latex(g)   # takes 3s but the result is cached
    sage: t = TikzPicture(s, standalone_options=["border=4mm"], usepackage=['tkz-graph'])
    sage: _ = t.pdf()    # not tested

If dot2tex Sage optional package and graphviz are installed, then the following
one liner works::

    sage: t = TikzPicture.from_graph(g)  # optional: dot2tex # long time (3s)

::

    sage: s = latex(transducers.GrayCode())
    sage: t = TikzPicture(s, usetikzlibrary=['automata'])
    sage: _ = t.pdf(view=False)  # long time (2s)

AUTHORS:

- Sébastien Labbé, initial version in slabbe-0.2.spkg, nov 2015.
"""
#*****************************************************************************
#       Copyright (C) 2015-2017 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import, print_function
from sage.misc.latex import have_pdflatex, have_convert, have_program
from sage.misc.temporary_file import tmp_filename
from sage.structure.sage_object import SageObject
import os

class TikzPicture(SageObject):
    r"""
    Creates a TikzPicture embedded in a LaTeX standalone document class.

    INPUT:

    - ``code`` -- string, tikzpicture code starting with ``r'\begin{tikzpicture}'``
      and ending with ``r'\end{tikzpicture}'``
    - ``standalone_options`` -- list of strings (default: ``[]``),
      latex document class standalone configuration options.
    - ``usepackage`` -- list of strings (default: ``['amsmath']``), latex
      packages.
    - ``usetikzlibrary`` -- list of strings (default: ``[]``), tikz libraries
      to use.
    - ``macros`` -- list of strings (default: ``[]``), stuff you need for the picture.
    - ``use_sage_preamble`` -- bool (default: ``False``), whether to include sage
      latex preamble and sage latex macros, that is, the content of
      :func:`sage.misc.latex.extra_preamble()`,
      :func:`sage.misc.latex.extra_macros()` and
      :func:`sage.misc.latex_macros.sage_latex_macros()`.

    EXAMPLES::

        sage: from slabbe import TikzPicture
        sage: g = graphs.PetersenGraph()
        sage: s = latex(g)
        sage: t = TikzPicture(s, standalone_options=["border=4mm"], usepackage=['tkz-graph'])
        sage: _ = t.pdf(view=False)   # long time (2s)

    Here are standalone configurations, packages, tikz libraries and macros you
    may want to set::

        sage: options = ['preview', 'border=4mm', 'beamer', 'float']
        sage: usepackage = ['nicefrac', 'amsmath', 'pifont', 'tikz-3dplot',
        ....:    'tkz-graph', 'tkz-berge', 'pgfplots']
        sage: tikzlib = ['arrows', 'snakes', 'backgrounds', 'patterns',
        ....:      'matrix', 'shapes', 'fit', 'calc', 'shadows', 'plotmarks',
        ....:      'positioning', 'pgfplots.groupplots', 'mindmap']
        sage: macros = [r'\newcommand{\ZZ}{\mathbb{Z}}']
        sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
        sage: t = TikzPicture(s, standalone_options=options, usepackage=usepackage, 
        ....:        usetikzlibrary=tikzlib, macros=macros)
        sage: _ = t.pdf(view=False)   # long time (2s)
    """
    def __init__(self, code, standalone_options=None, usepackage=['amsmath'],
            usetikzlibrary=None, macros=None, use_sage_preamble=False):
        r"""
        See the class documentation for full information.

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
        """
        self._code = code
        self._standalone_options = [] if standalone_options is None else standalone_options
        self._usepackage = usepackage
        self._usetikzlibrary = [] if usetikzlibrary is None else usetikzlibrary
        self._macros = [] if macros is None else macros
        if use_sage_preamble:
            from sage.misc.latex import _Latex_prefs
            for key in ['preamble', 'macros']:
                s = _Latex_prefs._option[key]
                if s: 
                    self._macros.append(s)
            from sage.misc.latex_macros import sage_latex_macros
            self._macros.extend(sage_latex_macros())

    @classmethod
    def from_graph(cls, graph, merge_multiedges=True,
            merge_label_function=tuple, **kwds):
        r"""
        Convert a graph to a tikzpicture using graphviz and dot2tex.

        .. NOTE::

            Prerequisite: dot2tex optional Sage package and graphviz must be
            installed.

        INPUT:

        - ``graph`` -- graph
        - ``merge_multiedges`` -- bool (default: ``True``), if the graph
          has multiple edges, whether to merge the multiedges into one
          single edge
        - ``merge_label_function`` -- function (default:``tuple``), a
          function to apply to each list of labels to be merged. It is
          ignored if ``merge_multiedges`` is not ``True`` or if the graph
          has no multiple edges.

        Other inputs are used for latex drawing with dot2tex and graphviz:

        - ``prog`` -- string (default: ``'dot'``) the program used for the
          layout corresponding to one of the software of the graphviz
          suite: 'dot', 'neato', 'twopi', 'circo' or 'fdp'.
        - ``edge_labels`` -- bool (default: ``True``)
        - ``color_by_label`` -- bool (default: ``False``)
        - ``rankdir`` -- string (default: ``'down'``)

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: g = graphs.PetersenGraph()
            sage: tikz = TikzPicture.from_graph(g) # optional dot2tex # long time (3s)
            sage: _ = tikz.pdf()      # not tested

        Using ``prog``::

            sage: tikz = TikzPicture.from_graph(g, prog='neato', color_by_label=True) # optional dot2tex # long time (3s)
            sage: _ = tikz.pdf()      # not tested

        Using ``rankdir``::

            sage: tikz = TikzPicture.from_graph(g, rankdir='right') # optional dot2tex # long time (3s)
            sage: _ = tikz.pdf()      # not tested

        Using ``merge_multiedges``::

            sage: alpha = var('alpha')
            sage: m = matrix(2,range(4)); m.set_immutable()
            sage: G = DiGraph([(0,1,alpha), (0,1,'a'), (0,2,9), (0,2,m)], multiedges=True)
            sage: tikz = TikzPicture.from_graph(G, merge_multiedges=True) # optional dot2tex
            sage: _ = tikz.pdf()      # not tested

        Using ``merge_multiedges`` with ``merge_label_function``::

            sage: fn = lambda L: LatexExpr(','.join(map(str, L)))
            sage: G = DiGraph([(0,1,'a'), (0,1,'b'), (0,2,'c'), (0,2,'d')], multiedges=True)
            sage: tikz = TikzPicture.from_graph(G, merge_multiedges=True,
            ....:               merge_label_function=fn) # optional dot2tex
            sage: _ = tikz.pdf()      # not tested

        """
        if merge_multiedges and graph.has_multiple_edges():
            from slabbe.graph import merge_multiedges
            graph = merge_multiedges(graph,
                    label_function=merge_label_function)

        default = dict(format='dot2tex', edge_labels=True,
                       color_by_label=False, prog='dot', rankdir='down')
        default.update(kwds)

        graph.latex_options().set_options(**default)
        tikz = graph._latex_()
        return TikzPicture(tikz, standalone_options=["border=4mm"])

    @classmethod
    def from_graph_with_pos(cls, graph, scale=1, merge_multiedges=True,
            merge_label_function=tuple):
        r"""
        Convert a graph with positions defined for vertices to a tikzpicture.

        INPUT:

        - ``graph`` -- graph (with predefined positions)
        - ``scale`` -- number (default:``1``), tikzpicture scale
        - ``merge_multiedges`` -- bool (default: ``True``), if the graph
          has multiple edges, whether to merge the multiedges into one
          single edge
        - ``merge_label_function`` -- function (default:``tuple``), a
          function to apply to each list of labels to be merged. It is
          ignored if ``merge_multiedges`` is not ``True`` or if the graph
          has no multiple edges.

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: g = graphs.PetersenGraph()
            sage: tikz = TikzPicture.from_graph_with_pos(g)

        ::

            sage: edges = [(0,0,'a'),(0,1,'b'),(0,1,'c')]
            sage: kwds = dict(format='list_of_edges', loops=True, multiedges=True)
            sage: G = DiGraph(edges, **kwds)
            sage: G.set_pos({0:(0,0), 1:(1,0)})
            sage: f = lambda label:','.join(label)
            sage: TikzPicture.from_graph_with_pos(G, merge_label_function=f)
            \documentclass[tikz]{standalone}
            \standaloneconfig{border=4mm}
            \usepackage{amsmath}
            \begin{document}
            \begin{tikzpicture}
            [auto,scale=1]
            % vertices
            \node (node_0) at (0, 0) {0};
            \node (node_1) at (1, 0) {1};
            % edges
            \draw[->] (node_0) -- node {b,c} (node_1);
            % loops
            \draw (node_0) edge [loop above] node {a} ();
            \end{tikzpicture}
            \end{document}

        TESTS::

            sage: edges = [(0,0,'a'),(0,1,'b'),(0,1,'c')]
            sage: kwds = dict(format='list_of_edges', loops=True, multiedges=True)
            sage: G = DiGraph(edges, **kwds)
            sage: TikzPicture.from_graph_with_pos(G)
            Traceback (most recent call last):
            ...
            ValueError: vertex positions need to be set first
        """
        pos = graph.get_pos()
        if pos is None:
            raise ValueError('vertex positions need to be set first')

        if merge_multiedges and graph.has_multiple_edges():
            from slabbe.graph import merge_multiedges
            graph = merge_multiedges(graph,
                    label_function=merge_label_function)

        keys_for_vertices = graph._keys_for_vertices()

        lines = []
        lines.append(r'\begin{tikzpicture}')
        lines.append(r'[auto,scale={}]'.format(scale))

        # vertices
        lines.append(r'% vertices')
        for u in graph.vertices():
            line = r'\node ({}) at {} {{{}}};'.format(keys_for_vertices(u),
                                                      pos[u], u)
            lines.append(line)

        # edges
        lines.append(r'% edges')
        arrow = '->' if graph.is_directed() else ''
        for (u,v,label) in graph.edges():
            if u == v:
                # loops are done below
                continue
            if label:
                line = r'\draw[{}] ({}) -- node {{{}}} ({});'.format(arrow,
                                                    keys_for_vertices(u),
                                                    label,
                                                    keys_for_vertices(v))
            else:
                line = r'\draw[{}] ({}) -- ({});'.format(arrow,
                                                    keys_for_vertices(u),
                                                    keys_for_vertices(v))
            lines.append(line)

        # loops
        lines.append(r'% loops')
        for (u,v,label) in graph.loop_edges():
            line = r'\draw ({}) edge [loop above] node {{{}}} ();'.format(
                                              keys_for_vertices(u), label)
            lines.append(line)

        lines.append(r'\end{tikzpicture}')
        tikz = '\n'.join(lines)
        return TikzPicture(tikz, standalone_options=["border=4mm"])

    @classmethod
    def from_poset(cls, poset, **kwds):
        r"""
        Convert a poset to a tikzpicture using graphviz and dot2tex.

        .. NOTE::

            Prerequisite: dot2tex optional Sage package and graphviz must be
            installed.

        INPUT:

        - ``poset`` -- poset
        - ``prog`` -- string (default: ``'dot'``) the program used for the
          layout corresponding to one of the software of the graphviz
          suite: 'dot', 'neato', 'twopi', 'circo' or 'fdp'.
        - ``edge_labels`` -- bool (default: ``True``)
        - ``color_by_label`` -- bool (default: ``False``)
        - ``rankdir`` -- string (default: ``'down'``)

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: P = posets.PentagonPoset()
            sage: tikz = TikzPicture.from_poset(P) # optional dot2tex # long time (3s)
            sage: tikz = TikzPicture.from_poset(P, prog='neato', color_by_label=True) # optional dot2tex # long time (3s)

        ::

            sage: P = posets.SymmetricGroupWeakOrderPoset(4)
            sage: tikz = TikzPicture.from_poset(P) # optional dot2tex # long time (4s)
            sage: tikz = TikzPicture.from_poset(P, prog='neato') # optional dot2tex # long time (4s)
        """
        graph = poset.hasse_diagram()
        return cls.from_graph(graph, **kwds)

    def _latex_file_header_lines(self):
        r"""
        EXAMPLES::

            sage: latex.extra_preamble('')
            sage: from slabbe import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s, standalone_options=["border=4mm"], usepackage=['tkz-graph'])
            sage: t._latex_file_header_lines()[:6]
            ['\\documentclass[tikz]{standalone}',
             '\\standaloneconfig{border=4mm}',
             '\\usepackage{tkz-graph}']
        """
        lines = []
        lines.append(r"\documentclass[tikz]{standalone}")
        for config in self._standalone_options:
            lines.append(r"\standaloneconfig{{{}}}".format(config))
        for package in self._usepackage:
            lines.append(r"\usepackage{{{}}}".format(package))
        for library in self._usetikzlibrary:
            lines.append(r"\usetikzlibrary{{{}}}".format(library))
        lines.extend(self._macros)
        return lines

    def _repr_(self):
        r"""
        Returns the first few and last few lines.

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: g = graphs.PetersenGraph()
            sage: s = latex(g)
            sage: t = TikzPicture(s, usepackage=['tkz-graph'])
            sage: t
            \documentclass[tikz]{standalone}
            \usepackage{tkz-graph}
            \begin{document}
            \begin{tikzpicture}
            \definecolor{cv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cfv0}{rgb}{1.0,1.0,1.0}
            \definecolor{clv0}{rgb}{0.0,0.0,0.0}
            \definecolor{cv1}{rgb}{0.0,0.0,0.0}
            ...
            ... 65 lines not printed (3695 characters in total) ...
            ...
            \Edge[lw=0.1cm,style={color=cv6v8,},](v6)(v8)
            \Edge[lw=0.1cm,style={color=cv6v9,},](v6)(v9)
            \Edge[lw=0.1cm,style={color=cv7v9,},](v7)(v9)
            %
            \end{tikzpicture}
            \end{document}
        """
        lines = self._latex_file_header_lines()
        lines.append(r"\begin{document}")
        L = self._code.splitlines()
        if len(L) <= 10:
            lines.extend(L)
        else:
            lines.extend(L[:5])
            lines.append('...')
            lines.append('... {} lines not printed ({} characters in total) ...'.format(len(L)-10, 
                                                           len(self._code)))
            lines.append('...')
            lines.extend(L[-5:])
        lines.append(r"\end{document}")
        return '\n'.join(lines)

    def __str__(self):
        r"""
        Returns the complete string.

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
            sage: print(t)
            \documentclass[tikz]{standalone}
            \usepackage{amsmath}
            \begin{document}
            \begin{tikzpicture}
            \draw (0,0) -- (1,1);
            \end{tikzpicture}
            \end{document}
        """
        lines = self._latex_file_header_lines()
        lines.append(r"\begin{document}")
        lines.append(self._code)
        lines.append(r"\end{document}")
        return '\n'.join(lines)

    def tikz_picture_code(self):
        r"""
        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
            sage: print(t.tikz_picture_code())
            \begin{tikzpicture}
            \draw (0,0) -- (1,1);
            \end{tikzpicture}
        """
        return self._code

    def pdf(self, filename=None, view=True):
        """
        Compiles the latex code with pdflatex and create a pdf file.

        INPUT:

        - ``filename`` -- string (default:``None``), the output filename. 
          If ``None``, it saves the file in a temporary directory.

        - ``view`` -- bool (default:``True``), whether to open the file in a
          pdf viewer. This option is ignored and automatically set to
          ``False`` if ``filename`` is not ``None``.

        OUTPUT:

            string, path to pdf file

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: V = [[1,0,1],[1,0,0],[1,1,0],[0,0,-1],[0,1,0],[-1,0,0],[0,1,1],[0,0,1],[0,-1,0]]
            sage: P = Polyhedron(vertices=V).polar()
            sage: s = P.projection().tikz([674,108,-731],112)
            sage: t = TikzPicture(s)
            sage: _ = t.pdf()    # not tested

        ::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.pdf')
            sage: _ = t.pdf(filename)   # long time (2s)

        ACKNOWLEDGEMENT:

            The code was adapted and taken from the module :mod:`sage.misc.latex.py`.
        """
        if not have_pdflatex():
            raise RuntimeError("PDFLaTeX does not seem to be installed. " 
                    "Download it from ctan.org and try again.")

        # set up filenames
        _filename_tex = tmp_filename('tikz_','.tex')
        with open(_filename_tex, 'w') as f:
            f.write(str(self))
        base, _filename_tex = os.path.split(_filename_tex)
        _filename, ext = os.path.splitext(_filename_tex)

        # subprocess stuff
        from subprocess import check_call, CalledProcessError, PIPE

        # running pdflatex
        cmd = ['pdflatex', '-interaction=nonstopmode', _filename_tex]
        cmd = ' '.join(cmd)
        try:
            check_call(cmd, shell=True, stdout=PIPE, stderr=PIPE, cwd=base)
        except (CalledProcessError, OSError):
            _filename_log = os.path.join(base, _filename+'.log')
            if os.path.exists(_filename_log):
                with open(_filename_log) as f:
                    print(f.read())
            else:
                print("Error: log file was not found")
            raise OSError("Error when running pdflatex (see log printed above).")
        _filename_pdf = os.path.join(base, _filename+'.pdf')

        # move the pdf into the good location
        if filename:
            filename = os.path.abspath(filename)
            cmd = ['mv', _filename_pdf, filename]
            check_call(cmd, stdout=PIPE, stderr=PIPE)
            return filename

        # open the tmp pdf
        elif view:
            from sage.misc.viewer import pdf_viewer
            cmd = [pdf_viewer(), _filename_pdf]
            cmd = ' '.join(cmd)
            check_call(cmd, shell=True, cwd=base, stdout=PIPE, stderr=PIPE)

        return _filename_pdf

    def png(self, filename=None, density=150, view=True):
        """
        Compiles the latex code with pdflatex and converts to a png file.

        INPUT:

        - ``filename`` -- string (default:``None``), the output filename. 
          If ``None``, it saves the file in a temporary directory.

        - ``density`` -- integer, (default: ``150``), horizontal and vertical
          density of the image

        - ``view`` -- bool (default:``True``), whether to open the file in a
          png viewer. This option is ignored and automatically set to
          ``False`` if ``filename`` is not ``None``.

        OUTPUT:

            string, path to png file

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: V = [[1,0,1],[1,0,0],[1,1,0],[0,0,-1],[0,1,0],[-1,0,0],[0,1,1],[0,0,1],[0,-1,0]]
            sage: P = Polyhedron(vertices=V).polar()
            sage: s = P.projection().tikz([674,108,-731],112)
            sage: t = TikzPicture(s)
            sage: _ = t.png()    # not tested

        ::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.png')
            sage: _ = t.png(filename)      # long time (2s)

        ACKNOWLEDGEMENT:

            The code was adapted and taken from the module :mod:`sage.misc.latex.py`.
        """
        if not have_convert():
            raise RuntimeError("convert (from the ImageMagick suite) does not "
                  "appear to be installed. Converting PDFLaTeX output to png "
                  "requires this program, so please install and try again. "
                  "Go to http://www.imagemagick.org to download it.")

        # subprocess stuff
        from subprocess import check_call, PIPE

        _filename_pdf = self.pdf(filename=None, view=False)
        _filename, ext = os.path.splitext(_filename_pdf)
        _filename_png = _filename+'.png'

        # convert to png
        cmd = ['convert', '-density',
               '{0}x{0}'.format(density), '-trim', _filename_pdf,
               _filename_png]
        cmd = ' '.join(cmd)
        check_call(cmd, shell=True, stdout=PIPE, stderr=PIPE)

        # move the png into the good location
        if filename:
            filename = os.path.abspath(filename)
            cmd = ['mv', _filename_png, filename]
            check_call(cmd, stdout=PIPE, stderr=PIPE)
            return filename

        # open the tmp png
        elif view:
            from sage.misc.viewer import png_viewer
            cmd = [png_viewer(), _filename_png]
            cmd = ' '.join(cmd)
            check_call(cmd, shell=True, stdout=PIPE, stderr=PIPE)

        return _filename_png

    def svg(self, filename=None, view=True):
        """
        Compiles the latex code with pdflatex and converts to a svg file.

        INPUT:

        - ``filename`` -- string (default:``None``), the output filename. 
          If ``None``, it saves the file in a temporary directory.

        - ``view`` -- bool (default:``True``), whether to open the file in
          a browser. This option is ignored and automatically set to
          ``False`` if ``filename`` is not ``None``.

        OUTPUT:

            string, path to svg file

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: V = [[1,0,1],[1,0,0],[1,1,0],[0,0,-1],[0,1,0],[-1,0,0],[0,1,1],[0,0,1],[0,-1,0]]
            sage: P = Polyhedron(vertices=V).polar()
            sage: s = P.projection().tikz([674,108,-731],112)
            sage: t = TikzPicture(s)
            sage: _ = t.svg()    # not tested

        ::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.svg')
            sage: _ = t.svg(filename)      # long time (2s)

        ACKNOWLEDGEMENT:

            The code was adapted and taken from the module :mod:`sage.misc.latex.py`.
        """
        if not have_program('pdf2svg'):
            raise RuntimeError("pdf2svg does not seem to be installed. " 
                    "Install it for example with ``brew install pdf2svg``"
                    " or ``apt-get install pdf2svg``.")

        # subprocess stuff
        from subprocess import check_call, PIPE

        _filename_pdf = self.pdf(filename=None, view=False)
        _filename, ext = os.path.splitext(_filename_pdf)
        _filename_svg = _filename+'.svg'

        # convert to svg
        cmd = ['pdf2svg', _filename_pdf, _filename_svg]
        cmd = ' '.join(cmd)
        check_call(cmd, shell=True, stdout=PIPE, stderr=PIPE)

        # move the svg into the good location
        if filename:
            filename = os.path.abspath(filename)
            cmd = ['mv', _filename_svg, filename]
            check_call(cmd, stdout=PIPE, stderr=PIPE)
            return filename

        # open the tmp svg
        elif view:
            from sage.misc.viewer import browser
            cmd = [browser(), _filename_svg]
            cmd = ' '.join(cmd)
            check_call(cmd, shell=True, stdout=PIPE, stderr=PIPE)

        return _filename_svg

    def tex(self, filename=None, include_header=True):
        """
        Writes the latex code to a file.

        INPUT:

        - ``filename`` -- string (default:``None``), the output filename.
          If ``None``, it saves the file in a temporary directory.
        - ``include_header`` -- bool (default:``True``) whether to include
          the header latex part. If ``False``, it prints only the
          tikzpicture part to the file.

        OUTPUT:

            string, path to tex file

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: V = [[1,0,1],[1,0,0],[1,1,0],[0,0,-1],[0,1,0],[-1,0,0],[0,1,1],[0,0,1],[0,-1,0]]
            sage: P = Polyhedron(vertices=V).polar()
            sage: s = P.projection().tikz([674,108,-731],112)
            sage: t = TikzPicture(s)
            sage: _ = t.tex()

        Write only the tikzpicture without header and begin/end document::

            sage: _ = t.tex(include_header=False)

        Write to a given filename::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.tex')
            sage: _ = t.tex(filename)

        """
        if filename is None:
            filename = tmp_filename('tikz_','.tex')
        else:
            filename = os.path.abspath(filename)

        if include_header:
            output = str(self)
        else:
            output = self.tikz_picture_code()

        with open(filename, 'w') as f:
            f.write(output)

        return filename

