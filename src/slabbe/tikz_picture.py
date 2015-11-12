# -*- coding: utf-8 -*-
r"""
Dealing with tikzpicture

Conversions to pdf, png.

EXAMPLES::

    sage: from slabbe import TikzPicture
    sage: g = graphs.PetersenGraph()
    sage: s = latex(g)
    sage: t = TikzPicture(s, standalone_configs=["border=3mm"], packages=['tkz-graph'])
    sage: t.pdf()    # not tested

::

    sage: t
    \documentclass[tikz]{standalone}
    \standaloneconfig{border=3mm}
    \usepackage{tkz-graph}
    \begin{document}
    \begin{tikzpicture}
    %
    \useasboundingbox (0,0) rectangle (5.0cm,5.0cm);
    %
    \definecolor{cv0}{rgb}{0.0,0.0,0.0}
    ...
    ... 68 lines not printed (3748 characters in total) ...
    ...
    \Edge[lw=0.1cm,style={color=cv6v8,},](v6)(v8)
    \Edge[lw=0.1cm,style={color=cv6v9,},](v6)(v9)
    \Edge[lw=0.1cm,style={color=cv7v9,},](v7)(v9)
    %
    \end{tikzpicture}
    \end{document}

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
from sage.misc.latex import have_pdflatex
from sage.misc.latex import have_convert
from sage.misc.temporary_file import tmp_filename
from sage.structure.sage_object import SageObject
import os

class TikzPicture(SageObject):
    def __init__(self, code, standalone_configs=[], packages=[],
            tikzlibraries=[], macros=[]):
        r"""
        INPUT:

        - ``code`` -- string, tikzpicture code
        - ``standalone_configs`` -- list of strings (default: ``[]``),
          standalone configuration options.
        - ``packages`` -- list of strings or ``'sage_preamble'`` (default:
          ``[]``), latex packages. If ``'sage_preamble'``, it is replaced by the
          sage latex preamble.
        - ``tikzlibraries`` -- list of strings (default: ``[]``), tikz libraries
          to use.
        - ``macros`` -- list of strings or ``'sage_latex_macros'`` (default:
          ``[]``), stuff you need for the picture. If ``'sage_latex_macros'``,
          it is replaced by the sage latex macros.

        POSSIBLE OPTIONS:

        Here are some standalone config you may want to use (see standalone
        documentation for more options)::

            standalone_configs = ['preview', 'border=3mm', 'varwidth', 'beamer',
            'float']

        Here are some packages you may want to load::

            packages = ['nicefrac', 'amsmath', 'pifont', 'tikz-3dplot',
            'tkz-graph', 'tkz-berge', 'pgfplots']

        Here are some tikzlibraries you may want to load::

            tikzlibraries = ['arrows', 'snakes', 'backgrounds', 'patterns',
            'matrix', 'shapes', 'fit', 'calc', 'shadows', 'plotmarks',
            'positioning', 'pgfplots.groupplots', 'mindmap']

        Here are some macros you may want to set::

            macros = [r'\newcommand{\ZZ}{\mathbb{Z}}']

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: g = graphs.PetersenGraph()
            sage: s = latex(g)
            sage: t = TikzPicture(s, standalone_configs=["border=3mm"], packages=['tkz-graph'])
            sage: _ = t.pdf()    # not tested

        ::

            sage: from slabbe.matrix_cocycle import cocycles
            sage: B = cocycles.Brun()
            sage: s = B.tikz_n_cylinders(1, scale=3)
            sage: t = TikzPicture(s)
        """
        self._code = code
        self._standalone_configs = standalone_configs
        self._tikzlibraries = tikzlibraries
        if packages == 'sage_preamble':
            from sage.misc.latex import _Latex_prefs
            self._packages = _Latex_prefs._option['preamble']
        else:
            self._packages = packages
        if macros == 'sage_latex_macros':
            from sage.misc.latex_macros import sage_latex_macros
            self._macros = sage_latex_macros()
        else:
            self._macros = macros

    def _latex_file_header_lines(self):
        r"""
        """
        lines = []
        lines.append(r"\documentclass[tikz]{standalone}")
        for config in self._standalone_configs:
            lines.append(r"\standaloneconfig{{{}}}".format(config))
        for package in self._packages:
            lines.append(r"\usepackage{{{}}}".format(package))
        for library in self._tikzlibraries:
            lines.append(r"\usetikzlibrary{{{}}}".format(library))
        lines.extend(self._macros)
        return lines

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: g = graphs.PetersenGraph()
            sage: s = latex(g)
            sage: t = TikzPicture(s, packages=['tkz-graph'])
            sage: t
            \documentclass[tikz]{standalone}
            \usepackage{tkz-graph}
            \begin{document}
            \begin{tikzpicture}
            %
            \useasboundingbox (0,0) rectangle (5.0cm,5.0cm);
            %
            \definecolor{cv0}{rgb}{0.0,0.0,0.0}
            ...
            ... 68 lines not printed (3748 characters in total) ...
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
        """
        lines = self._latex_file_header_lines()
        lines.append(r"\begin{document}")
        lines.append(self._code)
        lines.append(r"\end{document}")
        return '\n'.join(lines)

    def pdf(self, filename=None, view=True):
        """
        Compiles the latex code with pdflatex and create a pdf file.

        INPUT:

        - ``filename`` -- string (default:``None``), the output filename. 
          If ``None``, a temporary file is created.

        - ``view`` -- bool (default:``True``), whether to open the file in a
          pdf viewer. This option is ignored if ``filename`` is not ``None``.

        OUTPUT:

            pdf file

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: from slabbe.matrix_cocycle import cocycles
            sage: B = cocycles.Brun()
            sage: s = B.tikz_n_cylinders(1, scale=3)
            sage: t = TikzPicture(s)
            sage: _ = t.pdf()    # not tested

        ::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.pdf')
            sage: _ = t.pdf(filename)
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
        cmd = ['sage-native-execute', 'pdflatex', '-interaction=nonstopmode', _filename_tex]
        cmd = ' '.join(cmd)
        try:
            check_call(cmd, shell=True, stdout=PIPE, stderr=PIPE, cwd=base)
        except (CalledProcessError, OSError):
            _filename_log = os.path.join(base, _filename+'.log')
            if os.path.exists(_filename_log):
                with open(_filename_log) as f:
                    print(f.read())
            else:
                print "Error: log file was not found"
            raise OSError("Error when running pdflatex (see log printed above).")
        _filename_pdf = os.path.join(base, _filename+'.pdf')

        # move the pdf into the good location
        if filename:
            filename = os.path.abspath(filename)
            cmd = ['sage-native-execute', 'mv', _filename_pdf, filename]
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
          If ``None``, it opens the file in a viewer. Otherwise, it just
          saves the file into location ``filename``.

        - ``density`` -- integer, (default: ``150``), horizontal and vertical
          density of the image

        - ``view`` -- bool (default:``True``), whether to open the file in a
          png viewer. This option is ignored if ``filename`` is not ``None``.

        OUTPUT:

            png file

        EXAMPLES::

            sage: from slabbe import TikzPicture
            sage: from slabbe.matrix_cocycle import cocycles
            sage: B = cocycles.Brun()
            sage: s = B.tikz_n_cylinders(1, scale=3)
            sage: t = TikzPicture(s)
            sage: _ = t.png()    # not tested

        ::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.png')
            sage: _ = t.png(filename)
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
        cmd = ['sage-native-execute', 'convert', '-density',
               '{0}x{0}'.format(density), '-trim', _filename_pdf,
               _filename_png]
        cmd = ' '.join(cmd)
        check_call(cmd, shell=True, stdout=PIPE, stderr=PIPE)

        # move the png into the good location
        if filename:
            filename = os.path.abspath(filename)
            cmd = ['sage-native-execute', 'mv', _filename_png, filename]
            check_call(cmd, stdout=PIPE, stderr=PIPE)
            return filename

        # open the tmp png
        elif view:
            from sage.misc.viewer import png_viewer
            cmd = [png_viewer(), _filename_png]
            cmd = ' '.join(cmd)
            check_call(cmd, shell=True, stdout=PIPE, stderr=PIPE)

        return _filename_png

