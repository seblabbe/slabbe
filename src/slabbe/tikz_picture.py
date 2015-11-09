# -*- coding: utf-8 -*-
r"""
Dealing with tikzpicture

EXAMPLES::

    sage: g = graphs.PetersenGraph()
    sage: s = latex(g)
    sage: t = TikzPicture(s, ['tkz-graph'])
    sage: t.pdf()

::

    sage: t
    \documentclass[tikz]{standalone}
    \usepackage{tikz}
    \usepackage{tkz-graph}
    \begin{document}
    \begin{tikzpicture}
    %
    \useasboundingbox (0,0) rectangle (5.0cm,5.0cm);
    %
    \definecolor{cv0}{rgb}{0.0,0.0,0.0}
    ...
    ... 68 lines not printed (3748 characters) ...
    ...
    \Edge[lw=0.1cm,style={color=cv6v8,},](v6)(v8)
    \Edge[lw=0.1cm,style={color=cv6v9,},](v6)(v9)
    \Edge[lw=0.1cm,style={color=cv7v9,},](v7)(v9)
    %
    \end{tikzpicture}
    \end{document}

.. TODO::

    - png output
    - larger crop (Petersen Graph)
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

class TikzPicture(SageObject):
    def __init__(self, code, packages=[], tikzlibraries=[], commands=[]):
        r"""
        INPUT:

        - ``code`` -- string, tikzpicture code
        - ``packages`` -- list or ``'sage_preamble'`` (default: ``[]``),
          latex packages. If ``'sage_preamble'``, it is replaced by the
          sage latex preamble.
        - ``tikzlibraries`` -- list (default: ``[]``), tikz libraries to
          use.
        - ``commands`` -- list or ``'sage_latex_macros'`` (default: ``[]``), 
          stuff you need for the picture. If ``'sage_latex_macros'``, it is
          replaced by the sage latex macros.

        Here are some packages you may want to load:

            packages = ['nicefrac', 'amsmath', 'pifont', 'tikz-3dplot',
            'tkz-graph', 'tkz-berge', 'pgfplots']

        Here are some tikzlibraries you may want to load:

            tikzlibraries = ['arrows', 'snakes', 'backgrounds', 'patterns',
            'matrix', 'shapes', 'fit', 'calc', 'shadows', 'plotmarks',
            'positioning', 'pgfplots.groupplots', 'mindmap']

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: s = latex(g)
            sage: t = TikzPicture(s, ['tkz-graph'])
            sage: t.pdf()

        ::

            sage: from slabbe.matrix_cocycle import cocycles
            sage: B = cocycles.Brun()
            sage: s = B.tikz_n_cylinders(1, scale=3)
            sage: t = TikzPicture(s)

        """
        self._code = code
        if packages == 'sage_preamble':
            from sage.misc.latex import _Latex_prefs
            self._packages = _Latex_prefs._option['preamble']
        else:
            self._packages = packages
        self._tikzlibraries = tikzlibraries
        if commands == 'sage_latex_macros':
            from sage.misc.latex_macros import sage_latex_macros
            self._commands = sage_latex_macros()
        else:
            self._commands = commands

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: s = latex(g)
            sage: t = TikzPicture(s, ['tikz', 'tkz-graph'])
            sage: t
            \documentclass{standalone}
            \usepackage{tikz}
            \usepackage{tkz-graph}
            \begin{document}
            \begin{tikzpicture}
            %
            \useasboundingbox (0,0) rectangle (5.0cm,5.0cm);
            %
            \definecolor{cv0}{rgb}{0.0,0.0,0.0}
            ...
            ... 68 lines not printed (3748 characters total) ...
            ...
            \Edge[lw=0.1cm,style={color=cv6v8,},](v6)(v8)
            \Edge[lw=0.1cm,style={color=cv6v9,},](v6)(v9)
            \Edge[lw=0.1cm,style={color=cv7v9,},](v7)(v9)
            %
            \end{tikzpicture}
            \end{document}

        """
        lines = self._latex_file_head_lines()
        L = self._code.splitlines()
        if len(L) <= 10:
            lines.extend(L)
        else:
            lines.extend(L[:5])
            lines.append('...')
            lines.append('... {} lines not printed ({} characters total) ...'.format(len(L)-10, 
                                                           len(self._code)))
            lines.append('...')
            lines.extend(L[-5:])
        lines.append(r"\end{document}")
        return '\n'.join(lines)

    def latex_file(self):
        r"""
        """
        lines = self._latex_file_head_lines()
        lines.append(self._code)
        lines.append(r"\end{document}")
        return '\n'.join(lines)

    def _latex_file_head_lines(self):
        r"""
        """
        lines = []
        lines.append(r"\documentclass[tikz]{standalone}")
        for package in self._packages:
            lines.append(r"\usepackage{{{}}}".format(package))
        for library in self._tikzlibraries:
            lines.append(r"\usetikzlibrary{{{}}}".format(library))
        lines.extend(self._commands)
        lines.append(r"\begin{document}")
        return lines

    def pdf(self, filename=None):
        """
        Compiles the latex code with pdflatex and create a pdf file.

        INPUT:

        - ``filename`` -- string (default:``None``), the output filename. 
          If ``None``, it opens the file in a viewer. Otherwise, it just
          saves the file into location ``filename``.

        OUTPUT:

            pdf file

        EXAMPLES::

            sage: from slabbe.matrix_cocycle import cocycles
            sage: B = cocycles.Brun()
            sage: s = B.tikz_n_cylinders(1, scale=3)
            sage: t = TikzPicture(s)
            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.pdf')
            sage: t.pdf(filename)   # random
            Creation of file /Users/slabbe/.sage/temp/Computer.local/368/tempgHgyRj.pdf

        """
        if not have_pdflatex():
            raise RuntimeError("PDFLaTeX does not seem to be installed. " 
                    "Download it from ctan.org and try again.")

        # set up filenames
        _filename_tex = tmp_filename('tikz_','.tex')
        with open(_filename_tex, 'w') as f:
            f.write(self.latex_file())
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

        # open or move the pdf into the good location
        if filename is None:
            from sage.misc.viewer import pdf_viewer
            viewer = pdf_viewer()
            cmd = [viewer, base+'/'+_filename+'.pdf']
            cmd = ' '.join(cmd)
            check_call(cmd, shell=True, cwd=base, stdout=PIPE, stderr=PIPE)
        else:
            filename = os.path.abspath(filename)
            cmd = ['sage-native-execute', 'mv', _filename+'.pdf', filename]
            check_call(cmd, stdout=PIPE, stderr=PIPE, cwd=base)
            print "Creation of file {}".format(filename)

    def png(self, filename=None):
        """
        Compiles the latex code with pdflatex and converts to a png file.

        INPUT:

        - ``filename`` -- string (default:``None``), the output filename. 
          If ``None``, it opens the file in a viewer. Otherwise, it just
          saves the file into location ``filename``.

        OUTPUT:

            png file

        EXAMPLES::

            sage: from slabbe.matrix_cocycle import cocycles
            sage: B = cocycles.Brun()
            sage: s = B.tikz_n_cylinders(1, scale=3)
            sage: t = TikzPicture(s)
            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.png')
            sage: t.png(filename)   # random
            Creation of file /Users/slabbe/.sage/temp/Computer.local/368/tempgHgyRj.pdf

        """
        if not have_pdflatex():
            raise RuntimeError("PDFLaTeX does not seem to be installed. " 
                    "Download it from ctan.org and try again.")
        if not have_convert():
            raise RuntimeError("convert (from the ImageMagick suite) does not "
                  "appear to be installed. Converting PDFLaTeX output to png "
                  "requires this program, so please install and try again. "
                  "Go to http://www.imagemagick.org to download it.")

        # set up filenames
        _filename_tex = tmp_filename('tikz_','.tex')
        with open(_filename_tex, 'w') as f:
            f.write(self.latex_file())
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

        # convert to png
        density=150
        density = int(1.4 * density / 1.3)
        cmd = ['sage-native-execute', 'convert', '-density',
               '{0}x{0}'.format(density), '-trim', _filename + '.pdf',
               filename + '.png']
        cmd = ' '.join(cmd)
        check_call(cmd, shell=True, stdout=PIPE, stderr=PIPE, cwd=base)

        # open or move the pdf into the good location
        if filename is None:
            from sage.misc.viewer import pdf_viewer
            viewer = pdf_viewer()
            cmd = [viewer, base+'/'+_filename+'.png']
            cmd = ' '.join(cmd)
            check_call(cmd, shell=True, cwd=base, stdout=PIPE, stderr=PIPE)
        else:
            filename = os.path.abspath(filename)
            cmd = ['sage-native-execute', 'mv', _filename+'.png', filename]
            check_call(cmd, stdout=PIPE, stderr=PIPE, cwd=base)
            print "Creation of file {}".format(filename)



