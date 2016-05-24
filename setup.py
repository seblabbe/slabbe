from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import Cython.Compiler.Options
import os
from sage.env import sage_include_directories

Cython.Compiler.Options.old_style_globals = True

ext_modules = [        Extension('slabbe.kolakoski_word_pyx',
            sources = [os.path.join('slabbe','kolakoski_word_pyx.pyx')],
            ),
        Extension('slabbe.mult_cont_frac',
            sources = [os.path.join('slabbe','mult_cont_frac.pyx')],
            include_dirs=sage_include_directories())]

setup(name='slabbe',
    version=open("VERSION").read().strip(),
    description="Sebastien Labbe's Research code",
    long_description=open('README.rst').read(),
    author='Sebastien Labbe',
    author_email='slabbe@ulg.ac.be',
    install_requires=['cython','cysignals'],
    url='http://www.slabbe.org/',
    license = "GPL v2",
    packages=['slabbe'],
    ext_modules=cythonize(ext_modules),
)

