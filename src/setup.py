from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import Cython.Compiler.Options
import os

Cython.Compiler.Options.old_style_globals = True

setup(name='slabbe',
	version='0.1',
	description="Sebastien Labbe's Research code",
	author='Sebastien Labbe',
	author_email='labbe@liafa.univ-paris-diderot.fr',
	url='http://www.liafa.univ-paris-diderot.fr/~labbe/',
    license = "GPL v2",
	packages=['slabbe', 
              'slabbe.digital_geometry'
             ],
    #ext_modules=[
    #Extension('slabbe.mcf',
    #sources = [os.path.join('slabbe','mcf.pyx')]),
    #], 
    #cmdclass = {'build_ext': build_ext}
)

