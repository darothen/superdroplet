"""
Following http://stackoverflow.com/questions/28301931/how-to-profile-cython-functions-line-by-line
and parcel model setup to build a line_profiler enabled module.

"""

from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize
from Cython.Distutils import build_ext

# from Cython.Compiler.Options import directive_defaults
# directive_defaults['linetrace'] = True
# directive_defaults['binding'] = True

import numpy as np

extensions = [
    Extension("*", ["*.pyx", ],
              include_dirs=[np.get_include(), ".", ],
              define_macros=[
                  ('CYTHON_TRACE', '1'),
                  ('CYTHON_TRACE_NOGIL', '1'),
              ])
]

setup(
    name = 'Superdroplet wrapper',
    ext_modules = cythonize(extensions),
    cmdclass = {'build_ext': build_ext, },
)

# setup(
#     name = 'Superdroplet wrapper',
#     ext_modules = cythonize("sd_cy.pyx")
# )
