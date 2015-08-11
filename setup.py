from distutils.core import setup
from Cython.Build import cythonize
import numpy as np

setup(
    name = 'Superdroplet wrapper',
    ext_modules = cythonize("sd_cy.pyx")
)
