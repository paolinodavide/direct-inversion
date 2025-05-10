from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("gr_iteration.pyx", language_level = "3", annotate = True)
)

