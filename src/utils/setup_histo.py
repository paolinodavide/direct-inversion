from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("gr_function.pyx", language_level = "3", annotate = True)
)

