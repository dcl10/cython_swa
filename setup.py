from setuptools import setup
from Cython.Build import cythonize

setup(
    name="SWA app",
    ext_modules=cythonize("*.pyx"),
)