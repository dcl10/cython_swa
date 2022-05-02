from setuptools import setup
from Cython.Build import cythonize

setup(
    name="cython_swa",
    packages=["cython_swa"],
    ext_modules=cythonize("cython_swa/swa.py"),
    zip_safe=False,
)