from Cython.Build import cythonize
from setuptools.extension import Extension
from setuptools import setup
import numpy

extensions = [
    Extension(
        "cytfmra",
        ["cytfmra.pyx"])]

setup(
    name = "cytfmra",
    version = "0.0.1",
    author = "Stefan Hendricks",
    author_email = "stefan.hendricks@awi.de",
    ext_modules = cythonize(extensions),
    include_dirs = [numpy.get_include()]
)