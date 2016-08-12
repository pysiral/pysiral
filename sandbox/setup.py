from Cython.Build import cythonize
from setuptools.extension import Extension
from setuptools import setup

extensions = [
    Extension(
        "cyfindpeaks",
        ["cyfindpeaks.pyx"])]

setup(
    ext_modules = cythonize(extensions)
)