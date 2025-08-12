# -*- coding: utf-8 -*-

from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext


def include_numpy():
    import numpy as np
    try:
        numpy_include = np.get_include()
    except AttributeError:
        numpy_include = np.get_numpy_include()
    return numpy_include


setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(
        [
            Extension(
                "pysiral.retracker.cytfmra",
                ["src/pysiral/retracker/cytfmra.pyx"]
            )
        ],
        compiler_directives={'language_level': "3"}
    ),
    include_dirs=[include_numpy()]
)
