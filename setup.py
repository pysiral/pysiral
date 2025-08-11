# -*- coding: utf-8 -*-

import numpy
from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext


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
    include_dirs=[numpy.get_include()]
)
