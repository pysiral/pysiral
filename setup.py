# -*- coding: utf-8 -*-

import re
from pathlib import Path

import numpy
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from setuptools import find_packages, setup
from setuptools.extension import Extension

# Get the readme
readme = Path('README.md').read_text()
# Get the licence
license_text = Path('LICENSE').read_text()
# Get the version
mypackage_root_dir = Path(__file__).absolute().parent
version_file_path = mypackage_root_dir / "pysiral" / 'VERSION'
with open(str(version_file_path)) as version_file:
    version = version_file.read().strip()

# cythonized extensions go here
# TODO: Autodetect cython files
extensions = [
    Extension("pysiral.retracker.cytfmra", ["pysiral/retracker/cytfmra.pyx"])
]

# Package requirements
with open("requirements.txt") as f:
    requirements_content = f.read().splitlines()
install_requires = [r for r in requirements_content if r.find("git+")]
dependency_links = [r for r in requirements_content if not r.find("git+")]


# NOTE: links to git repositories cause the built_ext process to crash.
# -> Thus we sanitize the requirements string here.
for i, requirement in enumerate(requirements_content):
    m = re.search(r'/(.+?)\.git', requirement)
    if m:
        package_name = m[1].split("/")[-1]
        requirements_content[i] = package_name

setup(
    name='pysiral',
    version=version,
    description='PYthon Sea Ice Radar ALtimetry toolbox',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    long_description=readme,
    author='Stefan Hendricks',
    author_email='stefan.hendricks@awi.de',
    url='https://github.com/pysiral/pysiral',
    license=license_text,
    install_requires=install_requires,
    dependency_links=dependency_links,
    packages=find_packages(exclude=('tests', 'docs')),
    include_package_data=True,
    scripts=['bin/pysiral_cfg_update.py',
             'bin/pysiral_cfg_setdir.py',
             'bin/pysiral-l1preproc.py',
             'bin/pysiral-l2proc.py',
             'bin/pysiral-l2preproc.py',
             'bin/pysiral-l3proc.py'],
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(
        extensions,
        compiler_directives={'language_level': "3"}
    ),
    include_dirs=[numpy.get_include()]
)
