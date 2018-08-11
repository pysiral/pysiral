# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='pysiral',
    version='0.6.1',
    description='python sea ice radar altimetry processing library',
    long_description=readme,
    author='Stefan Hendricks',
    author_email='stefan.hendricks@awi.de',
    url='https://github.com/shendric/pysiral',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    scripts=['bin/pysiral-l1bpreproc', 'bin/pysiral-l2proc', 'bin/pysiral-l2preproc', 'bin/pysiral-l3proc'],
)