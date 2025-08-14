Installation
============


Supported Python versions
-------------------------

The support for python versions depends on pysiral versions.

+-------------------+------------------+
| pysiral version   | python version   |
+===================+==================+
| ``<0.7.3``        | ``2.7``          |
+-------------------+------------------+
| ``0.7.3-0.12``    | ``3.9``          |
+-------------------+------------------+
| ``>=0.13``        | ``3.10-3.13``    |
+-------------------+------------------+

Requirements
------------

- Python (recommended: 3.10 or later for latest version)
- C-compiler (required for compiling cython code)
- git version control


Installation from github
------------------------

pysiral is intended to be installed with python 3.10 or later directly
from github:

.. code-block:: shell

   pip install pysiral@git+https://github.com/pysiral/pysiral.git


Installation from source
------------------------

To install pysiral from source, clone the repository and run the setup script:

.. code-block:: shell

    git clone https://github.com/pysiral/pysiral.git
    cd pysiral
    pip install .


Specific Branches
-----------------

If a specific pysiral version or branch is required, it can be specified in the pip install command:

.. code-block:: shell

   pip install pysiral@git+https://github.com/pysiral/pysiral.git@production/awi-v2p6


For a full documentation of the syntax see the `pip documentation <pip_install_link>`__.

.. _pip_install_link: https://pip.pypa.io/en/stable/cli/pip_install/


.. warning::
    The branches of pysiral are used to ensure consistent product
    generation for different datasets. That means that users should
    be aware that each branch may be based on a different pysiral
    and will have different features and functionality.


Optional Dependencies
---------------------

pysiral has some optional dependencies that are not installed by default. 
There are two categories of optional

- ``dev``: python packages used for development, e.g. testing, linting, debugging and building of the documentation. 
- ``ml``: python packages used for using machine learning models (e.g. `pytorch` and `xgboost`).

The optional dependencies can be installed by adding the corresponding extras when installing pysiral.
In the examples from above, installing all optional dependencies would look like this:

.. code-block:: shell

   pip install pysiral[dev,ml]@git+https://github.com/pysiral/pysiral.git

respectively:

.. code-block:: shell

    git clone https://github.com/pysiral/pysiral.git
    cd pysiral
    pip install .[dev,ml]


Local Code
----------

For development purpose it may be useful to install pysiral in editable mode. 
This step requires a manual compilation of Cython code (line 4) and is recommended to install
pysiral with all optional dependencies (line 3) and also running the test suite (line 5)

.. code-block:: shell
   :linenos:
   :emphasize-lines: 3,4,5

   git clone https://github.com/pysiral/pysiral.git
   cd pysiral
   pip install -e .[dev,ml]
   python setup.py build_ext --inplace
   python -m unittest discover -s tests -t tests





