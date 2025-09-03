pysiral Command Line Scripts
=============================

.. note:: 
    This section describes the usage of pysiral scripts
    for pysiral version v0.12 or later.

pysiral provides several command-line scripts for access to
the different processors. These are available after :ref:`installation<SectInstallation>` of pysiral in 
the python environment:

.. code-block:: bash

    pysiral <script_name> <args>

For local versions of pysiral, the scripts can be run directly from the source directory with python

.. code-block:: bash

    export PYTHONPATH="${PYTHONPATH}:<pysriral src dir>"
    python -m pysiral.scripts <script_name> <args>

or imported directly in python:

.. code-block:: python

    # from pysiral.scripts.<script_name> import <script_name>
    from pysiral.scripts.l1preproc import l1preproc

All scripts have implemented a ``-h/--help`` option that will
display a documentation of the scripts scope as well as 
required and optional input arguments. 

.. versionchanged:: v0.12

    The scripts have been updated to improve usability and 
    ease of access. Most notably, the scripts are now 
    provided as executables with improved documentation
    and more common use of short and long flags. 

.. seealso:: 

    - Scripts Module: :py:mod:`pysiral.scripts`.




