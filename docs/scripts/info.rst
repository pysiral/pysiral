Pysiral Information (info)
==========================

Information about the pysiral version and the current
location of the configuration files can be obtained
via: 

.. code-block:: 

    python -m pysiral.scripts info

respectively 

.. code-block:: 

    pysiral(.exe) info


The output contains the following information:

.. code-block:: 

    Python sea ice radar altimeter toolbox (pysiral, version: {pysiral_package_version})

         repository: {origin url}
      documentation: https://pysiral.readthedocs.io
        git version: {commit hash} on branch {branch name} (source: {git_info_source})
      configuration: {package_config_dir} ({config_dir_target)
                     {local_machine_def_file_pat}
 python interpreter: {python_executable}

with the following parameters

- ``pysiral_package_version``: The pysiral package version (`pysiral.__version__`)
- ``origin_url``: The repository from which the package was obtained (e.g. http://github.com/pysiral/pysiral, (`pysiral.__git_origin__`)
- ``commit_hash``: The commit hash of the current version (`pysiral.__git_version__`).
- ``branch_name``: The name of the current branch (`pysiral.__git_branch__`)
- ``git_info_source``: The source of the git information (e.g. `local git` or `package`, see note below)
- ``package_config_dir``: The directory containing the processer and output configuration files
- ``local_machine_def_file_pat``: The path to the local machine definition file
- ``python_executable``: The path to the Python interpreter

.. tip:: 
    For traceabiliy of climate data records, pysiral contains the 
    variables ``__git_version__``, ``__git_branch__`` and ``__git_origin__``
    that, if documented in the product files, allow the determination of the
    exact software version used. If the code is excuted from within a git repository, the content of these
    variables is determined at runtime. This is not possible for python package
    installations. therefore the pysiral package includes a small configuration file 
    (`git_version.toml`) that is written by a github autocommit action for
    pysiral branches used in production. 

.. note:: 
    It is planned to allow querying specific processor settings
    and there dependencies in future versions.
