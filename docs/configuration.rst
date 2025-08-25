Configuration
=============

Several aspects needs configuration in pysiral, including:

1. Location of the input data files for the local environment
2. Location of the directory for output products
3. Properties of the various auxiliary datasets
4. Definition of the pysiral processors 
5. Definition of the various output formats
6. Definition of the grids used in the Level-3 processor
7. The configuration of parallel processing 
8. Location of the required configuration files for the points above. 

Most of these aspects are specified in configuration files. The standard format
in pysiral is YAML, due to its machine and human readability. 

All configuration files with the exception of the local data configuration file are 
located in the ``resources/pysiral-cfg`` directory of the pysiral package. 

Local Data
----------

Data storage is unique to each user's environment and must be configured accordingly. 
The local data configuration file specifies the paths to the input data files 
(both radar altimeter and auxiliary dataset) and the output directory for the processed products.

The local data configuration file is named ``local_machine_def.yaml`` and is located as default in the 
``.pysiral-cfg`` directory in the user's home directory (``$HOME`` `for Windows Powershell/Linux and` 
``%userprofile%`` `for Windows Command Prompt`). 

The pysiral code base includes a 
`template configuration <https://github.com/pysiral/pysiral/blob/main/pysiral/resources/pysiral-cfg/templates/local_machine_def.yaml>`__
file in the ``resources/pysiral-cfg/template`` directory that users can modify for their local setup. 
Each entry needs to specify the root directory for the input files, which partly need to be stored
in a pre-defined sub-directory structure depending on the type of data. 


Processor Definition
--------------------

The processor definitions files are located in the ``resources/pysiral-cfg/proc/<processing_level>`` directory and are 
named ``<processor_id>.yaml``. These files are meant to be under revision control and should not be
modified directly. Instead, users should create their own copies of these files in a separate directory
and modify the copies for their local setup. The ``<processor_id>`` needs to be a unique identifier, to allow
automatic file location lookup in the code. 

The convention for the ``<processor_id>`` depends on the specific processing level but generally follows a 
hierarchical structure that contains main product properties, e.g. 

.. code-block::

    <product-line>_<platform>_<hemisphere>_<version>_<data-record-type>.yaml

(e.g.: ``cci_cryosat2_nh_v4p0_cdr.yaml``).


Output Definition
-----------------

The output definition files are located in the ``resources/pysiral-cfg/output/<processing_level>`` directory and are 
named according to the specific output format they define. Simular to the processor definitions files, users should 
create their own versions of these files in a separate directory and modify the copies for their local setup.

These files include the global attributes, variables and their attributes according to the 
`Climate & Forecast netCDF conventions <https://cfconventions.org>`__. pysiral uses as extension of the `.format()`
python function to fill defined variables in the configuration files. This is done that for example only
one output definition file is required for a certain product and the attributes for different radar altimeter
platforms of hemispheres are filled dynamically.


Auxiliary Dataset Definition
----------------------------

.. attention:: 

    The ``auxdata_def.yaml`` will be removed in future versions of pysiral. It was designed to hold 
    class configuration data that does not need to be duplicated in the various Level-2 processor 
    definition files. Its content will be moved to `pydantic` based configuration models that will 
    also provide means of configuration data validation. 


Custom Location of Configuration Files
--------------------------------------

As a default, the configuration files are located in the python package (either in the installed package or in the 
source code directory), with the exception of the `local_machine_def.yaml` file, which is located in the users home directory.

In order to allow running multiple pysiral installations with different configurations on the same platform, the 
configuration location can be specified per pysiral installation. 

There are two ways of specifying the location of the configuration files:

1. By modifying the file `PYSIRAL-CFG-LOC` in the main pysiral package directory. By default the file contains only the
   string ``PACKAGE`` indicating that the configuration files are located in the python package directory and the
   `local_machine_def.yaml` file is located in the user's home directory. The string `PACKAGE` can be
   replaced with an absolute path to a directory containing the configuration files. In this case also the
   `local_machine_def.yaml` file needs to be located in the specified directory. The content of the `pysiral-cfg`
   then need to be copied over to the specified directory.

2. A pysiral script is available to automate the copying of the configuration files to a new location:
   ``python "${env_path}"/bin/pysiral_cfg_setdir.py -create-config-dir "<target-directory>"``. Please note 
   that the script will not create a `local_machine_def.yaml` file. 


.. note::

    The scripts managing the configuration files will change in upcoming versions of pysiral, allowing
    a more inuitive and file control of the pysiral configuration.

Multiprocessing
---------------


