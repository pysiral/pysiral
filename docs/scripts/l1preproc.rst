Level-1 Pre-Processor (l1preproc)
=================================

The Level-1 Pre-Processor (l1preproc) scripts provides command-line access for 
pysiral functionality to generate pre-processed Level-1 data required for
the geophysical retrieval of sea ice thickness from radar altimeter data in
the Level-2 processor. 

Required input is the id oder file location of the Level-1 pre-processor definition, 
the selected time period, and the specific input data version since the processor
definition may be applicable to multiple input data versions.

.. seealso:: 
    - :ref:`Level-1 Pre-Processor` documentation
    - l1preproc code: :py:mod:`pysiral.scripts.l1prepoc`


Usage
-----

The Level-1 Pre-Processor can be executed via python or via the command-line
interface.

.. code-block:: bash

    python -m pysiral.scripts l1preproc --help

or

.. code-block:: bash

    pysiral(.exe) l1preproc --help


**Output**

.. code-block:: 

        usage: pysiral l1preproc [-h] [-E <month number> [<month number> ...]] [-H nh|sh|global]
                            [-p <platform_id>] [-s <source_dataset_id>]
                            [--multiprocessing | --no-multiprocessing] [-m <num_cores>]
                            <id|filepath> <processing period>

    The Level-1 Pre-Processor (l1preproc) is used to generate Level-1 files (l1p) from individual
    source radar altimeter files for a given period. Processing steps include the harmonization of
    data formats from various radar altimeter missions, the generation of continuous trajectories
    over the polar oceans, ingesting auxiliary data and the pre-computation of waveform shape
    parameters. Basis for the Level-1 Pre-Processor is a Level-1 Processor Definition file that
    defines the the radar altimeter input and processing parameters. The output is a set of Level-1
    files (l1p) that contain the harmonized data from the source files, which can be used for
    further processing in the Level-2 processor (l2proc).

    positional arguments:
        <id|filepath>       Identifier or file path to the Level-1 Pre-Processor definition file.
                            This file contains the settings for the Level-1 processor. The default
                            location for these files is `{pysiral-cfg-location}/proc/l1/`. The
                            identifier is the filename without the `.yaml` extension.
                            E.g.`cryosat2_pds_ipf1e_v1p2` will be resolved to `{pysiral-cfg-
                            location}/proc/l1/cryosat2_pds_ipf1e_v1p2.yaml`.
        <processing period>
                            Period definition for processing, given as a string in the format "YYYY-
                            MM[-DD][:YYYY-MM[-DD]]". If only one date is given, it will be
                            interpreted as a period (e.g., "2023-01" for January 2023 and
                            ("2023-01-01" for one day). If two colon-separated dates are given, they
                            will be interpreted as a start and end date or month.

    options:
        -h, --help          show this help message and exit
        -E <month number> [<month number> ...], --exclude-months <month number> [<month number> ...]
                            List of months to be excluded from processing, given as integers (1-12).
        -H nh|sh|global, --hemisphere nh|sh|global
                            Target hemisphere for processing. Options are 'global', 'nh', or 'sh'.
                            The latitude limit of the hemisphere is defined in the Level-1 pre-
                            processor settings file. If 'global' is selected, the processor will run
                            for both hemispheres, but still within the latitude limits.
        -p <platform_id>, --platform-id <platform_id>
                            Radar altimeter platform id as defined in pysiral (see `pysiral info
                            --platforms`). This option is required only if the processor
                            configuration file is applicable to multiple platforms (e.g. sentinel3a,
                            sentinel3b, etc.). If the processor configuration file is platform-
                            specific, this option has no effect.
        -s <source_dataset_id>, --source-dataset-id <source_dataset_id>
                            Identifier of the source dataset to be used for processing, summarizing
                            the platform, version, and timeliness information. The source dataset ID
                            must be specified in the local machine definition file ({pysiral-cfg-
                            location}/local_machine_def.yaml, e.g.
                            `root.l1b_repository.<platform>.<source_dataset_id>`)
        --multiprocessing, --no-multiprocessing
                            Flag to allow disabling multiprocessing. If set, the processor will run
                            in single-threaded mode. Default is True, meaning multiprocessing is
                            enabled. (also see option -m/--multiprocessing-num-cores)
        -m <num_cores>, --multiprocessing-num-cores <num_cores>
                            Set the number of CPU cores to be used for multiprocessing. If not set,
                            the default value (derived from `multiprocessing.cpu_count()`) will be
                            used. NOTE: In some managed environments, the default value is not
                            reliable, which may lead to performance issues. In this case, it is
                            recommended to set this value manually to the known number of available
                            CPU cores.

    For more information, see: https://pysiral.readthedocs.io


.. tip:: 

    Future versions of pysiral will support additional options to list and query
    specific processor settings and their dependencies via
    
    ``pysiral info --target l1p-settings``
