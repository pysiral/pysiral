.. _Auxiliary_Datasets:

Auxiliary Datasets
==================


The configuration file ``auxdata_def.yaml`` defines the auxiliary datasets used 
mainly in the Level-2 processor. 

.. _Auxiliary_ID:

Each dataset is associated with an **auxiliary dataset ID**: ``<category>:<dataset_id>``

This ID is used to link the register the auxiliary dataset in 
the configuration file ``auxdata_def.yaml`` (main pysiral configuration directory, 
`example <https://github.com/pysiral/pysiral/blob/main/pysiral/resources/pysiral-cfg/auxdata_def.yaml>`_).
with the data location specification in the ``local_machine_def.yaml``
(main pysiral configuration directory, 
`example <https://github.com/pysiral/pysiral/blob/main/pysiral/resources/pysiral-cfg/templates/local_machine_def.yaml>`_).
and the usage in a Level-2 processor configuration file.

+----------------------------------------+--------------------------------------------------------+
| Location                               |  Yaml Tag                                              |
+========================================+========================================================+
| ``auxdata_def.yaml``                   | ``root.<category>.<dataset_id>``                       |
+----------------------------------------+--------------------------------------------------------+
| ``local_machine_def.yaml`` (optional)  | ``root.auxdata_repository.<category>.<dataset_id>``    |
+----------------------------------------+--------------------------------------------------------+
| Level-2 processor definition           | ``root.auxdata.<category>.<dataset_id>``               |
+----------------------------------------+--------------------------------------------------------+

.. note:: 
   An auxiliary dataset only needs an entry in ``local_machine_def.yaml`` 
   if it depends on actual data files.

The full list of supported auxiliary datasets is available in this section. 


.. toctree::
   :maxdepth: 1
   :caption: Dataset Categories

   icecharts
   mdt
   ml
   mss
   region
   sic
   sitype
   snow




