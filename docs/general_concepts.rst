General Concepts 
================

Processing Levels
-----------------


Key Terms
---------

- **Product Line**: A product line refers to the source of an processing
  algorithm spanning several processing levels. 


Dataset Identifiers
-------------------

Within pysiral, each dataset, either primary radar altimeter data, or
auxiliary dataset has a identifier string that is used to select 
appropriate code classes and data files. 

.. attention:: 
    The strict use of dataset identifiers will be gradually introduced
    to pysiral in the coming versions. The feature is intended to 
    simplify the workflow and scripts calls. The description in this section
    describes the target definition for dataset identifiers.

Primary Radar Altimeter Datasets
++++++++++++++++++++++++++++++++


    ``ra_source:{platform_id}:{origin}:{timeliness}:{version_or_baseline}[:{qualifier}]``

Example:

    ``ra_source:cryosat2:esa_pds:nrt:E001``


Auxiliary Dataset Identifier
++++++++++++++++++++++++++++

    ``aux:{auxiliary_type}:{dataset_id}:{timeliness}:{version_or_baseline}[:{qualifier}]``

Example

   ``aux:siconc:osi408:nrt:B``


Level-1P Dataset Identifier
+++++++++++++++++++++++++++

    ``l1p:{platform_id}:{origin}:{timeliness}:{hemisphere}:{version}[:{qualifier}]``


Example

    ``l1p:cryosat2:esa-pds:nrt:nh:v1p3``


Product Identifier
++++++++++++++++++

    ``{processing_level}:{product_line}:{platforms_of_missions}:{timeliness}:{hemisphere}:{version}[:{qualifier}]``

Examples

    - ``l2i:awi:cryosat2:nrt:nh:v2p6``
    - ``l3s:awi:cryosat2_sentinel3a_sentinel3b:rep:sh:v2p6``