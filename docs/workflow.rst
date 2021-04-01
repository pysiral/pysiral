Processing Workflow
===================

General workflows
-----------------

TBD

Cryo-TEMPO workflow
-------------------

.. mermaid::

   graph TD
        classDef auxdata fill:#ccc
        id1[(CS2 L1b)] --> id2[[L1p processing]]
        id1a[(FES2014b)]:::auxdata -.-> id2
        id2 --> id3[(L1p)]
        id3 --> id4[[L2 sea ice processing]]
        id3 --> id5[[L2 polar ocean processing]]
        id8[(AUX)]:::auxdata -.-> id4
        id8 -.-> id5
        id4 --> id6[(L2 SI)]
        id5 --> id7[(L2 PO)]

The Cryo-TEMPO processing firstly processes CS2 L1b data to pysiral L1p data, replacing the tide information with FES2014b as part of the processing. The processing then splits. A specific L2 processor is used for each of the sea ice and polar ocean themes, and produces a different L2 product for each theme. 

The commands used are:

* | L1p : ``pysiral-l1preproc.py -l1p-settings cryosat2_pds_ipf1d_fes2014 -source-repo-id baseline_d``
  | ``-start 2014 04 01 -stop 2014 04 01``
* | L2 SI : ``pysiral-l2proc.py -l2-settings esa_cryosat2_cryotempo_si -l2-output l2i_esa_cryotempo_si``
  | ``-input-version baseline_d -start 2014 04 01 -stop 2014 04 01``
* | L2 PO : ``pysiral-l2proc.py -l2-settings esa_cryosat2_cryotempo_po -l2-output l2i_esa_cryotempo_po``
  | ``-input-version baseline_d -start 2014 04 01 -stop 2014 04 01``
