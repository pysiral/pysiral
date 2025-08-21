What is pysiral?
================

pysiral is the python sea ice radar altimetry toolbox. It has been intially 
developed as a unified processing toolbox for estimating sea ice thickness
from Envisat and CryoSat-2 satellite radar altimeter data for the ESA Climate
Change Initiative on Sea Ice in 2018. 

Since then the functionality and support for additional radar altimeter 
platforms has been extended. The core concept of pysiral is to unify the 
various input formats from radar altimeter data providers and to define 
processing workflows for higher processing levels (Level-2 trajectories
and Level-3 spatio-temporal grids) based on procesor configuration 
files. The design of this approach is based on the needs for both the 
flexibility of prototyping and stability for an operational 
environment. 

Pysiral is currenlty used operationally for the generation of several 
satellite remote sensing products by different institutions: 

- `ESA Climate Change Initiative - Sea Ice Thickness Climate Data Record <https://climate.esa.int/en/projects/sea-ice/>`_
- `ESA CryoSat-2 thematic sea ice and ocean prodcuts (CryoTEMPO) <http://cryosat.mssl.ucl.ac.uk/tempo/>`_
- `Copernicus Climate Change Services - Sea Ice Thickness (Interim) Climate Data Record <https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-sea-ice-thickness>`_
- `AWI CryoSat-2 sea ice thickness product <https://spaces.awi.de/display/SIRAL/Sea+Ice+Thickness+from+Satellite+Radar+Altimetry>`_
  feeding into the `ESA SMOS - CryoSat-2 L4 sea ice thickness data prodcut <https://earth.esa.int/eogateway/catalog/smos-cryosat-l4-sea-ice-thickness>`_

