Snow on Sea Ice (snow)
======================


Warren99 & AMSR2 Monthly Climatology (clim_w99amsr2)
----------------------------------------------------

The monthly Warren99 and AMSR2 snow climatology provides monthly snow depth and density values
for the Arctic winter months from October through April. 

.. note:: 

    This auxiliary dataset only provives snow information for sea ice in the Northern Hemisphere.


Configuration Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

The :term:`Auxiliary Dataset ID` is ``snow:clim_w99amsr2`` and the configuration 
in the auxiliary data definition of the Level-2 processor definition file is:


.. code-block:: yaml

    - snow:
        name: clim_w99amsr2
        options:
            daily_scaling: True
            fyi_correction_factor: 0.5
            exception_on_error: False

**Options**

- ``daily_scaling``: Flag, when set to `True`, will scale the monthly climatology to daily values by
  linear interpolation of the monthly values based on the current day. The monthly mean values are 
  defined as the climatology on the 15th day of the respective month, with the exception of October (1st day)
  and April (last day). 
- ``fyi_correction_factor``: Factor to correct the snow depth values for :term:`FYI` (requires `sea_ice_type`
  as input. E.g. a factor of `0.5` halves the climatological snow depth values.
- ``exception_on_error``: Flag, when set to `True`, will raise an exception if the auxiliary dataset is not available
  for the requested time period. If set to `False`, the processor will continue without the auxiliary dataset and
  all output is set to `NaN`. An error will be logged in this case.

Data Variables
^^^^^^^^^^^^^^

The following variables are returned from this dataset:

- ``snow_depth``: Snow depth on sea ice in meters
- ``snow_depth_uncertainty``: Uncertainty estimate for snow depth in meters
- ``snow_density``: Snow density on sea ice :math:`\,^{kg}\!/\!_{m^3}`
- ``snow_density_uncertainty``: Uncertainty estimate for snow density in :math:`\,^{kg}\!/\!_{m^3}`



Data Source and Storage
^^^^^^^^^^^^^^^^^^^^^^^

Data files are available on the ftp of the Alfred Wegener Institute and need to be downloaded into the
directory defined in the ``local_machine_def.yaml`` file under the key ``auxdata_repository.snow.climatology_w99amsr2``.

.. code-block::

    ftp://ftp.awi.de/sea_ice/auxiliary/snow_on_sea_ice/w99_amsr2_merged

(FTP settings: `anonymous login, encryption disabled`)