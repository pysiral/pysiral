Snow on Sea Ice (snow)
======================


Warren99 & AMSR2 Monthly Climatology (clim_w99amsr2)
----------------------------------------------------

The monthly Warren99 and AMSR2 snow climatology provides monthly snow depth and density values
for the Arctic winter months from October through April. 



Configuration Parameters
^^^^^^^^^^^^^^^^^^^^^^^^


.. code-block:: yaml

    - snow:
        name: clim_w99amsr2
        options:
            daily_scaling: True
            fyi_correction_factor: 0.5
            exception_on_error: False


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