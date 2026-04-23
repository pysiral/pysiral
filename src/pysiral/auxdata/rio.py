# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan

Module created for FMI version of pysiral
"""

import numpy as np

from pysiral.auxdata import AuxdataBaseClass


class RIO(AuxdataBaseClass):

    def __init__(self, *args, **kwargs):
        super(RIO, self).__init__(*args, **kwargs)
        self._rio_pc1 = None
        self._rio_pc2 = None
        self._rio_pc3 = None
        self._rio_pc4 = None
        self._rio_pc5 = None
        self._rio_pc6 = None
        self._rio_pc7 = None
        self._rio_1asuper = None
        self._rio_1a = None
        self._rio_1b = None
        self._rio_1c = None
        self._rio_no_ice_class = None

    def get_l2_track_vars(self, l2):

        # Get the data
        self._get_requested_date(l2)
        self._get_data(l2)
        rio = self._get_rio_track(l2)

        # Register the rio code as auxiliary data variables
        self.register_auxvar("rio_pc1", "rio_pc1", rio[0], None)
        self.register_auxvar("rio_pc2", "rio_pc2", rio[1], None)
        self.register_auxvar("rio_pc3", "rio_pc3", rio[2], None)
        self.register_auxvar("rio_pc4", "rio_pc4", rio[3], None)
        self.register_auxvar("rio_pc5", "rio_pc5", rio[4], None)
        self.register_auxvar("rio_pc6", "rio_pc6", rio[5], None)
        self.register_auxvar("rio_pc7", "rio_pc7", rio[6], None)
        self.register_auxvar("rio_1asuper", "rio_1asuper", rio[7], None)
        self.register_auxvar("rio_1a", "rio_1a", rio[8], None)
        self.register_auxvar("rio_1b", "rio_1b", rio[9], None)
        self.register_auxvar("rio_1c", "rio_1c", rio[10], None)
        self.register_auxvar("rio_no_ice_class", "rio_no_ice_class", rio[11], None)

    def _get_requested_date(self, l2):
        """ Use first timestamp as reference, date changes are ignored """
        year = l2.track.timestamp[0].year
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self._requested_date = [year, month, day]

    def ice_chart_to_rio(self, icechart, ice_class="none", summer=False):
        rio_list = []
        for i in np.arange(len(icechart['CT'])):
            RV_A = self.sigrid3_to_rv(icechart['SA'][i])
            RV_B = self.sigrid3_to_rv(icechart['SB'][i])
            RV_C = self.sigrid3_to_rv(icechart['SC'][i])
            conc_A = self.sigrid3_to_conc(icechart['CA'][i])
            conc_B = self.sigrid3_to_conc(icechart['CB'][i])
            conc_C = self.sigrid3_to_conc(icechart['CC'][i])
            conc_T = self.sigrid3_to_conc(icechart['CT'][i])
            vec_RV = [RV_A[ice_class]*conc_A, RV_B[ice_class]*conc_B, RV_C[ice_class]*conc_C, (10.0-conc_T)*3]
            if np.isnan(RV_A['PC1']):
                rio_list.append(np.nan)
            else:
                rio_list.append(np.nansum(vec_RV))
        return rio_list

    @staticmethod
    def sigrid3_to_conc(int_c):

        # Concentration intervals are handled by taking the median
        if int_c == 0:
            return 0.
        elif int_c in [1, 2]:
            return 0.5
        elif int_c == 10:
            return 1.
        elif int_c == 12:
            return 1.5
        elif int_c in [13, 20]:
            return 2.
        elif int_c == 23:
            return 2.5
        elif int_c in [24, 30]:
            return 3.
        elif int_c == 34:
            return 3.5
        elif int_c in [35, 40]:
            return 4.
        elif int_c == 45:
            return 4.5
        elif int_c in [46, 50]:
            return 5.
        elif int_c == 56:
            return 5.5
        elif int_c in [57, 60]:
            return 6.
        elif int_c == 67:
            return 6.5
        elif int_c in [68, 70]:
            return 7.
        elif int_c == 78:
            return 7.5
        elif int_c in [79, 80]:
            return 8.
        elif int_c == 89:
            return 8.5
        elif int_c in [81, 90]:
            return 9.
        elif int_c == 91:
            return 9.5
        elif int_c == 92:
            return 10.
        elif int_c == 247:  # WH
            return 10.
        else:
            return np.nan

    @staticmethod
    def sigrid3_to_rv(ice_type):

        RV = {}

        # ICE-FREE (IMO) ICE-FREE (SG3)
        if ice_type == 0:
            RV['PC1'] = 3
            RV['PC2'] = 3
            RV['PC3'] = 3
            RV['PC4'] = 3
            RV['PC5'] = 3
            RV['PC6'] = 3
            RV['PC7'] = 3
            RV['1ASuper'] = 3
            RV['1A'] = 3
            RV['1B'] = 3
            RV['1C'] = 3
            RV['NO ICE CLASS'] = 3

        # NEW ICE (IMO) NEW ICE, NILAS, ICE RIND <10cm (SG3)
        elif ice_type in [81, 82]:
            RV['PC1'] = 3
            RV['PC2'] = 3
            RV['PC3'] = 3
            RV['PC4'] = 3
            RV['PC5'] = 3
            RV['PC6'] = 2
            RV['PC7'] = 2
            RV['1ASuper'] = 2
            RV['1A'] = 2
            RV['1B'] = 2
            RV['1C'] = 2
            RV['NO ICE CLASS'] = 1

        # GREY WHITE ICE (IMO) YOUNG ICE 10-30cm, GREY-WHITE ICE 15-30cm (SG3)
        elif ice_type in [83, 85]:
            RV['PC1'] = 3
            RV['PC2'] = 3
            RV['PC3'] = 3
            RV['PC4'] = 3
            RV['PC5'] = 3
            RV['PC6'] = 2
            RV['PC7'] = 2
            RV['1ASuper'] = 2
            RV['1A'] = 2
            RV['1B'] = 1
            RV['1C'] = 0
            RV['NO ICE CLASS'] = -1

        # GREY ICE (IMO) GREY ICE 10-15cm (SG3)
        elif ice_type == 84:
            RV['PC1'] = 3
            RV['PC2'] = 3
            RV['PC3'] = 3
            RV['PC4'] = 3
            RV['PC5'] = 3
            RV['PC6'] = 2
            RV['PC7'] = 2
            RV['1ASuper'] = 2
            RV['1A'] = 2
            RV['1B'] = 2
            RV['1C'] = 1
            RV['NO ICE CLASS'] = 0

        # THICK FY ICE (IMO) FY ICE 30-200cm (SG3)
        elif ice_type == 86:
            RV['PC1'] = 2
            RV['PC2'] = 2
            RV['PC3'] = 2
            RV['PC4'] = 1
            RV['PC5'] = 0
            RV['PC6'] = -1
            RV['PC7'] = -2
            RV['1ASuper'] = -2
            RV['1A'] = -3
            RV['1B'] = -4
            RV['1C'] = -5
            RV['NO ICE CLASS'] = -6

        # THIN FY ICE 2nd STAGE (IMO) THIN FY ICE 30-70cm, THIN FY STAGE 2 50-70cm (SG3)
        elif ice_type in [87, 89]:
            RV['PC1'] = 2
            RV['PC2'] = 2
            RV['PC3'] = 2
            RV['PC4'] = 2
            RV['PC5'] = 2
            RV['PC6'] = 1
            RV['PC7'] = 1
            RV['1ASuper'] = 1
            RV['1A'] = 0
            RV['1B'] = -1
            RV['1C'] = -2
            RV['NO ICE CLASS'] = -3

        # THIN FIRST YEAR 1st STAGE (IMO) THIN FIRST YEAR STAGE 1 30-50cm(SG3)
        elif ice_type == 88:
            RV['PC1'] = 2
            RV['PC2'] = 2
            RV['PC3'] = 2
            RV['PC4'] = 2
            RV['PC5'] = 2
            RV['PC6'] = 2
            RV['PC7'] = 1
            RV['1ASuper'] = 2
            RV['1A'] = 1
            RV['1B'] = 0
            RV['1C'] = -1
            RV['NO ICE CLASS'] = -2

        # MEDIUM FY ICE (IMO) MEDIUM FY ICE 70-120cm (SG3)
        elif ice_type == 91:
            RV['PC1'] = 2
            RV['PC2'] = 2
            RV['PC3'] = 2
            RV['PC4'] = 2
            RV['PC5'] = 1
            RV['PC6'] = 0
            RV['PC7'] = -1
            RV['1ASuper'] = -1
            RV['1A'] = -2
            RV['1B'] = -3
            RV['1C'] = -4
            RV['NO ICE CLASS'] = -5

        # THICK FY ICE (IMO) THICK FY ICE >120cm (SG3)
        elif ice_type == 93:
            RV['PC1'] = 2
            RV['PC2'] = 2
            RV['PC3'] = 2
            RV['PC4'] = 1
            RV['PC5'] = 0
            RV['PC6'] = -1
            RV['PC7'] = -2
            RV['1ASuper'] = -2
            RV['1A'] = -3
            RV['1B'] = -4
            RV['1C'] = -5
            RV['NO ICE CLASS'] = -6

        # HEAVY MY ICE (IMO) OLD ICE, MY ICE, GLACIER ICE (SG3)
        elif ice_type in [95, 97, 98]:
            RV['PC1'] = 1
            RV['PC2'] = 0
            RV['PC3'] = -1
            RV['PC4'] = -2
            RV['PC5'] = -2
            RV['PC6'] = -3
            RV['PC7'] = -3
            RV['1ASuper'] = -4
            RV['1A'] = -5
            RV['1B'] = -6
            RV['1C'] = -8
            RV['NO ICE CLASS'] = -8

        # SECOND YEAR ICE (IMO) SECOND YEAR ICE (SG3)
        elif ice_type == 96:
            RV['PC1'] = 2
            RV['PC2'] = 1
            RV['PC3'] = 1
            RV['PC4'] = 0
            RV['PC5'] = -1
            RV['PC6'] = -2
            RV['PC7'] = -3
            RV['1ASuper'] = -3
            RV['1A'] = -4
            RV['1B'] = -5
            RV['1C'] = -6
            RV['NO ICE CLASS'] = -7

        # DOES NOT EXIST
        elif ice_type == 247:
            RV['PC1'] = 0
            RV['PC2'] = 0
            RV['PC3'] = 0
            RV['PC4'] = 0
            RV['PC5'] = 0
            RV['PC6'] = 0
            RV['PC7'] = 0
            RV['1ASuper'] = 0
            RV['1A'] = 0
            RV['1B'] = 0
            RV['1C'] = 0
            RV['NO ICE CLASS'] = 0

        # Backup
        else:
            RV['PC1'] = np.nan
            RV['PC2'] = np.nan
            RV['PC3'] = np.nan
            RV['PC4'] = np.nan
            RV['PC5'] = np.nan
            RV['PC6'] = np.nan
            RV['PC7'] = np.nan
            RV['1ASuper'] = np.nan
            RV['1A'] = np.nan
            RV['1B'] = np.nan
            RV['1C'] = np.nan
            RV['NO ICE CLASS'] = np.nan

        return RV

    def _get_data(self, l2):
        """ Hopefully gets data assigned already for l2 """

        self._data_ct = l2.ic_ct
        self._data_ca = l2.ic_ca
        self._data_cb = l2.ic_cb
        self._data_cc = l2.ic_cc
        self._data_sa = l2.ic_sa
        self._data_sb = l2.ic_sb
        self._data_sc = l2.ic_sc

        self.icechart = {
            'CT': self._data_ct,
            'CA': self._data_ca,
            'CB': self._data_cb,
            'CC': self._data_cc,
            'SA': self._data_sa,
            'SB': self._data_sb,
            'SC': self._data_sc
        }

        # self._current_date = self._requested_date

    def _get_rio_track(self, l2):

        # IDEA: Use dictionaries: data[iceclass] = ...

        data_pc1 = []
        data_pc2 = []
        data_pc3 = []
        data_pc4 = []
        data_pc5 = []
        data_pc6 = []
        data_pc7 = []
        data_1asuper = []
        data_1a = []
        data_1b = []
        data_1c = []
        data_no_ice_class = []

        iceclasses = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', '1ASuper', '1A', '1B', '1C', 'NO ICE CLASS']
        icechart = self.icechart
        for str_iceclass in iceclasses:
            rio_list = self.ice_chart_to_rio(icechart, str_iceclass)
            if str_iceclass == 'PC1':
                data_pc1 = rio_list
            elif str_iceclass == 'PC2':
                data_pc2 = rio_list
            elif str_iceclass == 'PC3':
                data_pc3 = rio_list
            elif str_iceclass == 'PC4':
                data_pc4 = rio_list
            elif str_iceclass == 'PC5':
                data_pc5 = rio_list
            elif str_iceclass == 'PC6':
                data_pc6 = rio_list
            elif str_iceclass == 'PC7':
                data_pc7 = rio_list
            elif str_iceclass == '1ASuper':
                data_1asuper = rio_list
            elif str_iceclass == '1A':
                data_1a = rio_list
            elif str_iceclass == '1B':
                data_1b = rio_list
            elif str_iceclass == '1C':
                data_1c = rio_list
            elif str_iceclass == 'NO ICE CLASS':
                data_no_ice_class = rio_list
            else:
                print("Daaaamn, you shouln't be here!", str_iceclass)

        return [data_pc1, data_pc2, data_pc3, data_pc4, data_pc5, data_pc6, data_pc7,
                data_1asuper, data_1a, data_1b, data_1c, data_no_ice_class]
