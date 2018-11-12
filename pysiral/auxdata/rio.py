# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan

Module created for FMI version of pysiral
"""

from pysiral.auxdata import AuxdataBaseClass
import numpy as np


class RIO(AuxdataBaseClass):

    def __init__(self):
        super(RIO, self).__init__()
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

    def _get_along_track_rio(self, l2):

        # Get the data
        self._get_requested_date(l2)
        self._get_data(l2)
        rio = self._get_rio_track(l2)

        # Register the rio code as auxiliary data variables
        self.register_auxvar("data_pc1", "data_pc1", rio[0], None)
        self.register_auxvar("data_pc2", "data_pc2", rio[1], None)
        self.register_auxvar("data_pc3", "data_pc5", rio[2], None)
        self.register_auxvar("data_pc4", "data_pc4", rio[3], None)
        self.register_auxvar("data_pc5", "data_pc5", rio[4], None)
        self.register_auxvar("data_pc6", "data_pc6", rio[5], None)
        self.register_auxvar("data_pc7", "data_pc7", rio[6], None)
        self.register_auxvar("data_1asuper", "data_1asuper", rio[7], None)
        self.register_auxvar("data_1a", "data_1a", rio[8], None)
        self.register_auxvar("data_1b", "data_1b", rio[9], None)
        self.register_auxvar("data_1c", "data_1c", rio[10], None)
        self.register_auxvar("data_no_ice_class", "data_no_ice_class", rio[11], None)

    def _get_requested_date(self, l2):
        """ Use first timestamp as reference, date changes are ignored """
        year = l2.track.timestamp[0].year
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self._requested_date = [year, month, day]
        
    def IceChartToRIO(self, icechart, ice_class="none", summer=False):
        RIO = list()

        for i in np.arange(len(icechart['CT'])):
            RV_A = self.SIGRID3toRV(icechart['SA'][i])
            RV_B = self.SIGRID3toRV(icechart['SB'][i])
            RV_C = self.SIGRID3toRV(icechart['SC'][i])
            conc_A =self. SIGRID3toConc(icechart['CA'][i])
            conc_B = self.SIGRID3toConc(icechart['CB'][i])
            conc_C = self.SIGRID3toConc(icechart['CC'][i])
            conc_T = self.SIGRID3toConc(icechart['CT'][i])
            vec_RV = [RV_A[ice_class]*conc_A, RV_B[ice_class]*conc_B, RV_C[ice_class]*conc_C, (10.0-conc_T)*3]
            if np.isnan(RV_A['PC1']):
                RIO.append(np.nan)
            else:
                RIO.append(np.nansum(vec_RV))

        return RIO

    def SIGRID3toConc(self, int_C):
        #Concentration intervals are handled by taking the median

        if int_C == 0:
            conc = 0.
        elif (int_C == 1) or (int_C == 2):
            conc = 0.5
        elif int_C == 10:
            conc = 1.
        elif int_C == 12:
            conc = 1.5
        elif (int_C == 13) or (int_C == 20):
            conc = 2.
        elif int_C == 23:
            conc = 2.5
        elif (int_C == 24) or (int_C == 30):
            conc = 3.
        elif int_C == 34:
            conc = 3.5
        elif (int_C == 35) or (int_C == 40):
            conc = 4.
        elif int_C ==  45:
            conc = 4.5
        elif (int_C == 46) or (int_C == 50):
            conc = 5.
        elif int_C == 56:
            conc = 5.5
        elif (int_C == 57) or (int_C == 60):
            conc = 6.
        elif int_C == 67:
            conc = 6.5
        elif (int_C == 68) or (int_C == 70):
            conc = 7.
        elif int_C == 78:
            conc = 7.5
        elif (int_C == 79) or (int_C == 80):
            conc = 8.
        elif int_C == 89:
            conc = 8.5
        elif (int_C == 81) or (int_C == 90):
            conc = 9.
        elif int_C == 91:
            conc = 9.5
        elif int_C == 92:
            conc = 10.
        elif int_C == 247: #WH
            conc = 10.
        else:
            conc = np.nan
        return conc

    def SIGRID3toRV(self, ice_type):
        RV = dict()

        if (ice_type == 0) or (ice_type ==00): #ICE-FREE (IMO) ICE-FREE (SG3)
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
        elif (ice_type == 81) or (ice_type ==82): #NEW ICE (IMO) NEW ICE, NILAS, ICE RIND <10cm (SG3)
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
        elif (ice_type == 83) or (ice_type ==85): #GREY WHITE ICE (IMO) YOUNG ICE 10-30cm, GREY-WHITE ICE 15-30cm (SG3)
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
        elif ice_type == 84: #GREY ICE (IMO) GREY ICE 10-15cm (SG3)
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
        elif ice_type == 86: #THICK FY ICE (IMO) FY ICE 30-200cm (SG3)
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
        elif (ice_type == 87) or (ice_type ==89): #THIN FY ICE 2nd STAGE (IMO) THIN FY ICE 30-70cm, THIN FY STAGE 2 50-70cm (SG3)
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
        elif ice_type == 88: #THIN FIRST YEAR 1st STAGE (IMO) THIN FIRST YEAR STAGE 1 30-50cm(SG3)
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
        elif ice_type == 91: #MEDIUM FY ICE (IMO) MEDIUM FY ICE 70-120cm (SG3)
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
        elif ice_type == 93: #THICK FY ICE (IMO) THICK FY ICE >120cm (SG3)
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
        elif (ice_type == 95) or (ice_type == 97) or (ice_type == 98): #HEAVY MY ICE (IMO) OLD ICE, MY ICE, GLACIER ICE (SG3)
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
        elif ice_type == 96: #SECOND YEAR ICE (IMO) SECOND YEAR ICE (SG3)
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
        elif ice_type == 247: #DOES NOT EXIST
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
        if self._requested_date == self._current_date:
            # Data already loaded, nothing to do... right?
            return
            self._data_ct = l2.ic_ct
            self._data_ca = l2.ic_ca
            self._data_cb = l2.ic_cb
            self._data_cc = l2.ic_cc
            self._data_sa = l2.ic_sa
            self._data_sb = l2.ic_sb
            self._data_sc = l2.ic_sc

            self.icechart = dict()
            self.icechart['CT'] = self._data_ct
            self.icechart['CA'] = self._data_ca
            self.icechart['CB'] = self._data_cb
            self.icechart['CC'] = self._data_cc
            self.icechart['SA'] = self._data_sa
            self.icechart['SB'] = self._data_sb
            self.icechart['SC'] = self._data_sc

            self._current_date = self._requested_date

    def _get_rio_track(self, l2):

        # IDEA: Use dictionaries: data[iceclass] = ...
        data_pc1 = list()
        data_pc2 = list()
        data_pc3 = list()
        data_pc4 = list()
        data_pc5 = list()
        data_pc6 = list()
        data_pc7 = list()
        data_1asuper = list()
        data_1a = list()
        data_1b = list()
        data_1c = list()
        data_no_ice_class = list()
        
        iceclasses = ['PC1','PC2','PC3','PC4','PC5','PC6','PC7','1ASuper','1A','1B','1C','NO ICE CLASS']
        icechart = self.icechart
        for str_iceclass in iceclasses:
            RIO = self.IceChartToRIO(icechart, str_iceclass)
            if str_iceclass =='PC1':
                data_pc1 = RIO
            elif str_iceclass =='PC2':
                data_pc2 = RIO
            elif str_iceclass =='PC3':
                data_pc3 = RIO
            elif str_iceclass =='PC4':
                data_pc4 = RIO
            elif str_iceclass =='PC5':
                data_pc5 = RIO
            elif str_iceclass =='PC6':
                data_pc6 = RIO
            elif str_iceclass =='PC7':
                data_pc7 = RIO
            elif str_iceclass =='1ASuper':
                data_1asuper = RIO
            elif str_iceclass =='1A':
                data_1a = RIO
            elif str_iceclass =='1B':
                data_1b = RIO
            elif str_iceclass =='1C':
                data_1c = RIO
            elif str_iceclass =='NO ICE CLASS':
                data_no_ice_class = RIO
            else:
                print "Daaaamn, you shouln't be here!", str_iceclass
		
     
    	return [data_pc1, data_pc2, data_pc3, data_pc4, data_pc5, data_pc6, data_pc7,
                data_1asuper, data_1a, data_1b, data_1c, data_no_ice_class]

