# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:57:56 2016

@author: Stefan

Module created for FMI version of pysiral
"""
from pysiral.auxdata import AuxdataBaseClass
#from pysiral.iotools import ReadNC

import scipy.ndimage as ndimage
#from pyproj import Proj
import pyproj
import numpy as np
import os

class ICBaseClass(AuxdataBaseClass):

    def __init__(self):
        super(ICBaseClass, self).__init__()
        self._msg = ""

    def get_along_track_sic(self, l2):
        sic, msg = self._get_along_track_sic(l2)
        return sic, msg


class IC(ICBaseClass):

    def __init__(self):
        super(IC, self).__init__()
        self._data_ct = None
        self._data_ca = None
        self._data_cb = None
        self._data_cc = None
        self._data_sa = None
        self._data_sb = None
        self._data_sc = None
        self._current_date = [0, 0, 0]
        self._requested_date = [-1, -1, -1]
        self.error.caller_id = self.__class__.__name__

    @property
    def year(self):
        return "%04g" % self._requested_date[0]

    @property
    def month(self):
        return "%02g" % self._requested_date[1]

    @property
    def day(self):
        return "%02g" % self._requested_date[2]

    def _initialize(self):
        pass

    def _get_along_track_sic(self, l2):
        self._msg = ""
        self._get_requested_date(l2)
        self._get_data(l2)
        if not self.error.status:
            ic_sics = self._get_sic_ic_track(l2)
            return ic_sics, self._msg
        else:
            return None, self._msg
        
    def _get_requested_date(self, l2):
        """ Use first timestamp as reference, date changes are ignored """
        year = l2.track.timestamp[0].year
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self._requested_date = [year, month, day]
        
    def _open_tif(self, path):
        """ Open the tif file"""
        import gdal
        G = gdal.Open(str(path))
        #print 'icechart G', G, type(G)
        #print 'PATH TO DATA', str(path)
      
        try:
            data = G.ReadAsArray()
        except AttributeError:
            data = np.nan
        G = None
        return data

    def _get_data(self, l2):
        """ Loads file from local repository only if needed """
        if self._requested_date == self._current_date:
            # Data already loaded, nothing to do
            self._msg = "IC: tif already present"
            return
        
        paths,timedelta = self._get_local_repository_filename(l2)
        print 'IN ICECHART _GET_DATA, PATHS 0', paths[0]
        #print os.path.isfile(paths[0])
        
        # Validation
        if not os.path.isfile(paths[0]):
            self._msg = "IC: File not found: %s " % paths[0]
            self.error.add_error("auxdata_missing_sic", self._msg)
            return
        
        self._data_ct = self._open_tif(paths[0])
        self._data_ca = self._open_tif(paths[1])
        self._data_cb = self._open_tif(paths[2])
        self._data_cc = self._open_tif(paths[3])
        self._data_sa = self._open_tif(paths[4])
        self._data_sb = self._open_tif(paths[5])
        self._data_sc = self._open_tif(paths[6])
        
        # to flag the CT values =-9 (corresponding to no data in the tif)
        flagged = np.where(self._data_ct == -9)
        self._data_ct[flagged] = np.nan
        self._data_ca[flagged] = np.nan
        self._data_cb[flagged] = np.nan
        self._data_cc[flagged] = np.nan
        self._data_sa[flagged] = np.nan
        self._data_sb[flagged] = np.nan
        self._data_sc[flagged] = np.nan
        
        # This step is important for calculation of image coordinates
        #self._data.ice_conc = np.flipud(self._data.ice_conc)
        #self._data_ct = np.flipud(self._data_ct)
        #self._data_ca = np.flipud(self._data_ca)
        #self._data_cb = np.flipud(self._data_cb)
        #self._data_cc = np.flipud(self._data_cc)
        #self._data_sa = np.flipud(self._data_sa)
        #self._data_sb = np.flipud(self._data_sb)
        #self._data_sc = np.flipud(self._data_sc)
        
        self._msg = "IC: Loaded IC file: %s (and corresponding CABC, SABC" % paths[0]
        self._current_date = self._requested_date

    def _get_local_repository_filename(self, l2):
        
        import datetime
        time = datetime.datetime(int(self.year),int(self.month),int(self.day))
        
        path = self._local_repository
        
        other_icevars = ['CA','CB','CC','SA','SB','SC']
        for delta in [0,-1,1,-2,2,-3,3,4,-4,5,-5,6,-6,7,-7]:
            fnames = []
            datestring = (time + datetime.timedelta(delta)).strftime('%Y%m%d')
            year = datestring[:4]
            month = datestring[4:6]
            str_filename = path+'/'+str(year)+'/'+str(month)+'/'+"merged_"+datestring+"_CT.tif"
            
            if os.path.isfile(str_filename):
                timedelta = delta
                fnames.append(str_filename)
                for o in other_icevars:
                    fnames.append(path+'/'+str(year)+'/'+str(month)+'/'+'merged_'+datestring+'_'+o+'.tif')
                return fnames, timedelta
                break
            else:
                fnames = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
                timedelta = np.nan
        return fnames, timedelta

    def XXYYGrids(self):
        ### Returns 2 km EASE2 grid XX and YY values

        import numpy as np

        vec_Y = np.arange(5400000, -5400000, -2000) - 1000
        vec_X = np.arange(-5400000, 5400000, 2000) + 1000
        [XX, YY] = np.meshgrid(vec_X, vec_Y)
        return XX,YY

    def _get_sic_ic_track(self, l2):
        # Convert grid/track coordinates to grid projection coordinates
        kwargs = self._options[l2.hemisphere].projection
        #p = Proj(**kwargs)
        data_ct = list()
        data_ca = list()
        data_cb = list()
        data_cc = list()
        data_sa = list()
        data_sb = list()
        data_sc = list()
        
        EASE_proj = pyproj.Proj('+proj=laea +lat_0=90 +lon_0=0 +ellps=WGS84 +datum=WGS84 +units=m')
        latlon_proj = pyproj.Proj('+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs')
        
        vec_x, vec_y =  pyproj.transform(latlon_proj, EASE_proj, l2.track.longitude, l2.track.latitude)
        
        # Convert track projection coordinates to image coordinates
        XX, YY = self.XXYYGrids()
        xOrigin = XX[0][0]
        yOrigin = YY[0][0]
        pixelWidth = XX[0][1]-XX[0][0]
        pixelHeight = YY[1][0]-YY[0][0]
        
        for coords in zip(vec_x,vec_y):
            x = coords[0]
            y = coords[1]
            xOffset = int((x - xOrigin) / pixelWidth)
            yOffset = int((y - yOrigin) / pixelHeight)
            
            stuf =  self._data_ct[yOffset][xOffset]
            dct = str(self._data_ct[yOffset][xOffset]) # AMANDINE: why this conversion to string?
            dca = self._data_ca[yOffset][xOffset]
            dcb = self._data_cb[yOffset][xOffset]
            dcc = self._data_cc[yOffset][xOffset]
            dsa = self._data_sa[yOffset][xOffset]
            dsb = self._data_sb[yOffset][xOffset]
            dsc = self._data_sc[yOffset][xOffset]
            
            # No value in the tif corresponds to -9
            if (dca == -9) and (dcb == -9) and (dcc == -9):
                #print "ICECHART, all is -9"
                dca = dct
            data_ct.append(dct)
            data_ca.append(dca)
            data_cb.append(dcb)
            data_cc.append(dcc)
            data_sa.append(dsa)
            data_sb.append(dsb)
            data_sc.append(dsc)
            
        return [data_ct,data_ca,data_cb,data_cc,data_sa,data_sb,data_sc]

class ICA(ICBaseClass):

    def __init__(self):
        super(ICA, self).__init__() #MUOKS20170517
        self._data_ct = None
        self._data_ca = None
        self._data_cb = None
        self._data_cc = None
        self._data_sa = None
        self._data_sb = None
        self._data_sc = None
        self._current_date = [0, 0, 0]
        self._requested_date = [-1, -1, -1]
        self.error.caller_id = self.__class__.__name__

    @property
    def year(self):
        return "%04g" % self._requested_date[0]

    @property
    def month(self):
        return "%02g" % self._requested_date[1]

    @property
    def day(self):
        return "%02g" % self._requested_date[2]

    def _initialize(self):
        pass

    def _get_along_track_sic(self, l2):
        self._msg = ""
        self._get_requested_date(l2)
        self._get_data(l2)
        ic_sics = self._get_sic_ic_track(l2)
        return ic_sics, self._msg

    def _get_requested_date(self, l2):
        """ Use first timestamp as reference, date changes are ignored """
        year = l2.track.timestamp[0].year
        month = l2.track.timestamp[0].month
        day = l2.track.timestamp[0].day
        self._requested_date = [year, month, day]

    def _open_tif(self, path):
        """ Open the tif file"""
        import gdal
        G = gdal.Open(str(path))
        # NOTE August 2017 exception added for missing AARI years
        try:
            data = G.ReadAsArray()
        except AttributeError:
            data = None
        G = None
        return data

    def _get_data(self, l2):
        """ Loads file from local repository only if needed """
        if self._requested_date == self._current_date:
            # Data already loaded, nothing to do
            self._msg = "ICA: tif already present"
            return
        paths,timedelta = self._get_local_repository_filename(l2)
        #print 'AARI PATHS 0', paths[0]
        
        # Validation
        if not os.path.isfile(paths[0]):
            self._msg = "ICA: File not found: %s " % paths[0]
            self.error.add_error("auxdata_missing_sic", self._msg)
            return
        
        self._data_ct = self._open_tif(paths[0])
        self._data_ca = self._open_tif(paths[1])
        self._data_cb = self._open_tif(paths[2])
        self._data_cc = self._open_tif(paths[3])
        self._data_sa = self._open_tif(paths[4])
        self._data_sb = self._open_tif(paths[5])
        self._data_sc = self._open_tif(paths[6])
        flagged = np.where(self._data_ct == -5)
        # NOTE August 2017 exception added for missing AARI years
        try:
            self._data_ct[flagged] = np.nan
            self._data_ca[flagged] = np.nan
            self._data_cb[flagged] = np.nan
            self._data_cc[flagged] = np.nan
            self._data_sa[flagged] = np.nan
            self._data_sb[flagged] = np.nan
            self._data_sc[flagged] = np.nan
        except TypeError:
            #print 'No AARI icechart for this day'
            pass
        # This step is important for calculation of image coordinates
        #self._data.ice_conc = np.flipud(self._data.ice_conc)
        #self._data_ct = np.flipud(self._data_ct)
        #self._data_ca = np.flipud(self._data_ca)
        #self._data_cb = np.flipud(self._data_cb)
        #self._data_cc = np.flipud(self._data_cc)
        #self._data_sa = np.flipud(self._data_sa)
        #self._data_sb = np.flipud(self._data_sb)
        #self._data_sc = np.flipud(self._data_sc)
        self._msg = "ICA: Loaded IC file: %s (and corresponding CABC, SABC" % paths[0]
        self._current_date = self._requested_date

    def _get_local_repository_filename(self, l2):
        
        import datetime
        time = datetime.datetime(int(self.year),int(self.month),int(self.day))
        
        path = self._local_repository
                
        other_icevars = ['CA','CB','CC','SA','SB','SC']
        for delta in [0,-1,1,-2,2,-3,3,4,-4,5,-5,6,-6,7,-7]:
            fnames = []
            datestring = (time + datetime.timedelta(delta)).strftime('%Y%m%d')
            year = datestring[:4]
            month = datestring[4:6]
            str_filename = path+'/'+str(year)+'/'+str(month)+'/'+"merged_aari_"+datestring+"_CT.tif"
            if os.path.isfile(str_filename):
                timedelta = delta
                fnames.append(str_filename)
                for o in other_icevars:
                    fnames.append(path+'/'+str(year)+'/'+str(month)+'/'+'merged_aari_'+datestring+'_'+o+'.tif')
                return fnames, timedelta
                break
            else:
                fnames = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
                timedelta = np.nan
        return fnames, timedelta

    def XXYYGrids(self):
        ### Returns 2 km EASE2 grid XX and YY values

        import numpy as np

        vec_Y = np.arange(5400000, -5400000, -2000) - 1000
        vec_X = np.arange(-5400000, 5400000, 2000) + 1000
        [XX, YY] = np.meshgrid(vec_X, vec_Y)
        return XX,YY
    
    def _get_sic_ic_track(self, l2):
        # Convert grid/track coordinates to grid projection coordinates
        kwargs = self._options[l2.hemisphere].projection
        #p = Proj(**kwargs)
        data_ct = list()
        data_ca = list()
        data_cb = list()
        data_cc = list()
        data_sa = list()
        data_sb = list()
        data_sc = list()

        EASE_proj = pyproj.Proj('+proj=laea +lat_0=90 +lon_0=0 +ellps=WGS84 +datum=WGS84 +units=m')
        latlon_proj = pyproj.Proj('+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs')
        
        vec_x, vec_y =  pyproj.transform(latlon_proj, EASE_proj, l2.track.longitude, l2.track.latitude)
        
        # Convert track projection coordinates to image coordinates
        XX, YY = self.XXYYGrids()
        xOrigin = XX[0][0]
        yOrigin = YY[0][0]
        pixelWidth = XX[0][1]-XX[0][0]
        pixelHeight = YY[1][0]-YY[0][0]

        for coords in zip(vec_x,vec_y):
            x = coords[0]
            y = coords[1]
            xOffset = int((x - xOrigin) / pixelWidth)
            yOffset = int((y - yOrigin) / pixelHeight)
            
            try:
                stuf =  self._data_ct[yOffset][xOffset]
                data_ct.append(str(self._data_ct[yOffset][xOffset]))
                data_ca.append((self._data_ca[yOffset][xOffset]))
                data_cb.append((self._data_cb[yOffset][xOffset]))
                data_cc.append((self._data_cc[yOffset][xOffset]))
                data_sa.append((self._data_sa[yOffset][xOffset]))
                data_sb.append((self._data_sb[yOffset][xOffset]))
                data_sc.append((self._data_sc[yOffset][xOffset]))
            except TypeError:
                #print 'No AARI icechart for this day'
                data_ct.append(None)
                data_ca.append(None)
                data_cb.append(None)
                data_cc.append(None)
                data_sa.append(None)
                data_sb.append(None)
                data_sc.append(None)

        return [data_ct,data_ca,data_cb,data_cc,data_sa,data_sb,data_sc]

def get_l2_sic_ic_handler(name):
    pyclass = globals().get(name, None)
    if pyclass is not None:
        return pyclass()
    else:
        return pyclass
