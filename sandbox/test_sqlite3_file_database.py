# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 23:37:56 2016

@author: Stefan

Testbed for creating a relational database that contains it is searchable
for periods, regions, etc. Thus not all files need to be parsed by the
altimetry processor, e.g. skip all files for the Antarctic for Artic run etc.

The database is supposed only to be implemented for the l1bdata nc files, since
most of the files have to be parsed and validation to get the information for
the table

Filenaming: "pysiral_l1bdata_inventory.db"

Main rows that need to be in the database

Table: SOURCE_FILES
filename         serves also as unique identifier (str)
year	           Year of start time (int)
month            Month of start time (int)
week             Week number of start time (int)
day              Day of start time (int)
orbit            orbit number (int)
cycle            cycle number (int)
day_overlap      flag that indicating data coverage on the next day (bool)
region_id        id for region (see table REGION) (int)
radar_mode_id    id for radar mode (see table RADAR_MODE) (int)
mission_id       id for mission (see table MISSION) (int)
version          version id of altimeter source data (str)
path             relative link to file (str)



Table: REGION
region_id (int), full_name (str)
0: arctic
1: antarctic
2: global

Table: RADAR_MODE
radar_mode_id (int), full_name (str)
0: LRM
1: SAR
2: SIN

Table: MISSION_ID
mission_id (int), full_name (str)
0: ers1
1: ers2
2: envisat
3: icesat
4: cryosat2
5: altika
6: sentinel3a
7: sentinel3b
8: icesat2

Other information:

Table: PATH_INCLUDED  (list of directories that have been search)
directory (str)
"""

import os

start_dir = r"D:\awi\altim\data\altimetry\ers1\sgdr"
pattern   = "*.log"

for dirpath, dirnames, filenames in os.walk(start_dir):
    print dirpath
    for filename in filenames:
        print " "+filename
