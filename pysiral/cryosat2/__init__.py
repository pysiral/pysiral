# -*- coding: utf-8 -*-

"""
plib A collection of python libraries
"""

__all__ = ["functions", "iotools", "l1b_mds_def", "l1bfile", "preproc"]


def cs2_procstage2timeliness(processing_stage_str):
    """
    A small helper function that translates the various processing stage definition into
    a pysiral conformal timeliness definition
    :param processing_stage_str: A String from the CryoSat-2 product metadata
    :return: the 3-character pysiral timeliness code
    """
    timeliness_dct = {"n": "nrt", "o": "ntc", "r": "rep", "l": "rep", "test": "tds", "offl": "rep",
                      "nrt_": "nrt"}
    return timeliness_dct.get(processing_stage_str.lower(), "ukn")