# -*- coding: utf-8 -*-

"""
plib A collection of python libraries
"""

from loguru import logger

__all__ = ["functions", "iotools", "cs2_procstage2timeliness"]


def cs2_procstage2timeliness(processing_stage_str: str) -> str:
    """
    A small helper function that translates the various processing stage definition into
    a pysiral conformal timeliness definition
    :param processing_stage_str: A String from the CryoSat-2 product metadata
    :return: the 3-character pysiral timeliness code
    """

    # Range of option found in CryoSat-2 L1 metadata:
    timeliness_dct = {"n": "nrt", "nrt_": "nrt",                              # near-real time
                      "o": "ntc",                                             # non-time critical
                      "test": "tds",                                          # test data set
                      "r": "rep", "l": "rep", "offl": "rep", "lta_": "rep"}   # reprocessed

    # Raise a warning if not found since timeliness tag will end up in filename
    if processing_stage_str.lower() not in timeliness_dct.keys():
        msg = "Timeliness metadata tag ({}) not in range of expected values: [{}]"
        msg = msg.format(processing_stage_str.lower(), ", ".join(timeliness_dct.keys()))
        logger.warning(msg)

    return timeliness_dct.get(processing_stage_str.lower(), "ukn")
