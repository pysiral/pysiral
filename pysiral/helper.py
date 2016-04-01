# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 18:04:43 2015

@author: Stefan
"""

from dateutil import parser as dtparser
import numpy as np


def parse_datetime_str(dtstr):
    """ Converts a time string to a datetime object using dateutils """
    return dtparser.parse(dtstr)


def get_first_array_index(array, value):
    """ Get the index in array of the first occurance of ``value`` """
    try:
        index = list(array).index(value)
    except:
        index = None
    return index


def get_last_array_index(array, value):
    """ Get the index in array of the last occurance of ``value`` """
    listarray = list(array)
    try:
        index = (len(listarray) - 1) - listarray[::-1].index(value)
    except:
        index = None
    return index


def rle(inarray):
    """
    run length encoding. Partial credit to R rle function.
    Multi datatype arrays catered for including non Numpy
    returns: tuple (runlengths, startpositions, values)

    from: http://stackoverflow.com/questions/1066758/find-length-of-sequences-
                 of-identical-values-in-a-numpy-array
    """
    ia = np.array(inarray)                   # force numpy
    n = len(ia)
    if n == 0:
        return (None, None, None)
    else:
        y = np.array(ia[1:] != ia[:-1])      # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)    # must include last element posi
        z = np.diff(np.append(-1, i))        # run lengths
        p = np.cumsum(np.append(0, z))[:-1]  # positions
        return(z, p, ia[i])
