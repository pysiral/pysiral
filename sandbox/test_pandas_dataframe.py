# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 20:43:36 2015

@author: Stefan
"""

from pandas import DataFrame, Series

d = [{'a': 1, 'b': 2}, {'a': 5, 'b': 10}]
df = DataFrame(d)

print df["a"]
