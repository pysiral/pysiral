# -*- coding: utf-8 -*-
"""
Created on Mon Sep 05 18:58:18 2016

@author: Stefan
"""

import multiprocessing
import time

def funSquare(num):
    return num ** 2

if __name__ == '__main__':
    t0 = time.clock()
    pool = multiprocessing.Pool()
    results = pool.map(funSquare, range(10))
    t1 = time.clock()
    print(t1-t0)