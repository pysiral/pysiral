# -*- coding: utf-8 -*-
"""
Created on Wed Sep 07 19:23:11 2016

@author: Stefan
"""


#import timeit
#
#setup = """
#import numpy as np
#min = 0
#max = 100
#n = 100
#"""
#
#print timeit.timeit('b=np.linspace(min, max, n)', setup=setup)
#
#print timeit.timeit('b=np.arange(n)*float((max-min))/float(n-1)', setup=setup)



import numpy as np
min = 0
max = 100
n = 100


print np.linspace(min, max, n)

b = np.arange(n)*float((max-min))/float(n-1)
print b
print np.arange(n).dtype

# print timeit.timeit('b=bn.nanmax(a)', setup=setup2)