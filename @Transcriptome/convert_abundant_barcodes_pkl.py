#!/usr/bin/python

import pickle
import sys

data = pickle.load( open(sys.argv[1], "rb") )
f=open(sys.argv[2], 'w+')

for k,v in data.iteritems():
    f.write("%s,%s,%d\n" % (k, v[0], v[1]))

f.close()
