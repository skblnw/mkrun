#!/usr/bin/env python
import sys
import getopt
import numpy
import os.path
import math

#### Read input from command line ####
input_file_name=''

try:
        myopts, args = getopt.getopt(sys.argv[1:],"i:")
except getopt.GetoptError, err:
        print "Error: ", str(err)
        print "Usage: ./generate_hist.py -i dist.dat"
        sys.exit(2)

for o, a in myopts:
        if o == "-i":
                input_file_name = a
        else:
                print 'Error'
                print "Usage: ./generate_hist.py -i dist.dat"
                sys.exit(0)

#### Check options are set ####
if input_file_name == '':
        print 'Error: cannot find dist.dat'
        sys.exit(2)
elif os.path.isfile(input_file_name) == False:
	print 'Error: cannot find',input_file_name
	sys.exit(2)


# Fill numpy array:
length=0
with open(input_file_name, 'r') as f:
	for line in f:
		length+=1

# In-built here: distance data is second column!
i=0
dist=numpy.zeros((length,1))
with open(input_file_name, 'r') as f:
	for line in f:
		dist[i]=float(line.split()[1])
		i+=1

# Make histogram of data
hist,bin=numpy.histogram(dist,bins=100,range=(30,62))

# Print out
for b, h in zip(bin[1:],hist):
	print b, h

