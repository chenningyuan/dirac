#!/usr/bin/env python
#
# script to extract scalar data from line plot
#
# (trond saue Apr 29 2010)
import string
import math
from optparse import OptionParser
#-----------------------------------------------------------------------
def lines_from_file(file):
    f = open(file, 'r')
    s = f.readlines()
    f.close()
    return s
#-------------------------------------------------------------------------------
usage = "usage: %prog [options] infile arg"
parser = OptionParser()

parser.add_option('--outfile','-o',dest='file1')

(options, args) = parser.parse_args()
#print args
#print options
if len(args) != 1:
        parser.error("incorrect number of arguments")
#---------------------------------------------------------------------
for arg in args:
   f1=open(arg,'r')  
   line=f1.readline()  #read first line
   entry=line.split()
   x0=float(entry[0])
   y0=float(entry[1])
   z0=float(entry[2])
   l=[]
   l.append(float(0))
   l.append(float(entry[3]))
   print '%18.10e%18.10e' % tuple(l)
   while 1:
      line=f1.readline()
      if not line:
         break
      entry=line.split()
      line=f1.readline().split()
      dx=(float(entry[0])-x0)**2
      dy=(float(entry[1])-y0)**2
      dz=(float(entry[2])-z0)**2
      h=math.sqrt(dx+dy+dz)
      l=[]
      l.append(h)
      l.append(float(entry[3]))
      print '%18.10e%18.10e' % tuple(l)
   f1.close()



