#!/usr/bin/env python
#
# script to extract scalar data from line plot
#
# syntax: ./add_subtract_2d.py --file1 "file1.2d" --factor "-1.0" --file2 "file2.2d"
#         will send file1.2d minus file2.2d to stdout
#
# (trond saue Apr 29 2010)
import string
from optparse import OptionParser

#-------------------------------------------------------------------------------
usage = "usage: %prog [options] arg1"
parser = OptionParser()

parser.add_option('--file1',  dest='file1')

(options, args) = parser.parse_args()
print args
print options

optparse


