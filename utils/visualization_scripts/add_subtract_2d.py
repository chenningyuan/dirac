#!/usr/bin/env python
#
# script to add or subtract 2D plots generated with DIRAC
#
# syntax: ./add_subtract_2d.py --file1 "file1.2d" --factor "-1.0" --file2 "file2.2d"
#         will send file1.2d minus file2.2d to stdout
#
# (radovan bast 2008-05-22)

import string
from optparse import OptionParser

#-------------------------------------------------------------------------------

parser = OptionParser()

parser.add_option('--file1',  dest='file1')
parser.add_option('--file2',  dest='file2')
parser.add_option('--factor', dest='factor')

(options, args) = parser.parse_args()

file1  = str(options.file1)
file2  = str(options.file2)
factor = float(str(options.factor))

#-------------------------------------------------------------------------------

f = open(file1, 'r')
file1_s = f.readlines()
f.close()

f = open(file2, 'r')
file2_s = f.readlines()
f.close()

#-------------------------------------------------------------------------------

for i in range(len(file1_s)):

    if (file1_s[i] == '') or ('plane' in file1_s[i]):

       print file1_s[i]

    else:

       l = []

       l.append(float(file1_s[i].split()[0]))
       l.append(float(file1_s[i].split()[1]))
       l.append(float(file1_s[i].split()[2]))
       l.append(float(file1_s[i].split()[3]) + factor*float(file2_s[i].split()[3]))
       l.append(float(file1_s[i].split()[4]) + factor*float(file2_s[i].split()[4]))
       l.append(float(file1_s[i].split()[5]) + factor*float(file2_s[i].split()[5]))
       l.append(float(file1_s[i].split()[6]) + factor*float(file2_s[i].split()[6]))

       print '%18.10e%18.10e%18.10e%18.10e%18.10e%18.10e%18.10e' % tuple(l)
