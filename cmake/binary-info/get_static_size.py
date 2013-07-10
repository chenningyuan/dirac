#!/usr/bin/env python
#
# get static size of binary
# written by: Radovan Bast <bast@kth.se>
#             Miroslav Ilias <miroslav.ilias@umb.sk>

import subprocess
import sys
import os

nr_top_consumers = 20

if len(sys.argv) == 1:
    print('usage: %s binary' % sys.argv[0])
    sys.exit(1)
else:
    binary = sys.argv[1]

if not os.path.isfile(binary):
    print('binary %s does not exist' % binary)
    sys.exit(1)

command = "size %s 2> /dev/null | grep %s" % (binary, binary)
process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
output = process.stdout.read().strip()

print('total static size of %s:' % binary)

total_size = float(output.split()[2])/(1024*1024)
print('%10.2f MB' % total_size)

command = "readelf -s %s 2> /dev/null || echo 'cannot be determined, readelf not present'" % binary
process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

l = []
for line in process.stdout.readlines():
    if (line.find('DEFAULT') != -1) and (line.find('FUNC') != -1 or line.find('OBJECT') != -1) and (line.find('UND') == -1):
        size = int(line.split()[2], 16) # get the size in hex number
        name = line.split()[7]          # extract the name
        type = line.split()[-5]         # extract type, necessary to identify variable
        l.append([size, name, type])
l.sort(reverse=True)

print('\ntop %i consumers:' % nr_top_consumers)
for i in range(nr_top_consumers):
    ii = 2*i + 1
    print('%10.2f MB %20s  (%8s)' % (float(l[ii][0])/(1024*1024), l[ii][1], l[ii][2]))
