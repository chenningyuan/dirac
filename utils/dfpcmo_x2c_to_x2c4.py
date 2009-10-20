#!/usr/bin/env python

# radovan.bast@uit.no (2009--10--20):

# script to convert X2C DFPCMO to X2C4 DFPCMO
# this is useful if you have a converged X2C DFPCMO
# and need this as initial guess for a 4C (DC) calculation

# you need to do: 1. convert X2C DFPCMO to X2C4 DFPCMO
#                 2. run X2C4 from X2C4 DFPCMO, save DFCOEF
#                 3. restart 4C calculation from DFCOEF

# example: file is 2c_dfpcmo and there are large+small = 95 basis functions
#          ./dfpcmo_x2c_to_x2c4.py 2c_dfpcmo 95 > 4c_dfpcmo

import sys
import string

def write(output, f):
    output += '%22.16f' % f
    if len(output.split()) == 6:
       print output
       output = ''
    return output

file_name = sys.argv[1]
ao_total  = int(sys.argv[2])

f = open(file_name, 'r')
s = f.read()
f.close()

s = string.split(s, '\n')

print s[0]

mo       = int(s[1].split()[2])
ao_large = int(s[1].split()[3])

ao_small = ao_total - ao_large

nr_lines = (mo*ao_large*4)//6
if (mo*ao_large*4)%6 > 0:
   nr_lines += 1

print '   1         0%10i%10i' % (mo, ao_large + ao_small)

print s[2]

array = []
for i in range(3, nr_lines + 3):
    for w in s[i].split()[0:]:
        array.append(float(w))

i = 0
output = ''
for iz in [1, 2, 3, 4]:
    for k in range(mo):
        for l in range(ao_large):
            output = write(output, array[i])
            i += 1
        for l in range(ao_small):
            output = write(output, 0.0)

if output != '':
   print output

for i in range(nr_lines + 3, len(s) - 1):
    print s[i]
