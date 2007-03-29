#!/usr/bin/env python

import sys, string, os, re
from common import *

#-------------------------------------------------------------------------------

open_page = '''
<html>
<title> error report for dirac test runs </title>
<body text="#ffffff" bgcolor="#000000">
<pre>
'''

summary_file = 'summary.htm'

os.system('rm -f ' + summary_file)
to_file(open_page, summary_file)

for atom in d_atom:
    for sym in d_symtxt:
        if sym != 'c1':
           for field in d_operator:
               for orientation in d_orientation:
                   for quantity in d_quantity[field]:

                       file1 = 'plot_' + field + d_orientation[orientation] + '_' + quantity + '_' + atom + '_' + sym
                       file2 = 'plot_' + field + d_orientation[orientation] + '_' + quantity + '_' + atom + '_' + 'c1'

                       command = '../crep ' + file1.lower() + '.PLOT filter > ' + file1.lower() + '_deleteme'
                       os.system(command)
                       command = '../crep ' + file2.lower() + '.PLOT filter > ' + file2.lower() + '_deleteme'
                       os.system(command)
                       command = '../hdrep ' + file1.lower() + '_deleteme ' + file2.lower() + '_deleteme >> ' + summary_file
                       os.system(command)

os.system("echo '</body></html>' >> " + summary_file)
os.system('rm -f *deleteme*')
