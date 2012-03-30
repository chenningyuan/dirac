#!/usr/bin/env python

import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('--build',
                  type='string',
                  action='store',
                  default=None)
parser.add_option('--variable',
                  type='string',
                  action='store',
                  default=None)
parser.add_option('--enable',
                  action='store_true',
                  default=False)
(options, args) = parser.parse_args()

# define line which will be added to the header
add_line = '!#define ' + options.variable
if options.enable: add_line = '#define ' + options.variable

# by default the header is always updated
write_header = True

# check if the header should change
file = options.build + '/unit-tests/unit-test.h'
if os.path.isfile(file):
  f = open(file,'r')
  if f.readline() == add_line: write_header = False
  f.close()

# write new header
if write_header:
  dir = options.build + '/unit-tests'
  if not os.path.isdir(dir): os.makedirs(dir)
  f = open(file,'w')
  f.write('%s' % add_line)
  f.close()
