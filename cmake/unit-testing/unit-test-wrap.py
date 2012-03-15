#!/usr/bin/env python

import sys
import subprocess

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--build-dir',
                  type='string',
                  action='store',
                  default=None)
parser.add_option('--name',
                  type='string',
                  action='store',
                  default=None)
(options, args) = parser.parse_args()

command = 'cd %s/test_%s; ./%s.x' % (options.build_dir, options.name, options.name)

p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(out, err) = p.communicate()

if err == '':
    sys.exit(0)
else:
    sys.exit(1)
