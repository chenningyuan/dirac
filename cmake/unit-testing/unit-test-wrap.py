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

test_dir = '%s/unit-tests/%s' % (options.build_dir, options.name)

command = 'cd %s; ./%s.x' % (test_dir, options.name)

p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(out, err) = p.communicate()

f = open('%s/%s.log' % (test_dir, options.name), 'w')
f.write(out)
f.close()

f = open('%s/%s.err' % (test_dir, options.name), 'w')
f.write(err)
f.close()

if err == '':
    sys.exit(0)
else:
    sys.exit(1)
