#! /usr/bin/env python
#
# script to generate Makefile.depend
#
# based on f90dep.pl  written by Pawel Salek
#          mkdep90.py written by Jonas Juselius
#
# radovan bast 2008-05-07
#
# it assumes all included *.h   under ../include
#        and all used     *.mod under ../modules

import sys, string, re, os

use     = re.compile('\s*use\s+(\w+)', re.I)
module  = re.compile('\s*module(?!\s+procedure)\s+(\w+)', re.I)
include = re.compile('include/([\w\.]+)', re.I)

class dep_file_c:

      def __init__(self, name):
          self.file_name   = name
          self.object_name = ''
          ri               = string.rindex(name, '.')
          self.object_name = name[:ri] + '.o'
          self.module_name = ''

          self.include = {}
          self.use     = {}

      def add_include(self, u):
          self.include[u] = ''

      def add_use(self, u):
          self.use[u] = ''

#-------------------------------------------------------------------------------

cpp_enabled = 0

if os.system('cpp < /dev/null > /dev/null') == 0:
   cpp_prefix  = 'cpp '
   cpp_enabled = 1

if os.path.exists('/lib/cpp'):
   cpp_prefix  = '/lib/cpp '
   cpp_enabled = 1

if not cpp_enabled:
   print 'warning: C preprocessing not enabled'

cpp_l  = []
file_l = []

for i in range(len(sys.argv)):
    if i != 0:
       if '-D' in sys.argv[i]:
          cpp_l.append(sys.argv[i])
       else:
          file_l.append(sys.argv[i])

cpp_command = ''
if cpp_enabled:
   cpp_command += ' | ' + cpp_prefix
   for flag in cpp_l:
       cpp_command += flag + ' '
   cpp_command += ' -I../include'

#-------------------------------------------------------------------------------

dep_file_l = []

for file in file_l:

    dep_file_l.append(dep_file_c(file))

    s = os.popen('cat ' + file + cpp_command).readlines()

    for line in s:

        m = module.match(line)
        if m is not None:
           dep_file_l[-1].module_name = string.lower(m.group(1)) + '.mod'

        m = include.search(line)
        if m is not None:
           dep_file_l[-1].add_include(string.lower(m.group(1)))

        m = use.match(line)
        if m is not None:
           dep_file_l[-1].add_use(string.lower(m.group(1)) + '.mod')

#-------------------------------------------------------------------------------

include_path = os.path.abspath('..') + '/include/'
use_path     = os.path.abspath('..') + '/modules/'

for dep_file in dep_file_l:

    for key in dep_file.include.keys():
        if os.path.exists(include_path + key):
           dep_file.include[key] = include_path + key
        else:
           print >> sys.stdout, 'WARNING: missing dependency for file ', include_path + key, 'in', dep_file.file_name

    for key in dep_file.use.keys():
        if os.path.exists(use_path + key):
           dep_file.use[key] = use_path + key
        else:
           print >> sys.stdout, 'WARNING: missing dependency for module ', use_path + key, 'in', dep_file.file_name

#-------------------------------------------------------------------------------

# first round: module names

for dep_file in dep_file_l:
    print '%s: %s' % (dep_file.module_name, dep_file.object_name)

# second round: includes and uses

for dep_file in dep_file_l:
    for key in dep_file.include.keys():
        print '%s: %s' % (dep_file.object_name, dep_file.include[key])
    for key in dep_file.use.keys():
        print '%s: %s' % (dep_file.object_name, dep_file.use[key])
