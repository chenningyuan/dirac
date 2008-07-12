#! /usr/bin/env python
#
# script to generate Makefile.depend
#
# based on f90dep.pl  written by Pawel Salek
#          mkdep90.py written by Jonas Juselius
#
# Patrick Norman and Radovan Bast (Odense 2008-06-11)

debug_level = 0

import sys
import string
import re
import os
import fnmatch

use      = re.compile('\s*use\s+(\w+)',                         re.I)
module   = re.compile('\s*module(?!\s+procedure)\s+(\w+)',      re.I)
#include = re.compile('^#include\s+["<]*([a-zA-Z0-9-._]+)[">]', re.I)
include  = re.compile('^#.*include/([\w\.]+)',                  re.I)

class file_c:

      def __init__(self, name):
          self.file_name = name
          self.module    = {}
          self.include   = {}
          self.use       = {}

      def add_module(self, u):
          self.module[u] = ''

      def add_include(self, u):
          self.include[u] = ''

      def add_use(self, u):
          self.use[u] = ''

#-------------------------------------------------------------------------------

cpp_prefix = 'cpp -w -traditional-cpp -Iinclude '
for i in range(len(sys.argv)):
    if i > 0:
       if debug_level > 1:
          cpp_prefix += ''
       else:
          cpp_prefix += sys.argv[i] + ' '
cpp_prefix = string.replace(cpp_prefix, ',', ' ') + ' '

file_l = []
inc_d  = {}

for path, dirs, files in os.walk(os.getcwd()):

    for file in files:
        if fnmatch.fnmatch(file, '*.F'):
           file_l.append(file_c(path + '/' + file))
        if fnmatch.fnmatch(file, '*.F90'):
           file_l.append(file_c(path + '/' + file))
        if fnmatch.fnmatch(file, '*.h'):
           inc_d[file] = path + '/' + file

if debug_level > 0:
   f = open('Makefile.depend_debug', 'w')

pwd = os.getcwd()

for file in file_l:

    (dummy, stdout, stderr) = os.popen3(cpp_prefix + file.file_name)
    for line in stdout.readlines():

        if debug_level > 1:
           f.write('line: ' + line)

        m = include.search(line)
        if m:
           inc = string.lower(m.group(1))
           file.add_include(inc)
           file.include[inc] = inc_d[inc]

           if debug_level > 1:
              f.write('m.group(1): ' + m.group(1) + '\n')
              f.write('inc+inc_d[inc]: ' + inc + ' ' + inc_d[inc] + '\n')

#       this here should pick up *.h files that are not in the include directory
        if (pwd in line) and ('.h" ' in line):
           inc = string.split(string.split(line, '/')[-1], '"')[0]
           file.add_include(inc)
           file.include[inc] = inc_d[inc]

        m = use.match(line)
        if m:
           file.add_use(string.lower(m.group(1)))

        m = module.match(line)
        if m:
           file.add_module(string.lower(m.group(1)))

    if debug_level > 1:
       f.close()
       sys.exit()

if debug_level > 0:
   for file in file_l:
       f.write('file: ' + file.file_name + '\n')
       for key in file.module.keys():
           f.write('file.module: ' + key + '\n')
       for key in file.include.keys():
           f.write('file.include: ' + key + '\n')
       for key in file.use.keys():
           f.write('file.use: ' + key + '\n')
   f.close()

for file in file_l:
    for use in file.use:
        for file_lookup in file_l:
            for key in file_lookup.module.keys():
                if key == use:
                   file.use[use] = file_lookup.file_name

for file in file_l:
    without_suffix = file.file_name[:string.rindex(file.file_name, '.')]
    without_suffix_without_path = string.split(without_suffix, '/')[-1]

    for key in file.module.keys():
        print key + '.mod: ' + without_suffix_without_path + '.o'

for file in file_l:
    without_suffix = file.file_name[:string.rindex(file.file_name, '.')]
    without_suffix_without_path = string.split(without_suffix, '/')[-1]

    for key in file.use.keys():
        # if a file uses a module that is provided by exactly this file
        # then don't write this dependency
        if file.file_name != file.use[key]:
           print without_suffix_without_path + '.o: $(SOURCEDIR)/modules/' + key + '.mod'

    for key in file.include.keys():
        print without_suffix_without_path + '.o: $(SOURCEDIR)' + string.replace(file.include[key], os.getcwd(), '')
