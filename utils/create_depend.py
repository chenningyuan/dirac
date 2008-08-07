#! /usr/bin/env python
#
# script to generate Makefile.depend
#
# based on f90dep.pl  written by Pawel Salek
#          mkdep90.py written by Jonas Juselius
#
# Patrick Norman and Radovan Bast (Odense 2008-06-11)

import sys
import string
import re
import os
import fnmatch

use     = re.compile('\s*use\s+(\w+)',                    re.I)
module  = re.compile('\s*module(?!\s+procedure)\s+(\w+)', re.I)
include = re.compile('^#.*include/([\w\.]+)',             re.I)

class file_c:

      def __init__(self, path, file_name):
          self.path      = path
          self.file_name = file_name
          self.module    = {}
          self.include   = {}
          self.use       = {}

      def add_include(self, path, file_name):
          self.include[(path, file_name)] = ''

      def add_module(self, u, file_name):
          self.module[u] = file_name

      def add_use(self, u, file_name):
          self.use[u] = file_name

#-------------------------------------------------------------------------------

# directories that are to be excluded
exclude_directory_l = []
exclude_directory_l.append('hsmrcc_external')
exclude_directory_l.append('second_quant')

#-------------------------------------------------------------------------------

cpp_prefix = 'cpp -w -traditional-cpp -Iinclude '
for i in range(len(sys.argv)):
    if i > 0:
       cpp_prefix += sys.argv[i] + ' '

cpp_prefix = string.replace(cpp_prefix, ',', ' ')
cpp_prefix = string.replace(cpp_prefix, '-WF', ' ')

file_l = []
inc_d  = {}

pwd = os.getcwd()

for path, dir, files in os.walk(pwd):

    parse_this_path = 1
    for exclude_directory in exclude_directory_l:
        if exclude_directory in path:
           parse_this_path = 0

    if parse_this_path:
       for file in files:
           if fnmatch.fnmatch(file, '*.F') or fnmatch.fnmatch(file, '*.F90'):
              file_l.append(file_c(path, file))
           if fnmatch.fnmatch(file, '*.h') or fnmatch.fnmatch(file, '*.inc'):
              inc_d[(path, file)] = path + '/' + file

for file in file_l:
    (dummy, stdout, stderr) = os.popen3(cpp_prefix + file.path + '/' + file.file_name)
    for line in stdout.readlines():

#       *.h files
        m = include.search(line)
        if m:
           path = pwd + '/include'
           file_name = string.lower(m.group(1))
           file.add_include(path, file_name)
           file.include[(path, file_name)] = path + '/' + file_name

#       *.h files that are not under include/ and *.inc files
        if (pwd in line) and (('.h" ' in line) or ('.inc" ' in line)):
           between_quotes = string.split(line, '"')[1]
           file_name = string.split(between_quotes, '/')[-1]
           path = between_quotes[0:-len(file_name)-1]
           file.include[(path, file_name)] = path + '/' + file_name

#       *.mod files that are used
        m = use.match(line)
        if m:
           module_name = string.lower(m.group(1))
           file.add_use(module_name, file.file_name)

#       *.mod files that are provided
        m = module.match(line)
        if m:
           module_name = string.lower(m.group(1))
           file.add_module(module_name, file.file_name)

# write out list foo.mod: bar.o (where bar.o contains foo.mod)
for file in file_l:
    without_suffix = file.file_name[:string.rindex(file.file_name, '.')]
    for key in file.module.keys():
        print key + '.mod: ' + without_suffix + '.o'

# file_uses uses a module that is provided by file_provides
# figure out this dependency
for file_uses in file_l:
    for key_uses in file_uses.use.keys():
        for file_provides in file_l:
            for key_provides in file_provides.module.keys():
                if key_uses == key_provides:
                   file_uses.use[key_uses] = file_provides.file_name

for file in file_l:
    without_suffix = file.file_name[:string.rindex(file.file_name, '.')]

#   write out list foo.o: bar.mod (where foo.o uses bar.mod)
    for key in file.use.keys():
#       if a file uses a module that is provided by exactly this file
#       then don't write this dependency
        if file.file_name != file.use[key]:
           print without_suffix + '.o: $(SOURCEDIR)/modules/' + key + '.mod'

#   write out list foo.o: bar.h or bar.inc (where foo.o includes bar.h or bar.inc)
    for key in file.include.keys():
        print without_suffix + '.o: $(SOURCEDIR)' + string.replace(file.include[key], pwd, '')
