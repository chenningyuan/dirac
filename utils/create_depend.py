#! /usr/bin/env python
#
# script to generate Makefile.depend
#
# based on f90dep.pl  written by Pawel Salek
#          mkdep90.py written by Jonas Juselius
#
# Patrick Norman and Radovan Bast (2008-06-11)

import sys, string, re, os, fnmatch

use     = re.compile('\s*use\s+(\w+)', re.I)
module  = re.compile('\s*module(?!\s+procedure)\s+(\w+)', re.I)
#include = re.compile('^#include\s+["<]*([a-zA-Z0-9-._]+)[">]', re.I)
include = re.compile('^#.*include/([\w\.]+)', re.I)

class file_c:

      def __init__(self, name):
          self.file_name   = name

          self.module  = {}
          self.include = {}
          self.use     = {}

      def add_module(self, u):
          self.module[u] = ''

      def add_include(self, u):
          self.include[u] = ''

      def add_use(self, u):
          self.use[u] = ''

#-------------------------------------------------------------------------------


cpp_command = 'cpp -traditional-cpp -Iinclude '
for i in range(len(sys.argv)):
    if i > 0:
       cpp_command += sys.argv[i]
cpp_command += ' '

file_l = []
inc_d = {}


for path, dirs, files in os.walk(os.getcwd()):

    for file in files:
        if fnmatch.fnmatch(file, '*.F'):
           file_l.append(file_c(path + '/' + file))
        if fnmatch.fnmatch(file, '*.F90'):
           file_l.append(file_c(path + '/' + file))
        if fnmatch.fnmatch(file, '*.h'):
           inc_d[file] = path + '/' + file


for file in file_l:

    s = os.popen(cpp_command + file.file_name).readlines()

    for line in s:

        m = module.match(line)
        if m:
           file.add_module(string.lower(m.group(1)))


        m = include.search(line)
        if m:
           inc = string.lower(m.group(1))

           file.add_include(inc)
           file.include[inc] = inc_d[inc]

        m = use.match(line)
        if m:
           file.add_use(string.lower(m.group(1)))



# create dictionary use -> file.mod
for file in file_l:
    for use in file.use:
        for file_lookup in file_l:
            for key in file_lookup.module.keys():
                if  key == use:
                    file.use[use] = file_lookup.file_name
 






# print out
for file in file_l:
    bare_file_name = file.file_name[:string.rindex(file.file_name, '.')]

    for key in file.use.keys():
        print bare_file_name + '.o: ' + os.getcwd() + '/modules/' + key + '.mod'

    for key in file.include.keys():
        print bare_file_name + '.o: ' + file.include[key]
