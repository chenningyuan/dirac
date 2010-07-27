#! /usr/bin/env python

#radovan: this is under construction, use at own risk

from optparse import OptionParser
import sys


# dirac module class
# ==================

class dirac_module:
     """bla blub"""
     def __init__(self, flag, cppflag, message):
         self.depends_on  = []
         self.directories = []
         self.files       = []
         self.flag        = flag
         self.cppflag     = cppflag
         self.message     = message
         self.is_included = False
     def add_directory(self, x):
         self.directories.append(x)
     def add_file(self, x):
         self.files.append(x)
     def add_dependency(self, x):
         self.depends_on.append(x)
     def include(self):
         self.is_included = True
     def exclude(self):
         self.is_included = False

def callback_include(option, opt, value, parser):
    dirac_modules[str(option)[2:]].include()


# definition of modules
# =====================

dirac_modules = {}

dirac_modules['visual'] = dirac_module('visual', 'MOD_VISUAL', 'Include visualization module')
dirac_modules['visual'].add_directory('visual')
dirac_modules['visual'].add_file('somefile')
dirac_modules['visual'].add_dependency('matrixop')

dirac_modules['matrixop'] = dirac_module('matrixop', 'MOD_MATRIX_OPERATIONS', 'Include matrix operations module')
dirac_modules['matrixop'].add_directory('matrix_operations')


# define input options
# ======================

parser = OptionParser()

parser.add_option('--tarfile', dest='tarfile')

for module in dirac_modules:
    parser.add_option('--' + module,
                      action="callback",
                      callback=callback_include,
                      help=dirac_modules[module].message)
parser.add_option('-a',
                  '--adm',
                  action='store_true',
                  dest='adm',
                  default=False,
                  help="Include ALL modules (administrator)")


# parse input
# ===========

(options, args) = parser.parse_args()
if options.adm:
   for module in dirac_modules:
       dirac_modules[module].include()


# check dependencies among modules
# ================================

for module in dirac_modules:
    if dirac_modules[module].is_included:
       if len(dirac_modules[module].depends_on) > 0:
          for need_module in dirac_modules[module].depends_on:
              if not dirac_modules[need_module].is_included:
                 error_message = '%s depends on %s but %s is not included' % (module, need_module, need_module)
                 sys.exit(error_message)


# report on screen
# ================

for module in dirac_modules:
    if dirac_modules[module].is_included:
       print dirac_modules[module].message
   #else:
   #   print 'module %s is not included' % module


# prepare dontbuild
# =================


# touch Makefile.dep and Makefile.obj
# ===================================


# prepare Makefile.extobj and Makefile.moddef
# ===========================================



# export dirac to tar.gz file
# ===========================

if options.tarfile != None:
   print 'export', options.tarfile




