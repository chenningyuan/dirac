#! /usr/bin/env python

#radovan: this is under construction, use at own risk

from optparse import OptionParser
import sys


# dirac module class
# ==================

class dirac_module:
     """bla blub"""
     def __init__(self, flag, message):
         self.flag        = flag
         self.message     = message
         self.is_included = False
         self.depends_on  = []
         self.directories = []
         self.files       = []
         self.cppflags    = []
     def add_directory(self, x):
         self.directories.append(x)
     def add_file(self, x):
         self.files.append(x)
     def add_dependency(self, x):
         self.depends_on.append(x)
     def add_cppflag(self, x):
         self.cppflags.append(x)
     def include(self):
         self.is_included = True
     def exclude(self):
         self.is_included = False

def callback_include(option, opt, value, parser):
    dirac_modules[str(option)[2:]].include()


# definition of modules
# =====================

dirac_modules = {}

dirac_modules['visual'] = dirac_module('visual', 'Include visualization module')
dirac_modules['visual'].add_cppflag('MOD_VISUAL')
dirac_modules['visual'].add_directory('visual')
dirac_modules['visual'].add_file('somefile')
dirac_modules['visual'].add_dependency('matrixop')

dirac_modules['matrixop'] = dirac_module('matrixop', 'Include matrix operations module')
dirac_modules['matrixop'].add_cppflag('MOD_MATRIX_OPERATIONS')
dirac_modules['matrixop'].add_directory('matrix_operations')

dirac_modules['eri'] = dirac_module('eri', 'Include ERI integral routines')
dirac_modules['eri'].add_cppflag('MOD_ERI')

dirac_modules['densfit'] = dirac_module('densfit', 'Include density fitting module')
dirac_modules['densfit'].add_cppflag('MOD_DNF')

dirac_modules['largemolecules'] = dirac_module('largemolecules', 'Increases common blocks to handle more atoms and basis sets')
dirac_modules['largemolecules'].add_cppflag('LARGE_MOLECULES')

dirac_modules['fmm'] = dirac_module('fmm', 'Include fast multipole method')
dirac_modules['fmm'].add_cppflag('MOD_FMM')

dirac_modules['esr'] = dirac_module('esr', 'Include ESR property module')
dirac_modules['esr'].add_cppflag('MOD_ESR')

dirac_modules['xpp'] = dirac_module('xpp', 'Include PP excitation energy and STEX modules')
dirac_modules['xpp'].add_cppflag('MOD_XPP')
dirac_modules['xpp'].add_cppflag('MOD_STEX')

dirac_modules['xml'] = dirac_module('xml', 'Include XML input/output module')
dirac_modules['xml'].add_cppflag('MOD_XML')

dirac_modules['krmc'] = dirac_module('krmc', 'Include KR-MCSCF and LUCIAREL modules')
dirac_modules['krmc'].add_cppflag('MOD_KRMC')

dirac_modules['emb'] = dirac_module('emb', 'Include ab initio/DFT-embedding module')
dirac_modules['emb'].add_cppflag('MOD_EMB')

dirac_modules['openrsp'] = dirac_module('openrsp', 'Include OpenRSP')
dirac_modules['openrsp'].add_cppflag('MOD_OPENRSP')

dirac_modules['xcfun'] = dirac_module('xcfun', 'Include XCFun module')
dirac_modules['xcfun'].add_cppflag('MOD_XCFUN')

dirac_modules['cap'] = dirac_module('cap', 'Include CAP module')
dirac_modules['cap'].add_cppflag('MOD_CAP')

#fixme: in pamadm this flag triggers a lot of operations
#       so far this is not done here
#       (can they be avoided?)
dirac_modules['hscc'] = dirac_module('hscc', 'Include interface for the (external) HSMRCC code together with its objects')
dirac_modules['hscc'].add_cppflag('MOD_HSCC')

#fixme: in pamadm this flag triggers a lot of operations
#       so far this is not done here
#       (can they be avoided?)
dirac_modules['hsfs'] = dirac_module('hsfs', 'Include higher external Higher sector FSCC code')
dirac_modules['hsfs'].add_cppflag('MOD_HSFS')

#fixme: in pamadm this flag triggers a lot of operations
#       so far this is not done here
#       (can they be avoided?)
dirac_modules['hsfsdyn'] = dirac_module('hsfsdyn', 'Include higher external Higher sector FSCC code - dynamic allocation version.')
dirac_modules['hsfsdyn'].add_cppflag('HSFS_DYN_ALLOC')

dirac_modules['qdots'] = dirac_module('qdots', 'Include Quantum Dots functionality (artificial atom)')
dirac_modules['qdots'].add_cppflag('MOD_QDOTS')

dirac_modules['ihfscc'] = dirac_module('ihfscc', 'Include Intermediate Hamiltonian functionality to FSCC code')
dirac_modules['ihfscc'].add_cppflag('MOD_IHFSCC')

dirac_modules['amfi'] = dirac_module('amfi', 'Include the AMFI functionality into DIRAC')
dirac_modules['amfi'].add_cppflag('MOD_AMFI')

dirac_modules['london'] = dirac_module('london', 'Include London orbitals functionality for magnetic properties')
dirac_modules['london'].add_cppflag('MOD_LONDON')

dirac_modules['oof'] = dirac_module('oof', 'Include object oriented modularity of Fotran 90/95')
dirac_modules['oof'].add_cppflag('MOD_OOF')

dirac_modules['adc'] = dirac_module('adc', 'Include ADC electron propagator module')
dirac_modules['adc'].add_cppflag('MOD_ADC')

dirac_modules['ecp'] = dirac_module('ecp', 'Include ARGOS ECP module')
dirac_modules['ecp'].add_cppflag('MOD_ECP')

dirac_modules['srdft'] = dirac_module('srdft', 'Include srDFT module')
dirac_modules['srdft'].add_cppflag('MOD_SRDFT')
dirac_modules['srdft'].add_cppflag('MOD_KRMC')

dirac_modules['unreleased'] = dirac_module('unreleased', 'Include unreleased functionality (not covered by other flags)')
dirac_modules['unreleased'].add_cppflag('MOD_UNRELEASED')

dirac_modules['log'] = dirac_module('log', 'Will log interesting matrices in matlab format (very verbose!)')
dirac_modules['log'].add_cppflag('MOD_LOG')


# define input options
# ======================

parser = OptionParser()

parser.add_option('--tarfile', 
                  dest='tarfile',
                  help='Export or release DIRAC to TARFILE')

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
