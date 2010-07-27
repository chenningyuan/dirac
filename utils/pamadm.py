#! /usr/bin/env python

#radovan: this is under construction, use at own risk

from optparse import OptionParser
import sys
import os


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
         self.part_of_adm = True
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
     def take_out_of_adm(self):
         self.part_of_adm = False

def callback_include(option, opt, value, parser):
    dirac_modules[str(option)[2:]].include()


# list vital directories and files
# (these are always exported/released)
# ====================================

vital_directories = []
vital_directories.append('dirac')
vital_directories.append('prp')
vital_directories.append('abacus')
vital_directories.append('gp')
vital_directories.append('pdpack')
vital_directories.append('main')
vital_directories.append('relccsd')
vital_directories.append('moltra')
vital_directories.append('dirrci')
vital_directories.append('lucita')
vital_directories.append('utils')
vital_directories.append('dft')

vital_files = []
vital_files.append('configure')


# list directories and files that are not compiled as part of main executable
# ===========================================================================

directories_not_built = []
directories_not_built.append('utils')

files_not_built = []
files_not_built.append('utils/cube2pov.c')
files_not_built.append('test/MPI2_test.F')
files_not_built.append('relccsd/ccdriv_nomemget.F')
files_not_built.append('relccsd/ccints_nomemget.F')
files_not_built.append('eri/eri2par.F')
files_not_built.append('openrsp/aoresponse_visual.F90')
files_not_built.append('densfit/denfit_noddy.F')
files_not_built.append('lucita/unused.F')
files_not_built.append('dft/vrg.F')
files_not_built.append('eri/eriprg.F')
files_not_built.append('dft/fun-tester.c')
files_not_built.append('relccsd/ccmain.F')
files_not_built.append('relccsd/quit.F')
files_not_built.append('build/gen_fortran_interface_h.c')
files_not_built.append('dirrci/cimain.F')


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
dirac_modules['eri'].add_directory('eri')

dirac_modules['densfit'] = dirac_module('densfit', 'Include density fitting module')
dirac_modules['densfit'].add_cppflag('MOD_DNF')
dirac_modules['densfit'].add_dependency('eri')
dirac_modules['densfit'].add_directory('densfit')

dirac_modules['largemolecules'] = dirac_module('largemolecules', 'Increases common blocks to handle more atoms and basis sets')
dirac_modules['largemolecules'].add_cppflag('LARGE_MOLECULES')
dirac_modules['largemolecules'].take_out_of_adm()

dirac_modules['fmm'] = dirac_module('fmm', 'Include fast multipole method')
dirac_modules['fmm'].add_cppflag('MOD_FMM')
dirac_modules['fmm'].add_dependency('eri')
dirac_modules['fmm'].add_dependency('densfit')
dirac_modules['fmm'].add_directory('fmm')
dirac_modules['fmm'].take_out_of_adm()

dirac_modules['esr'] = dirac_module('esr', 'Include ESR property module')
dirac_modules['esr'].add_cppflag('MOD_ESR')
dirac_modules['esr'].add_directory('esr')

dirac_modules['xpp'] = dirac_module('xpp', 'Include PP excitation energy and STEX modules')
dirac_modules['xpp'].add_cppflag('MOD_XPP')
dirac_modules['xpp'].add_cppflag('MOD_STEX')
dirac_modules['xpp'].add_directory('xpp')

dirac_modules['xml'] = dirac_module('xml', 'Include XML input/output module')
dirac_modules['xml'].add_cppflag('MOD_XML')
dirac_modules['xml'].add_directory('gp/xml')
dirac_modules['xml'].take_out_of_adm()

dirac_modules['krmc'] = dirac_module('krmc', 'Include KR-MCSCF and LUCIAREL modules')
dirac_modules['krmc'].add_cppflag('MOD_KRMC')
dirac_modules['krmc'].add_directory('krmc')
dirac_modules['krmc'].add_directory('luciarel')

dirac_modules['emb'] = dirac_module('emb', 'Include ab initio/DFT-embedding module')
dirac_modules['emb'].add_cppflag('MOD_EMB')
dirac_modules['emb'].add_directory('embedding')

dirac_modules['openrsp'] = dirac_module('openrsp', 'Include OpenRSP')
dirac_modules['openrsp'].add_cppflag('MOD_OPENRSP')
dirac_modules['openrsp'].add_directory('openrsp')

dirac_modules['xcfun'] = dirac_module('xcfun', 'Include XCFun module')
dirac_modules['xcfun'].add_cppflag('MOD_XCFUN')
dirac_modules['xcfun'].add_directory('xcfun')
dirac_modules['xcfun'].take_out_of_adm()

dirac_modules['cap'] = dirac_module('cap', 'Include CAP module')
dirac_modules['cap'].add_cppflag('MOD_CAP')
dirac_modules['cap'].add_dependency('oof')
dirac_modules['cap'].add_directory('cap')
dirac_modules['cap'].take_out_of_adm()

#fixme: in pamadm this flag triggers a lot of operations
#       so far this is not done here
#       (can they be avoided?)
#       same for directories ... so complicating
dirac_modules['hscc'] = dirac_module('hscc', 'Include interface for the (external) HSMRCC code together with its objects')
dirac_modules['hscc'].add_cppflag('MOD_HSCC')
dirac_modules['hscc'].add_dependency('oof')
dirac_modules['hscc'].add_directory('hsmrcc')
dirac_modules['hscc'].take_out_of_adm()

#fixme: in pamadm this flag triggers a lot of operations
#       so far this is not done here
#       (can they be avoided?)
#       same for directories ... so complicating
dirac_modules['hsfs'] = dirac_module('hsfs', 'Include higher external Higher sector FSCC code')
dirac_modules['hsfs'].add_cppflag('MOD_HSFS')
dirac_modules['hsfs'].add_dependency('oof')
dirac_modules['hsfs'].add_dependency('hscc')
dirac_modules['hsfs'].add_directory('hsmrcc')
dirac_modules['hsfs'].take_out_of_adm()

#fixme: in pamadm this flag triggers a lot of operations
#       so far this is not done here
#       (can they be avoided?)
#       same for directories ... so complicating
dirac_modules['hsfsdyn'] = dirac_module('hsfsdyn', 'Include higher external Higher sector FSCC code - dynamic allocation version.')
dirac_modules['hsfsdyn'].add_cppflag('HSFS_DYN_ALLOC')
dirac_modules['hsfsdyn'].add_dependency('oof')
dirac_modules['hsfsdyn'].add_dependency('hsfs')
dirac_modules['hsfsdyn'].add_dependency('hscc')
dirac_modules['hsfsdyn'].add_directory('hsmrcc')
dirac_modules['hsfsdyn'].take_out_of_adm()

dirac_modules['qdots'] = dirac_module('qdots', 'Include Quantum Dots functionality (artificial atom)')
dirac_modules['qdots'].add_cppflag('MOD_QDOTS')
dirac_modules['qdots'].add_dependency('oof')
dirac_modules['qdots'].take_out_of_adm()

dirac_modules['ihfscc'] = dirac_module('ihfscc', 'Include Intermediate Hamiltonian functionality to FSCC code')
dirac_modules['ihfscc'].add_cppflag('MOD_IHFSCC')

dirac_modules['amfi'] = dirac_module('amfi', 'Include the AMFI functionality into DIRAC')
dirac_modules['amfi'].add_cppflag('MOD_AMFI')
dirac_modules['amfi'].add_directory('amfi')

dirac_modules['london'] = dirac_module('london', 'Include London orbitals functionality for magnetic properties')
dirac_modules['london'].add_cppflag('MOD_LONDON')

dirac_modules['oof'] = dirac_module('oof', 'Include object oriented modularity of Fotran 90/95')
dirac_modules['oof'].add_cppflag('MOD_OOF')
dirac_modules['oof'].add_directory('oof')
dirac_modules['oof'].take_out_of_adm()

dirac_modules['adc'] = dirac_module('adc', 'Include ADC electron propagator module')
dirac_modules['adc'].add_cppflag('MOD_ADC')
dirac_modules['adc'].add_directory('reladc')

dirac_modules['ecp'] = dirac_module('ecp', 'Include ARGOS ECP module')
dirac_modules['ecp'].add_cppflag('MOD_ECP')
dirac_modules['ecp'].add_directory('argosecp')
dirac_modules['ecp'].take_out_of_adm()

dirac_modules['srdft'] = dirac_module('srdft', 'Include srDFT module')
dirac_modules['srdft'].add_cppflag('MOD_SRDFT')
dirac_modules['srdft'].add_cppflag('MOD_KRMC')
dirac_modules['srdft'].add_directory('srdft')
dirac_modules['srdft'].take_out_of_adm()

dirac_modules['unreleased'] = dirac_module('unreleased', 'Include unreleased functionality (not covered by other flags)')
dirac_modules['unreleased'].add_cppflag('MOD_UNRELEASED')

dirac_modules['log'] = dirac_module('log', 'Will log interesting matrices in matlab format (very verbose!)')
dirac_modules['log'].add_cppflag('MOD_LOG')
dirac_modules['log'].take_out_of_adm()


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
                  help="Include (almost) all modules (administrator)")


# parse input
# ===========

(options, args) = parser.parse_args()
if options.adm:
   for module in dirac_modules:
       if dirac_modules[module].part_of_adm:
          dirac_modules[module].include()


# check dependencies among modules
# ================================

for module in dirac_modules:
    if dirac_modules[module].is_included:
       for need_module in dirac_modules[module].depends_on:
           if not dirac_modules[need_module].is_included:
              error_message = '%s depends on %s but %s is not included' % (module, need_module, need_module)
              sys.exit(error_message)


# report on screen
# ================

print 'Active modules:'
for module in dirac_modules:
    if dirac_modules[module].is_included:
       print dirac_modules[module].message
print


# prepare dontbuild
# =================

f = open('build/dontbuild', 'w')

for module in dirac_modules:
    if not dirac_modules[module].is_included:
       for directory in dirac_modules[module].directories:
           s = '%s/.*' % directory
           f.write(s + '\n')
for directory in directories_not_built:
    s = '%s/.*' % directory
    f.write(s + '\n')
for file in files_not_built:
    s = '%s' % file
    f.write(s + '\n')

f.close()


# touch Makefile.dep and Makefile.obj
# ===================================

os.system('touch -t 198001010000 build/Makefile.dep')
os.system('touch -t 198001010000 build/Makefile.obj')


# prepare Makefile.extobj and Makefile.moddef
# ===========================================

#fixme: so far doing nothing with Makefile.extobj
os.system('touch -t 198001010000 build/Makefile.extobj')

f = open('build/Makefile.moddef', 'w')

f.write('# file generated by the minions (pamadm) of Dirac - DO NOT EDIT!\n')
f.write('#\n')

flags = 'MODS='
for module in dirac_modules:
    if dirac_modules[module].is_included:
       s = '# %s' % dirac_modules[module].message
       for cppflag in dirac_modules[module].cppflags:
           s += ' -D%s' % cppflag
           flags += '-D%s ' % cppflag
       f.write(s + '\n')
f.write(flags + '\n')

# determine subversion url and revision
#fixme: so far no svn determ
f.write("DIRAC_SVN_URL='unknown Subversion URL'\n")
f.write("DIRAC_SNAP='unknown Subversion revision'\n")

f.close()


# create and distribute Makefile over active directories
# ======================================================



# export dirac to tar.gz file
# ===========================

if options.tarfile != None:
   print 'export', options.tarfile





