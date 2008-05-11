#!/usr/bin/env python
#
# experimental alternative for pamadm
# use at own risk
# 2008-05-10 radovan bast
#
# for syntax:
# ./alternative_pamadm.py -h
#
# presently it builds only the targets dirac.x and dirmpi.x
# so no utils, no standalones
# also not so nice: have to introduce a blacklist of objects
# which are not built
# otherwise there are conflicts
# the (ordered) list of directories which contain f90 modules is given by hand
# it's difficult to let the script evaluate mutual dependencies
#
# not sure about portability (tested only two machines)
#
# currently no release/export options
# no backup

#-------------------------------------------------------------------------------

import string, sys, os, fnmatch
from optparse import OptionParser

#-------------------------------------------------------------------------------

class dirac_module_c:

      def __init__(self):
          self.help        = ''
          self.flag        = ''
          self.activated   = False
          self.requires    = []
          self.directories = []

#-------------------------------------------------------------------------------
# rules

dirac_module_d = {}

name = 'dft'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include DFT module'
dirac_module_d[name].flag        = '-DMOD_DFT'
dirac_module_d[name].directories = ['dft']

name = 'esr'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include ESR property module'
dirac_module_d[name].flag        = '-DMOD_ESR'
dirac_module_d[name].directories = ['prp/esr']

name = 'eri'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include ERI integral routines'
dirac_module_d[name].flag        = '-DMOD_ERI'
dirac_module_d[name].directories = ['eri']

name = 'xpp'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include PP excitation energy and STEX modules'
dirac_module_d[name].flag        = '-DMOD_XPP -DMOD_STEX'
dirac_module_d[name].directories = ['prp/xpp']

name = 'fmm'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include fast multipole method'
dirac_module_d[name].flag        = '-DMOD_FMM'
dirac_module_d[name].directories = ['fmm']
dirac_module_d[name].requires    = ['eri', 'densfit']

name = 'densfit'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include density fitting module'
dirac_module_d[name].flag        = '-DMOD_DNF'
dirac_module_d[name].directories = ['densfit']
dirac_module_d[name].requires    = ['eri']

name = 'krmc'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include KR-MCSCF and LUCIAREL modules'
dirac_module_d[name].flag        = '-DMOD_KRMC'
dirac_module_d[name].directories = ['krmc', 'luciarel']

name = 'embedding'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include ab initio/DFT-embedding module'
dirac_module_d[name].flag        = '-DMOD_EMB'
dirac_module_d[name].directories = ['embedding']

name = 'adc'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include ADC electron propagator module'
dirac_module_d[name].flag        = '-DMOD_ADC'
dirac_module_d[name].directories = ['reladc']

name = 'ecp'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include ARGOS ECP module'
dirac_module_d[name].flag        = '-DMOD_ECP'
dirac_module_d[name].directories = ['argosecp']

name = 'cdft'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include SDFT/CDFT module'
dirac_module_d[name].flag        = '-DMOD_CDFT'
dirac_module_d[name].directories = ['cdft']

name = 'visual'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include visualization module'
dirac_module_d[name].flag        = '-DMOD_VISUAL'
dirac_module_d[name].directories = ['visual']

name = 'aoresponse'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include AO response module'
dirac_module_d[name].flag        = '-DMOD_AORESPONSE'
dirac_module_d[name].directories = ['aoresponse']

name = 'cap'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include CAP module'
dirac_module_d[name].flag        = '-DMOD_CAP'
dirac_module_d[name].directories = ['cap']
dirac_module_d[name].requires    = ['oof']

name = 'hscc'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include interface for the (external) HSMRCC code'
dirac_module_d[name].flag        = '-DMOD_HSCC'
dirac_module_d[name].directories = ['hsmrcc', 'hsmrcc/second_quant']
dirac_module_d[name].requires    = ['oof']

name = 'oof'
dirac_module_d[name]             = dirac_module_c()
dirac_module_d[name].help        = 'include object oriented modularity of Fotran 90/95'
dirac_module_d[name].flag        = '-DMOD_OOF'
dirac_module_d[name].directories = ['oof']

#-------------------------------------------------------------------------------
# directories which need sfmakedepend
# this is not so nice to maintain but we will go away
# from sfmakedepend anyway sooner or later

need_sfmakedepend_l = []

need_sfmakedepend_l.append('abacus')
need_sfmakedepend_l.append('gp')
need_sfmakedepend_l.append('utils')
need_sfmakedepend_l.append('embedding')
need_sfmakedepend_l.append('visual')
need_sfmakedepend_l.append('cap')
need_sfmakedepend_l.append('hsmrcc')
need_sfmakedepend_l.append('hsmrcc/second_quant')
need_sfmakedepend_l.append('oof')

#-------------------------------------------------------------------------------
# directories which contain f90 modules
# the order matters, here we need a sophisticated sctipt
# which figures out the mutual dependencies

contains_modules_l = []

contains_modules_l.append('abacus')
contains_modules_l.append('gp')
contains_modules_l.append('visual')
contains_modules_l.append('aoresponse')
contains_modules_l.append('embedding')
contains_modules_l.append('oof')
contains_modules_l.append('cap')
contains_modules_l.append('hsmrcc')

#-------------------------------------------------------------------------------
# parse arguments

parser = OptionParser()

def little_helper(option, opt_str, value, parser):
    dirac_module_d[opt_str[2:]].activated = True

for dirac_module in dirac_module_d:
    flag = '--' + dirac_module
    parser.add_option(flag, action = 'callback', callback = little_helper, help = dirac_module_d[dirac_module].help)

parser.add_option('--adm',   action = 'store_true', default = False, help = 'include ALL modules (administrator)')
parser.add_option('--ibm',   action = 'store_true', default = False, help = 'prepare build for IBM architecture')
parser.add_option('--large', action = 'store_true', default = False, help = 'increase commom blocks to handle more atoms and basis sets')

(options, args) = parser.parse_args()

#-------------------------------------------------------------------------------
# consistency

# modules activated by -adm
if options.adm:

   dirac_module_d['dft'       ].activated = True
   dirac_module_d['esr'       ].activated = True
   dirac_module_d['eri'       ].activated = True
   dirac_module_d['xpp'       ].activated = True
#  dirac_module_d['fmm'       ].activated = True
   dirac_module_d['densfit'   ].activated = True
   dirac_module_d['krmc'      ].activated = True
#  dirac_module_d['embedding' ].activated = True
#  dirac_module_d['adc'       ].activated = True
#  dirac_module_d['ecp'       ].activated = True
   dirac_module_d['cdft'      ].activated = True
   dirac_module_d['visual'    ].activated = True
#  dirac_module_d['aoresponse'].activated = True
#  dirac_module_d['cap'       ].activated = True
#  dirac_module_d['hscc'      ].activated = True
#  dirac_module_d['oof'       ].activated = True

# loop over all modules, if module activated, activate also list module.requires
for dirac_module in dirac_module_d:
    if dirac_module_d[dirac_module].activated:
       for required_module in dirac_module_d[dirac_module].requires:
           dirac_module_d[required_module].activated = True

#-------------------------------------------------------------------------------
# collect flags

flags    = ''
nr_flags = 0
for dirac_module in dirac_module_d:
    if dirac_module_d[dirac_module].activated:
       nr_flags += 1
       if nr_flags > 1:
          if options.ibm:
             flags += ','
          else:
             flags += ' '
       flags += dirac_module_d[dirac_module].flag

if options.large:
   if options.ibm:
      flags += ','
   else:
      flags += ' '
   flags += '-DLARGE_MOLECULES'

#-------------------------------------------------------------------------------
# collect directories

directory_l = ['dirac', 'prp', 'pdpack', 'abacus', 'amfi', 'gp', 'main', 'relccsd', 'moltra', 'dirrci', 'lucita', 'utils']
for dirac_module in dirac_module_d:
    if dirac_module_d[dirac_module].activated:
       for directory in dirac_module_d[dirac_module].directories:
           directory_l.append(directory)

#-------------------------------------------------------------------------------
# for directory in directory_l generate Makefile from Makefile.in

for directory in directory_l:

    f = open(directory + '/Makefile.in', 'r')
    s = f.read()
    f.close()

    f = open(directory + '/Makefile', 'w')
    f.write('MODS=' + flags + '\n' + s)
    f.close()
    print 'generated ' + directory + '/Makefile'

#-------------------------------------------------------------------------------
# collect objects

blacklist = []

blacklist.append('pdpack/arhpack.o')
blacklist.append('pdpack/gpblas.o')
blacklist.append('pdpack/linextra.o')
blacklist.append('pdpack/gp_dlapack.o')
blacklist.append('pdpack/eispack.o')
blacklist.append('pdpack/dge.o')
blacklist.append('pdpack/dsp.o')
blacklist.append('pdpack/printpkg.o')
blacklist.append('pdpack/gp_zlapack.o')
blacklist.append('prp/oldrvc.o')
blacklist.append('abacus/hercen.o')
blacklist.append('gp/graphic.o')
blacklist.append('relccsd/ccints_nomemget.o')
blacklist.append('relccsd/ccdriv_nomemget.o')
blacklist.append('relccsd/quit.o')
blacklist.append('lucita/doubledir.o')
blacklist.append('lucita/unused.o')
blacklist.append('utils/vibcal.o')
blacklist.append('utils/mx2fit.o')
blacklist.append('utils/twofit.o')
blacklist.append('utils/polfit.o')
blacklist.append('utils/qjaco.o')
blacklist.append('utils/rspread.o')
blacklist.append('utils/labread.o')
blacklist.append('eri/eri2par.o')
blacklist.append('eri/eriprg.o')
blacklist.append('densfit/denfit_noddy.o')
blacklist.append('utils/test_allocator.o')
blacklist.append('densfit/df_dofit.o')
blacklist.append('relccsd/ccunix.o')
blacklist.append('lucita/flush_stdout.o')
blacklist.append('utils/basis_set.o')
blacklist.append('dft/fun-tester.o')
blacklist.append('dft/fun-pw92.o')
blacklist.append('relccsd/ccmain.o')
blacklist.append('dirrci/cimain.o')

f77_object_l = []
f90_object_l = []
c_object_l   = []

for directory in directory_l:
    for file in os.listdir(directory):

        if fnmatch.fnmatch(file, '*.F'):
           object = directory + '/' + string.replace(file, '.F', '.o')
           if not object in blacklist:
              f77_object_l.append(object)

        if fnmatch.fnmatch(file, '*.F90'):
           object = directory + '/' + string.replace(file, '.F90', '.o')
           if not object in blacklist:
              f90_object_l.append(object)

        if fnmatch.fnmatch(file, '*.c'):
           object = directory + '/' + string.replace(file, '.c', '.o')
           if not object in blacklist:
              c_object_l.append(object)

#-------------------------------------------------------------------------------

makefile  = 'MODS=' + flags + '\n'
makefile += 'include ./Makefile.config' + '\n\n'

makefile += 'F77_OBJECTS = \\' + '\n'
for i in range(len(f77_object_l)):
    makefile += f77_object_l[i]
    if i == len(f77_object_l) - 1:
       makefile += '\n\n'
    else:
       makefile += ' \\\n'

makefile += 'F90_OBJECTS = \\' + '\n'
for i in range(len(f90_object_l)):
    makefile += f90_object_l[i]
    if i == len(f90_object_l) - 1:
       makefile += '\n\n'
    else:
       makefile += ' \\\n'

makefile += 'C_OBJECTS = \\' + '\n'
for i in range(len(c_object_l)):
    makefile += c_object_l[i]
    if i == len(c_object_l) - 1:
       makefile += '\n\n'
    else:
       makefile += ' \\\n'

#-------------------------------------------------------------------------------

makefile += 'PDLIB = -Lpdpack -lpdlib\n\n'
makefile += 'TOTAL_F_OBJECTS = $(F77_OBJECTS) $(F90_OBJECTS)\n\n'
makefile += 'F_CPPS = $(TOTAL_F_OBJECTS:.o=.i)\n\n'

makefile += 'dirac.x: ALL_F90MODS'
makefile += ' \\\n\tEXTRA_STUFF'
for directory in directory_l:
    if not directory == 'utils':
       directory_object = string.replace(directory, '/', '') + '_object'
       makefile += ' \\\n\t' + directory_object.upper()
makefile += '\n\t$(F77) $(LDFLAGS) -o $(INSTALLDIR)/dirac.x \\\n'
makefile += '\t$(F77_OBJECTS) $(F90_OBJECTS) $(C_OBJECTS) \\\n'
makefile += '\t$(PDLIB) $(LIB_LIST)\n\n'

makefile += 'pre: $(F_CPPS)\n\n'

makefile += 'ftnchek: pre\n'
makefile += '\tftnchek $(CHEKFLAGS) $(F_CPPS) > dirac.ftnchek\n\n'

makefile += 'dirmpi.x: ALL_F90MODS'
makefile += ' \\\n\tEXTRA_STUFF'
for directory in directory_l:
    if not directory == 'utils':
       directory_object = string.replace(directory, '/', '') + '_object'
       makefile += ' \\\n\t' + directory_object.upper()
makefile += '\n\t$(F77) $(LDFLAGS) -o $(INSTALLDIR)/dirac.x \\\n'
makefile += '\t$(F77_OBJECTS) $(F90_OBJECTS) $(C_OBJECTS) \\\n'
makefile += '\t$(PDLIB) $(LIB_LIST) $(MPI_LIB_PATH) $(MPI_LIB)\n\n'

#-------------------------------------------------------------------------------
# target clean

makefile += 'clean:\n'

for directory in directory_l:
    makefile += '\tcd ' + directory + '; $(MAKE) clean\n'

makefile += '\t$(RM) -f *~\n\n'

#-------------------------------------------------------------------------------
# target mod

makefile += 'mod:\n'

for directory in directory_l:
    if directory in contains_modules_l:
       makefile += '\tcd ' + directory + ' && $(MAKE) mod &&  cd .. ;\n'

#-------------------------------------------------------------------------------

makefile += 'ALL_F90MODS:\n'
makefile += '\t$(MAKE) mod\n\n'

#-------------------------------------------------------------------------------
# directory objects

for directory in directory_l:
    directory_object = string.replace(directory, '/', '') + '_object'
    makefile += directory_object.upper() + ':\n'
    makefile += '\tcd ' + directory + '; $(MAKE) all\n\n'

makefile += 'EXTRA_STUFF:\n'
makefile += '\tcd gp && $(MAKE) mxm && cd .. ;\n'
makefile += '\tcd abacus && $(MAKE) slave && cd .. ;\n'

#-------------------------------------------------------------------------------
# target depend

makefile += 'depend:\n'

for directory in directory_l:
    makefile += '\tcd ' + directory + '; touch Makefile.depend; $(MAKE) depend\n'
    if directory in need_sfmakedepend_l:
       makefile += '\tcd ' + directory + '; $(MAKE) f90dep\n'

#-------------------------------------------------------------------------------

f = open('Makefile', 'w')
f.write(makefile)
f.close()

print 'generated Makefile'

print '\nrun:'
print '  make depend'
print 'followed by:'
print '  make'
