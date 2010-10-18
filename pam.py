#! /usr/bin/env python

##################################################
#
# Python replacement of the pam script
#
#   Reads all data from Makefile.config  !
#
##################################################

import  sys, os, string, subprocess, optparse, distutils.util

# installation directory inserted in by the configure script
global install_dir
install_dir='/home/ilias/qch_work/qch_software/dirac/gfortran_noopt_ILP64_ownlibs/adm/trunk'

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class pam_enviroment:
  # ... class containing all important variables, methods for submitting
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def __init__(self):
  #
  #           ... initialize the global set of variables ...
  #
     # print level
    self.iprint=1
     ### contains installation directory ...
    self.dirac_home=install_dir
    self.parallel = False
    self.mem = 64
    self.slavemem = 64

  def read_pam_variables():
  # important - reads variables from Makefile.config
    print 'reading pam variables...'

#def notify_user() :
# replacement
#  print 'notify_user...'

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#class pam_run(pam_variables):
class pam_run:
#-----------------------------------------------------
#-----------------------------------------------------
  def __init__(self):
    print 'pam run initialize'

  def perform(self):
    print 'performing own pam run..'
#set_work_dir  
#set_dirac_environment
#print_header
#copy_files
#save_old_output
#setup_command
#distribute_work_dir
#copy_files_master
#execute_dirac
#profile_dirac

  def process_pam_arguments():
    print 'processing pam arguments....'
    from optparse import OptionParser
    usage="usage: %prog arg1 arg2..."
    parser = OptionParser(usage)
  # ... get in the print level 
    iprint_label='script verbose print level, default:'+str(1)
    parser.add_option("-i", "--iprint", type="int", nargs=1,  
                  action="store", dest="dest_iprint", default=1,
                  help=iprint_label, metavar="INT_NUM")
    (options, args) = parser.parse_args()


###########################################################################
###                   process the main function
###########################################################################
def main(*args):
   print "pam.py main function..."
   run = pam_run()
   run.perform()

if __name__ == '__main__':
   sys.exit(main(*sys.argv))


