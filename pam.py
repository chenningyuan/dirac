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
class pam_variables:
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
    print 'readin pam variables...'
#def notify_user() :
# replacement
#  print 'notify_user...'

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_pam_arguments():
#-----------------------------------------------------
#
#-----------------------------------------------------
   print 'processing pam arguments....'
# from optparse import OptionParser
# usage="usage: %prog arg1 arg2..."
# parser = OptionParser(usage)

def main(*args):
   print "pam.py main function..."
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

###########################################################################
###                   process the main function
###########################################################################
if __name__ == '__main__':
   sys.exit(main(*sys.argv))

