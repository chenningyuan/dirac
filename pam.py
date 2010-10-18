#! /usr/bin/env python

##################################################
#
# Python replacement of the pam script
#
#   Reads all data from Makefile.config  !
#
##################################################

import  sys, os, string, subprocess, optparse, distutils.util, shutil

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
    print 'pam_enviroment class initialize'
     # print level
    self.iprint=1
     ### contains installation directory ...
    self.dirac_home=install_dir
    self.basdir=''
    self.wrkdir='/tmp/ilias/xxx'
    self.wrkpath='/tmp/ilias'
    self.parallel = False
    self.mpi=1
    self.keepwrk=False
    self.dirac_exe=self.dirac_home
    self.mem = 64
    self.slavemem = 64
    self.full_restart=False
    self.debugger=''
    self.noarch=False
    self.input_molecule=''
    self.input_dirac=''

  #def read_pam_variables(self):
  # important - reads variables from Makefile.config
    #print 'reading pam variables...'

  def process_pam_arguments(self):
    print 'processing pam arguments....'

    def process_remaining_args(rem_args):
    # checks and process *.mol & *.inp arguments
      print rem_args
      print len(rem_args)
      if (len(rem_args) == 2):
     # add further checking on input files
        self.input_molecule=rem_args[0]
        self.input_dirac=rem_args[1]
      else:
        print 'wrong number of arguments !'
        sys.exit(-2)

    from optparse import OptionParser
    usage="usage: %prog arg1 arg2..."
    parser = OptionParser(usage)
  # ... get in the print level 
    this_label='script verbose print level, default:'+str(1)
    parser.add_option("-i", "--iprint", type="int", nargs=1,  
                  action="store", dest="dest_iprint", default=1,
                  help=this_label, metavar="INT_NUM")
  # ... get in the  mw
    this_label='memory (in megawords) for the master, default:'+str(self.mem)
    parser.add_option("--mw", type="int", nargs=1,  
                  action="store", dest="dest_mw", default=1,
                  help=this_label, metavar="INT_NUM")

  # ... get in the  nw
    this_label='memory (in megawords) for nodes, default:'+str(self.slavemem)
    parser.add_option("--nw", type="int", nargs=1,  
                  action="store", dest="dest_nw", default=1,
                  help=this_label, metavar="INT_NUM")

    # process all
    (options, args) = parser.parse_args()
    
    # now assign individual arguments
    if (options.dest_iprint is not None) :
     self.iprint=int(options.dest_iprint)
     if self.iprint > 7:
      print 'assigned iprint=',self.iprint

    if (options.dest_mw is not None) :
     self.iprint=int(options.dest_mw)
     if self.iprint > 7:
      print 'assigned mw=',self.mem

    # remaining arguments - *.mol and *.inp files
    
    process_remaining_args(args)
    
  def execute_dirac(self):
    print 'executing dirac ...'
    if not os.path.isdir(self.wrkdir):
      print 'creating directory'
      os.mkdir(self.wrkdir)
    else:
      print 'directory exists'
    #copy executable
    src=install_dir+'/'+'dirac.x'
    dest_dir=self.wrkdir+'/'+'dirac.x'
    shutil.copy(src, dest_dir)
    print 'dirac copied'
    #copy input files
    src=self.input_molecule
    dest=self.wrkdir+'/'+'MOLECULE.MOL'
    shutil.copy(src, dest)
    #
    src=self.input_dirac
    dest=self.wrkdir+'/'+'DIRAC.INP'
    shutil.copy(src, dest)
    for filename in os.listdir(self.wrkdir):
      print  filename
   
    # launch dirac
    command='cd '+self.wrkdir+' ; ./dirac.x'
    print command
    os.system(command)
   

#def notify_user() :
# replacement
#  print 'notify_user...'

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#class pam_run(pam_variables):
class pam_run(pam_enviroment):
#-----------------------------------------------------
#-----------------------------------------------------
  #pass
  #global pam_variables
  #pam_variables=pam_enviroment()

  def __init__(self):
    print 'pam_run class initialize'
    self.pam_variables=pam_enviroment()

  def perform(self):
    print 'performing own pam run..'
    self.pam_variables.process_pam_arguments()
    self.pam_variables.execute_dirac()

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
def main(*args):
   print "pam.py main function..."
   run = pam_run()
   run.perform()

if __name__ == '__main__':
   sys.exit(main(*sys.argv))


