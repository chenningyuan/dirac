#! /usr/bin/env python

##################################################
#
# Python replacement of the pam script
#
# Reads crucial data from Makefile.config
#
#  Written by Miro ILIAS,Banska Bystrica, 2010
#
#
##################################################

import  sys, os, string, subprocess, optparse, distutils.util, shutil, signal

# installation directory inserted in by the configure script
global install_dir, install_dir_default
#install_dir_default='/home/ilias/qch_work/qch_software/dirac/gfortran_noopt_ILP64_ownlibs/adm/trunk'
install_dir_default='/home/dv6-2145ec/Work/qch/software/dirac/gfortran_gcc_noopt_ownlibs_ilp32/adm_hsfsdyn_cap/trunk'
install_dir=''

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class workdir_class:
# class dealing with the work directory
  global install_dir
  def __init__(self):
    print 'work_direcory class initialize'
    self.workdir=''

    #self.workpath='/tmp/ilias'
    #self.workpath='/tmp/'+pwd.getpwuid(os.getuid())[0]
    self.workpath='/tmp/'+os.getlogin()
    #TODO: Windows user login name !!!
    self.keepwrk=False

    self.scratchdir=''
    self.iprint=0
    self.create_workdir_name=True

  def set_workdir_name(self,root_name):
  # sets the working directory name according to 
  # entering "root_name"
    print 'set the name of the work_direcory',root_name
    current_pid = os.getpid(); print 'pid=',current_pid
    workdir_name=root_name+"."+str(current_pid)
    self.workdir=self.workpath+"/"+workdir_name
    print 'workdir_name=',self.workdir

  def create_workdir(self):
  # creates the working directory
    if not os.path.isdir(self.workdir):
      print 'creating work directory '+self.workdir
      os.mkdir(self.workdir)
    else:
      print 'directory exists'

  def copy_file(self,src,dest):
  # ... copies file from the current dir "scr" to "dest"
    print 'copying file...',src,' to ',dest
    shutil.copy(src, dest)

  def list_workdir(self):
  # ... lists the content of the work-directory
  # more here http://www.daniweb.com/code/snippet216798.html
    for filename in os.listdir(self.workdir):
      print  filename,os.path.getsize(self.workdir+"/"+filename)

  def remove_workdir(self):
    print 'removes the work_direcory'
    self.list_workdir()

    if (not self.keepwrk):
      print 'deleting the work directory'
      shutil.rmtree(self.workdir)      
    else:
      print 'keeping the work directory'

  def copy_file_to_workdir(self,some_file):
    print 'copying file to the work_direcory'

  def get_file_from_workdir(self,some_file):
    print 'getting file to the work_direcory'


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class pam_enviroment(workdir_class):
  global install_dir
  # ... class containing all important variables, methods for submitting
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def __init__(self):
  #
  #           ... initialize the global set of variables ...
  #
    print 'pam_enviroment class initialize'

    #global install_dir

     # print level
    self.iprint=1
     ### contains installation directory ...
    self.dirac_home=install_dir
    self.basdir=''

    # ... object from the inherited "workdirectory" class 
    self.wrkdir_obj = workdir_class()
    # ... output file
    self.output_file=''

    self.parallel = False

    self.profile = False

    self.profiler = ''
    self.mpi=1

    self.dirac_exe=self.dirac_home

    self.mem = 64
    self.slavemem = 64

    self.full_restart=False
    self.debugger=''
    self.noarch=False
    # inpu files
    self.input_molecule=''
    self.input_dirac=''
    # output name from *inp *mol names
    self.output_name=''
    # output suffix default
    self.output_suffix="out"
    # output file...
    self.output_file=''

    #execution time-limit
    self.set_timeout_limit=False
    self.timeout=''
	
  def read_pam_variables(self):
   #important - reads variables from Makefile.config
    print 'reading pam variables from Makefile.config...'

  def process_pam_arguments(self):
    print 'processing pam arguments....'
    global install_dir

    def process_remaining_args(rem_args):
    # checks and process *.mol & *.inp arguments
    # which are the only remaining arguments after
    # picking up all others by flags
      print rem_args
      print len(rem_args)
      if (len(rem_args) == 2):
     # add further checking on input files
        self.input_molecule=rem_args[0]
        self.input_dirac=rem_args[1]
        self.set_output_name()
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
                  action="store", dest="dest_mw", default=64,
                  help=this_label, metavar="INT_NUM")

    # ... get in the  nw
    this_label='memory (in megawords) for nodes, default:'+str(self.slavemem)
    parser.add_option("--nw", type="int", nargs=1,  
                  action="store", dest="dest_nw", default=64,
                  help=this_label, metavar="INT_NUM")

    # ... get in the timeout limit
    this_label='limit dirac.x execution time (in seconds)'
    parser.add_option("--timeout", type="int", nargs=1,  
                  action="store", dest="dest_timeout",
                  help=this_label, metavar="INT_NUM")

    # .. read in the install_dir, what is global variable...
    this_label="pam.py\'s directory, default:"+install_dir_default
    parser.add_option("--install_dir", type="string", nargs=1,
                  action="store",dest="dest_install_dir",  
                  help=this_label,metavar="STRING")

    # ... read in timeout limit for dirac execution
       

    # process all
    (options, args) = parser.parse_args()
    
    # now assign individual arguments
    if (options.dest_iprint is not None) :
     self.iprint=int(options.dest_iprint)
     if self.iprint > 7:
      print 'assigned iprint=',self.iprint

    if (options.dest_mw is not None) :
     self.mem=int(options.dest_mw)
     if self.iprint > 7:
      print 'assigned mw=',self.mem

    if (options.dest_nw is not None) :
     self.memslave=int(options.dest_mw)
     if self.iprint > 7:
      print 'assigned mw=',self.memslave

    if (options.dest_install_dir is not None) :
     install_dir=str(options.dest_install_dir)
     print 'assigned install_dir=',install_dir
     if self.iprint > 7:
      print 'assigned install-dir=',install_dir
    elif (install_dir_default is None) :
     sys.exit("you must set the installation directory !")

    if (options.dest_timeout is not None) :
     self.timeout=int(options.dest_timeout)
     self.set_timeout_limit=True
     if self.iprint > 7:
      print 'assigned timeout=',self.timeout

    # remaining arguments - *.mol and *.inp files
    process_remaining_args(args)

  def set_wrkdir(self):
  # ... creates the workdir, copies wanted files
    if self.wrkdir_obj.create_workdir_name:
      self.wrkdir_obj.set_workdir_name(self.output_name)
    self.wrkdir_obj.create_workdir()
    print 'set_wrkdir: install_dir=',install_dir
    self.wrkdir_obj.copy_file(install_dir+'/'+'dirac.x',self.wrkdir_obj.workdir+'/'+'dirac.x')
    self.wrkdir_obj.copy_file(self.input_molecule,self.wrkdir_obj.workdir+'/'+'MOLECULE.MOL')
    self.wrkdir_obj.copy_file(self.input_dirac,self.wrkdir_obj.workdir+'/'+'DIRAC.INP')
   # TODO: copying more file through flag)

  def set_output_name(self):
  # set the output name string to be used later 
  # for naming output and work directory
    #mol_file_stripped=self.input_molecule.strip(".mol" )
    mol_file_stripped=self.input_molecule.replace(".mol","",1 )
    #inp_file_stripped=self.input_dirac.strip(".inp")
    inp_file_stripped=self.input_dirac.replace(".inp","",1)
    print mol_file_stripped, inp_file_stripped
    if (mol_file_stripped == inp_file_stripped):
      self.output_name=mol_file_stripped
    else:
      self.output_name=mol_file_stripped+"_"+inp_file_stripped
    print "output_name=",self.output_name

  def set_output_file(self):
  # sets the name of the output file by using "self.output_name"
    print 'setting output file..'
    self.output_file=self.output_name+"."+self.output_suffix
    print 'output file..',self.output_file
 
  def execute_dirac(self):
#
#todo: time limiting !
# http://efreedom.com/Question/1-366682/Limit-Execution-Time-Function-Call-Python
# google search: python limit execution time Windows
#
    def signal_handler(signum, frame):
      print 'Signal handler called with signal', signum
      print 'takes too long !'
      #raise Exception("Timed out!")
      #raise IOError("Take too long !")
      raise Error("Takes too long !")

    home_dir=os.getcwd()
    print 'changing directory from  '+home_dir+' to:'+self.wrkdir_obj.workdir
    # go to the working directory
    os.chdir(self.wrkdir_obj.workdir)

    f = open(home_dir+"/"+self.output_file, 'w',0)

    if (self.set_timeout_limit) :
      print 'signal stuff initiated, self.timeout=',self.timeout
      signal.signal(signal.SIGALRM, signal_handler)
      signal.alarm(self.timeout)   # 

    try:
      print 'executing dirac ...'
      p=subprocess.Popen('dirac.x', shell=False, stdin=subprocess.PIPE, stdout=f)
      print 'p.pid=',p.pid
      p.wait()
      f.close()
      signal.alarm(0)    
    except:
      print "dirac run interrupted !"
      p.kill()

    #os.system(command) > self.output_file

    print 'pam: end of executable run...'

    #subprocess.Popen(command,stdout=self.output_file)
    #self.output_file=subprocess.Popen(command)
    #self.output_file= subprocess.Popen([command]).communicate()[0]
    #f.write('xxx')
    #subprocess.Popen(command,stdout=f)
    #sys.stdout=f
    #f.write(os.system(command))

    #f_so.flush()

    #f_so.close()
    #f_se.close()
    
     #f_so.wait()
    # restore stdout
    #sys.stdout = save_stdout
    #sys.stderr = save_stderr

    #f_so.close()

    #sys.stdout.flush()
    #sys.stderr.flush()

    print 'end of executable run...'

    # do the profiling after the execution
    if self.profile: self.profile_dirac()

  def profile_dirac(self):
    print 'profiling...'
    #default output after finishing DIRAC
    GMON="gmon.out"
    #os.system(....)

  def remove_work_dir(self):
    print 'removing work-dir'
    self.wrkdir_obj.remove_workdir()

  def print_header(self):
    print '*****************************************'  
    print '**** OUTPUT FROM DIRAC PYTHON SCRIPT ****'  
    print '*****************************************'
    print '                                        '  
    print "   User@Host   : `whoami`@`hostname`"
    print "   Date        : `date`"
    print "   Command     : $0 $ARG_SAVE "
    print "   PID         : $$"
    print "   Input dir   : $LOC"
    print "   Work dir    :"

#def notify_user() :
# replacement
#  print 'notify_user...'

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class pam_run(pam_enviroment):
#-----------------------------------------------------
#-----------------------------------------------------
  def __init__(self):
    print 'pam_run class initialize'
    self.pam_variables=pam_enviroment()

  def perform(self):
    print 'performing own pam run..'
    self.pam_variables.process_pam_arguments()
    self.pam_variables.print_header()
    self.pam_variables.set_wrkdir()
    self.pam_variables.set_output_file()
    self.pam_variables.execute_dirac()
    self.pam_variables.remove_work_dir()


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
   dirac_run = pam_run()
   dirac_run.perform()

if __name__ == '__main__':
   sys.exit(main(*sys.argv))
