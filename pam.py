#! /usr/bin/env python

##################################################
#
# Python replacement of the pam script
#
# Reads crucial data from Makefile.config
#
# Written by Miro ILIAS,Banska Bystrica, October, November, 2010
#
##################################################

import os
import sys
import string 
import subprocess
import optparse
import distutils.util 
import shutil
import signal
import stat
import time


# installation directory inserted in by the configure script
global install_dir, install_dir_default
install_dir_default=''
#install_dir_default='/home/dv6-2145ec/Work/qch/software/dirac/gfortran_gcc_noopt_ownlibs_ilp32/adm_hsfsdyn_cap/trunk'
install_dir_default='/home/ilias/qch_work/qch_software/dirac/gfortran_noopt_ILP64_ownlibs/adm/trunk'
install_dir=install_dir_default

def append_dirac_py_modules(dirac_install_dir):
# append directory with common dirac python-modules and imports them
  modules_directory = dirac_install_dir + '/utils'
  if (not os.access(modules_directory,os.R_OK) ):
    print 'wrong directory or directory not readable !'
    sys.exit("wrong directory !")
  sys.path.append(modules_directory)
  # import common # does not work..
  #from common import *
  from common import is_windows
  print 'common dirac modules loaded !'
  print 'is_windows=',is_windows()

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class workdir_class:
###########################################################################################
#
#    class dealing with the work directory only
#
###########################################################################################
  global install_dir
  def __init__(self):
    print 'work_direcory class initialize'
    self.workdir=''

    #self.workpath='/tmp/ilias'
    #self.workpath='/tmp/'+pwd.getpwuid(os.getuid())[0]
    self.workpath='/tmp/'+os.getlogin()
    #TODO: Windows user login name !!!
    self.keepwrk=False

    self.cleanwrkdir=False

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
    print 
    print 'Content of the work directory',os.getcwd(),':'
    print '-------------------------------------------------------------------------------' 
    #print 'file............size'
    for filename in os.listdir(self.workdir):
      filesize_kB=float(os.path.getsize(self.workdir+"/"+filename))/1024
      # create a dictionary to hold file info
      file_name = filename; file_stats = os.stat(file_name)
      file_info = {
             'fname': file_name,
             'fsize': file_stats [stat.ST_SIZE],
             'f_lm': time.strftime("%m/%d/%Y %I:%M:%S %p",time.localtime(file_stats[stat.ST_MTIME])),
             'f_la': time.strftime("%m/%d/%Y %I:%M:%S %p",time.localtime(file_stats[stat.ST_ATIME])),
             'f_ct': time.strftime("%m/%d/%Y %I:%M:%S %p",time.localtime(file_stats[stat.ST_CTIME]))
      } 
      print "%20s"%filename,".....%10.4f"%filesize_kB,' kB'," last accessed = %(f_la)s" % file_info
      if (self.iprint > 7) :
        print
        print "file name = %(fname)s" % file_info
        print "file size = %(fsize)s bytes" % file_info
        print "last modified = %(f_lm)s" % file_info
        print "last accessed = %(f_la)s" % file_info
        print "creation time = %(f_ct)s" % file_info
        print
        if stat.S_ISDIR(file_stats[stat.ST_MODE]):
          print "This a directory"
        else:
          print "This is not a directory"
        print
        print "A closer look at the os.stat(%s) tuple:" % file_name
        print file_stats
        print
        print "The above tuple has the following sequence:"
        print """st_mode (protection bits), st_ino (inode number),
           st_dev (device), st_nlink (number of hard links),
           st_uid (user ID of owner), st_gid (group ID of owner),
           st_size (file size, bytes), st_atime (last access time, seconds since epoch),
           st_mtime (last modification time), st_ctime (time of creation, Windows)"""

  def remove_workdir(self):
    #print 'removes the work_direcory'
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

     # print level
    self.iprint=1
     ### contains installation directory ...
    self.dirac_home=install_dir

    ## basis set directory
    #self.basdir=self.dirac_home+"/basis:"+self.dirac_home+"/basis_dalton"
    self.basdir=''

    # ... object from the inherited "workdirectory" class 
    self.wrkdir_obj = workdir_class()
    # ... output file
    self.output_file=''

    # parallel stuff 
    self.parallel=False
    self.mpirun=''
    self.mpi=1

    self.dirac_cmd=''

    self.profile = False

    self.profiler = ''

    self.dirac_exe=self.dirac_home
 
    ### put/get files ###
    self.put_files=''
    self.get_files=''

    self.mem = 64
    self.maxmem=''
    self.slavemem = 64

    self.full_restart=False

    # debugging
    self.debug=False
    self.debugger=''

    # input directory
    self.input_directory='' 
    # input files
    self.input_molecule=''
    self.input_dirac=''
    # output name from *inp *mol names
    self.output_name=''
    # output suffix default
    self.output_suffix="out"
    # output file...
    self.output_file=''

    # do not backup files
    self.nobackup=False

    # default - do archivation of files
    self.no_archive=False
    # list of possible files for the archivation
    self.archive_files=['MOLECULE.MOL','MOLECULE.XYZ','DIRAC.INP','DFCOEF','*.cube','AOPROPER', \
                        'AOMOMAT','BSSMAT','matexport.m','DFCOEF.BEST']

    # all files from the archive are unpacked, it is assumed that archive file is copied into workdir 
    self.fullrestart=False
 
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
    usage="usage: %prog [options] [*.mol] [*.inp]"
    parser = OptionParser(usage)

  # ... get in the print level 
    this_label='script verbose print level, default:'+str(1)
    parser.add_option("-i", "--iprint", type="int", nargs=1,  
                  action="store", dest="dest_iprint", default=1,
                  help=this_label, metavar="INT_NUM")

  # ... archivation of files
    this_label='Do not create a .tgz archive for the calculation; default - create.'
    parser.add_option("--noarch",
                  action="store_true", dest="dest_noarch", 
                  help=this_label)

  # ... full-restart  ...
    this_label='All files from ARCHIVE.tgz are unpacked, not just DFCOEF. \
                This allows you to omit the molecule and input files. '

    parser.add_option("--full-restart",
                  action="store_true", dest="dest_fullrestart", 
                  help=this_label)

  # ... output suffix 
    this_label='file suffix for output; default:'+self.output_suffix
    parser.add_option("--run", type="string", nargs=1,  
                  action="store", dest="dest_run", 
                  help=this_label, metavar="STRING")

  # ... get in the "mw"
    this_label='memory (in megawords) for the master, default:'+str(self.mem)
    parser.add_option("--mw", type="int", nargs=1,  
                  action="store", dest="dest_mw", default=64,
                  help=this_label, metavar="INT_NUM")

    # ... get in the "nw" ...
    this_label='memory (in megawords) for nodes, default:'+str(self.slavemem)
    parser.add_option("--nw", type="int", nargs=1,  
                  action="store", dest="dest_nw", default=64,
                  help=this_label, metavar="INT_NUM")

    # ... get in the "aw"
    this_label='maximum memory(in megawords) that can be allocated dynamically (per allocation) '
    parser.add_option("--aw", type="int", nargs=1,  
                  action="store", dest="dest_aw",
                  help=this_label, metavar="INT_NUM")

    # ... get in the parallel stuff, TODO: add default from Makefile.config
    this_label="mpi launcher, default..."
    parser.add_option("--mpirun", type="string", nargs=1,
                  action="store",dest="dest_mpirun",  
                  help=this_label,metavar="STRING")

    this_label='run parallel with n nodes, default: serial run'
    parser.add_option("--mpi", type="int", nargs=1,  
                  action="store", dest="dest_mpi", 
                  help=this_label, metavar="INT_NUM")

    # ... get in the timeout limit (in seconds)
    this_label='limit dirac.x execution time (in seconds)'
    parser.add_option("--timeout", type="int", nargs=1,  
                  action="store", dest="dest_timeout",
                  help=this_label, metavar="INT_NUM")

    # .. read in the install_dir, what is the basic global variable...
    this_label="pam.py\'s directory, default:"+install_dir_default
    parser.add_option("--install_dir", type="string", nargs=1,
                  action="store",dest="dest_install_dir",  
                  help=this_label,metavar="STRING")

    # .. read in the install_dir, what is the basic global variable...
    this_label="basis set directory, default:"+self.basdir
    parser.add_option("--basdir", type="string", nargs=1,
                  action="store",dest="dest_basdir",  
                  help=this_label,metavar="STRING")

    # ... no backup 
    this_label="do not backup files (default: backup them)"
    parser.add_option("--nobackup", 
                  action="store_true",dest="dest_nobackup",  
                  help=this_label)

    # .. get files from the work dir
    this_label="get files from the working directory"
    parser.add_option("--get",type="string", nargs=1,
                  action="store",dest="dest_get",
                  help=this_label,metavar="STRING")

    # .. copy files into the work dir
    this_label="copy files into the working directory"
    parser.add_option("--put",type="string", nargs=1,
                  action="store",dest="dest_put",
                  help=this_label,metavar="STRING")

    # ... clean workdir before calculation
    this_label="clean the directory before calculations (default: do not clean)"
    parser.add_option("--clean", 
                  action="store_true",dest="dest_clean",  
                  help=this_label)

    # ... keep workdir after calculations
    this_label="keep the work directory after calculations (default: remove it)"
    parser.add_option("--keepwrk", 
                  action="store_true",dest="dest_keepwrk",  
                  help=this_label)

    # ... keep workdir after calculations
    this_label="run in debugger (default: do not run)"
    parser.add_option("--debug", 
                  action="store_true",dest="dest_debug", 
                  help=this_label)

    this_label="debugger (full path), default: "+self.debugger
    parser.add_option("--debugger",
                  action="store",dest="dest_debugger",
                  help=this_label, metavar="STRING")

    # process all
    (options, args) = parser.parse_args()
    
    # now assign individual arguments
    if (options.dest_iprint is not None) :
     self.iprint=int(options.dest_iprint)
     if self.iprint > 7:
      print 'assigned iprint=',self.iprint

    # now assign individual arguments
    if (options.dest_run is not None) :
     self.output_suffix=options.dest_run
     if self.iprint > 7:
      print 'assigned suffix=',self.output_suffix


    if (options.dest_mw is not None) :
     self.mem=int(options.dest_mw)
     if self.iprint > 7:
      print 'assigned mw=',self.mem

    if (options.dest_nw is not None) :
     self.memslave=int(options.dest_nw)
     if self.iprint > 7:
      print 'assigned mw=',self.memslave

    if (options.dest_aw is not None) :
     self.maxmem=int(options.dest_aw)
     if self.iprint > 7:
      print 'assigned maxmem=',self.maxmem

    if (options.dest_mpi is not None) :
     self.mpi=int(options.dest_mpi)
     # todo: check number of mpi's wrt to computer architecture
     self.parallel=True
     if self.iprint > 7:
      print 'parallel TRUE !',options.dest_mpi
      print 'assigned mpi=',self.mpi

    if (options.dest_mpirun is not None):
     self.mpirun=options.dest_mpirun
     # todo: check it
     self.parallel = True 
     if self.iprint > 7:
      print 'parallel TRUE  2!'
      print 'assigned mpirun=',self.mpirun

    if (options.dest_install_dir is not None) :
     install_dir=str(options.dest_install_dir)
     if self.iprint > 7:
      print 'assigned install-dir=',install_dir
    elif (install_dir_default is None) :
     sys.exit("you must set the installation directory !")

    if (options.dest_basdir is not None) :
     self.basdir=str(options.dest_basdir)
     if self.iprint > 7:
      print 'assigned basdir=',self.basdir
    else:
     #self.basdir=self.dirac_home+"/basis:"+self.dirac_home+"/basis_dalton"
     self.basdir=install_dir+"/basis:"+install_dir+"/basis_dalton"
     if self.iprint > 7:
      print 'default basdir=',self.basdir

    if (options.dest_timeout is not None) :
      self.timeout=int(options.dest_timeout)
      self.set_timeout_limit=True
      if self.iprint > 7:
        print 'assigned timeout=',self.timeout

    if (options.dest_nobackup is not None):
      self.nobackup=options.dest_nobackup      
      if self.iprint > 7:
        print 'assigned nobackup=',self.nobackup

    if (options.dest_noarch is not None):
      self.no_archive=options.dest_noarch
      if self.iprint > 7:
        print 'assigned no_archive=',self.no_archive

    if (options.dest_put is not None):
      self.put_files=options.dest_put
      if self.iprint > 7:
        print 'assigned put_files=',self.put_files

    if (options.dest_get is not None):
      self.get_files=options.dest_get
      if self.iprint > 7:
        print 'assigned get_files=',self.get_files

    if (options.dest_fullrestart is not None):
      self.fullrestart=dest_fullrestart
      if self.iprint > 7:
        print 'assigned full-restart=',self.fullrestart

    if (options.dest_clean is not None) :
      self.wrkdir_obj.cleanwrkdir=options.dest_clean
      if self.iprint > 7:
        print 'assigned wrkdir_obj.cleanwrkdir==',self.wrkdir_obj.cleanwrkdir

    if (options.dest_keepwrk is not None) :
      self.wrkdir_obj.keepwrk=options.dest_keepwrk
      if self.iprint > 7:
        print 'assigned wrkdir_obj.keepwrk =',self.wrkdir_obj.keepwrk

    if (options.dest_debug is not None) :
      self.debug=True
      if self.iprint > 7:
        print 'do debugging,',self.debug

    if (options.dest_debugger is not None) :
      self.debugger=options.dest_debugger
      if self.iprint > 7:
        print 'assigned debugger,',self.debugger

    #################################################################
    # VERY IMPORTANT: append external dirac modules to the path dir
    #################################################################
    append_dirac_py_modules(install_dir)

    # assign input directory ...
    self.input_directory=os.getcwd()

    # remaining arguments - *.mol and *.inp files
    process_remaining_args(args)

    if (self.iprint > 7) :
      print 'Parallel ?:',self.parallel

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

  def copy_files(self):
  # copies files from home to the workdir
    if (self.put_files is not None):
        for file_name in self.put_files:
          #full_file_name=self.wrkdir_obj.workdir+"/"+file_name
          full_file_name=self.input_directory+"/"+file_name
          if ( os.path.isfile(full_file_name) ):
            file_name_wrkdir=self.wrkdir_obj.workdir+"/"+file_name
            self.wrkdir_obj.copy_file(full_file_name,file_name_wrkdir)


  def set_output_file(self):
  # sets the name of the output file by using "self.output_name"
    #print 'setting output file..'
    self.output_file=self.output_name+"."+self.output_suffix
    if (self.iprint > 7):
      print 'the output file..',self.output_file
 
  def set_dirac_environment(self):
  # export important variables for DIRAC
  # TODO: implement for WINDOWS

    """ basis set directory """
    os.putenv("BASDIR",self.basdir)


    dirwrk  =  str( self.mem     *1000000 )
    if (self.parallel) :
      os.putenv("DIRPAR","1")
      dirnod  =  str( self.memslave*1000000 )

    if (self.maxmem is not None): 
      dirmaxmem=str( self.maxmem*1000000 )
      os.putenv("DIRMAX",dirmaxmem)
 
    #print dirwrk,dirwrk
    os.putenv("DIRWRK", dirwrk)
    if (self.parallel) :
      os.putenv("DIRNOD", dirnod)
 
  # Next line is for IBM only but will not hurt other systems
    os.putenv("XLFRTEOPTS","'namelist=old'")

    # control printout ... giving None
    #print "DIRWRK=",os.getenv("$DIRWRK")

  def copy_files(self):
  # before running DIRAC copy files into working directory
    print 'copying files into working directory...'
    if (self.fullrestart):
      print 'restarting from copied ARCHIVE file '

  def set_dirac_cmd(self):
  # set the dirac command
    if (self.parallel):
      # DIRAC_CMD="$DIRAC_CMD -np $MPI_NPROC -wd $WRK"
      # and/if  $DIRAC_CMD -machinefile $MACH_FILE $DIRAC_CMD
      self.dirac_cmd=self.mpirun+" -np "+str(self.mpi)+" -wd "+self.wrkdir_obj.workdir
      self.dirac_cmd=self.dirac_cmd+' dirac.x'
    else:
      if (self.debug):
       #TODO: adapt it for Windows
        self.dirac_cmd=self.debugger+' '+'dirac.x | tee -a '+os.getcwd()+"/"+self.output_file
        #self.dirac_cmd=self.debugger+' '+'dirac.x'
      else:
        self.dirac_cmd='dirac.x'

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def execute_dirac(self):
  #############################################################################################
  #
  # Excecute dirac.x, either in serial or parallel mode
  #
  #############################################################################################
    def signal_handler(signum, frame):
    # ... controling the time of the execution
      print 'Signal handler called with signal', signum
      print 'takes too long !'
      #raise Exception("Timed out!")
      #raise IOError("Take too long !")
      raise Error("Takes too long !")

    print 'changing directory from Input dir to Work dir'
    # go to the working directory
    os.chdir(self.wrkdir_obj.workdir)

    if (self.nobackup is False):
      from common import backup_file
      backup_file(some_file=self.input_directory+"/"+self.output_file,upto=None)

    if (self.set_timeout_limit) :
      print 'signal stuff initiated, self.timeout=',self.timeout
      signal.signal(signal.SIGALRM, signal_handler)
      signal.alarm(self.timeout)  # 

    if (self.parallel) :
      print 'executing parallel dirac'
    else :
      print 'executing serial dirac'

    try:
      if (not self.debug):
        print 'executing dirac ...',self.dirac_cmd
        output_file_full=self.input_directory+"/"+self.output_file
        f=open(output_file_full,'w',0)
        #p=subprocess.Popen(self.dirac_cmd, shell=False, stdin=subprocess.PIPE, stdout=f)
        p=subprocess.Popen(self.dirac_cmd,shell=True,stdin=subprocess.PIPE,stdout=f)
        #print 'p.pid=',p.pid
        p.wait()
        f.close()
      else:
        print 'debugging dirac...',self.dirac_cmd
        os.system(self.dirac_cmd) # this works...
        #p=subprocess.Popen(self.dirac_cmd, shell=True, stdin=subprocess.PIPE, stdout=f)

      signal.alarm(0)    
    except:
      print "dirac launching failed with command ",self.dirac_cmd

    print 'pam: end of executable run...'

    # do the profiling after the execution
    if self.profile: self.profile_dirac()

  def profile_dirac(self):
    print 'profiling...'
    #default output after finishing DIRAC
    GMON="gmon.out"
    #os.system(....)

  def save_files(self):
    # ... save files in various manners from the work directory 
    from common import backup_file
    if (not self.no_archive):
      # set name of the archive file
      archive_name=self.input_molecule+"_"+self.input_dirac+".tgz"
      print 'creating archive file ',archive_name
      archive_full_name=self.input_directory+"/"+archive_name
      print 'full archive file ',archive_full_name

      # select existing files in self.archive_files for the archive
      existing_archive_files=''
      for archfile in self.archive_files:
        if (os.path.isfile(archfile)):
           existing_archive_files=existing_archive_files+' '+archfile
      print 'existing archive files:',existing_archive_files
      # backup 
      if (self.nobackup is False):
         backup_file(some_file=archive_full_name,upto=None)
      #archive_command='tar c '+self.archive_files+' | gzip > '+archive_full_name
      archive_command='tar c '+existing_archive_files+' | gzip > '+archive_full_name
      # own packing (Linux)
      os.system(archive_command)

      # take care of the rest - get files from workdir to home dir
      if (self.get_files is not None):
        for file_name in self.get_files:
          full_file_name=self.wrkdir_obj.workdir+"/"+file_name
          if ( os.path.isfile(full_file_name) ):
            file_name_home=self.input_directory+"/"+file_name
            self.wrkdir_obj.copy_file(full_file_name,file_name_home)


  def remove_work_dir(self):
    print 'removing work-dir'
    self.wrkdir_obj.remove_workdir()

  def print_header(self):
    import datetime
    from common import get_username,get_hostname
    print '*****************************************'  
    print '**** OUTPUT FROM DIRAC PYTHON SCRIPT ****'  
    print '*****************************************'
    print '                                        '  
    print "   User@Host   : ",get_username(),"@",get_hostname()
    print "   Date        : ",datetime.datetime.now()
    print "   Command     : $0 $ARG_SAVE "
    print "   PID         :",os.getpid()
    print "   Input dir   :",self.input_directory
    print "   Work dir    :",self.wrkdir_obj.workdir

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
    self.pam_variables.set_wrkdir()
    self.pam_variables.copy_files()
    self.pam_variables.set_output_file()
    self.pam_variables.print_header()
    self.pam_variables.set_dirac_environment()
    self.pam_variables.set_dirac_cmd()
    self.pam_variables.execute_dirac()
    self.pam_variables.save_files()
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
