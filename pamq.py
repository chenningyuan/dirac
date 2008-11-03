#! /usr/bin/env python
#
#  script for submitting DIRAC pam job into a queue system
#
# writen by Miro ILIAS,2008 after Luuk's pamq script
#
# further improvement still in progress
#

import  sys, os, string, optparse

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class PamQVariables:

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def __init__(self):
  #
  #  initialize set of variables ...
  #

    # print level
    self.iprint=1

    self.pam_arguments=''
    self.submit_job=True

    # ... specify the default job name
    self.job_name='queue_job.'+str(os.getpid())

    self.job_inp='input_test'

    self.job_out='queue_job.'+str(os.getpid())+'.log_out'

    #self.queue_system='pbs'
    self.queue_system=''

   #  self.queue_name='default'
   # self.queue_name='workq'
   # self.queue_name='express'
    self.queue_name=''
    self.queue_name_specified=False

    # time 
    self.queue_time='0:20'
    self.queue_time_specified=False

    # ... is pure serial job (no mpi added to pam)
    self.is_serial=False

   # self.nnodes=1
    self.nnodes=0

   # self.proc_per_node=4
   # self.proc_per_node=1
    self.proc_per_node=0

   # self.pam_script_name="pam"
   # self.pam_script_dir="~/bin"
   # self.pam_script=self.pam_script_dir + "/" + self.pam_script_name
    self.pam_script=''

   # self.queue_script="pamq."+str(os.getpid())+".batch_script"
    self.queue_script="queue_job."+str(os.getpid())+".batch_script"

    # home directory ...
    self.home_dir=os.getenv("PWD")

    # ... configuration file with selected defaults ...
    self.prepare_config_file=False

    #-------------------------------------------------------------------
    #   ... list of implemented queue system for the script:
    #-------------------------------------------------------------------
    self.implemented_queue_systems=['pbs','bsub','loadleveler']

    # ... set up file with defaults
    self.pamq_defaults_file=str(os.getenv("HOME"))+'/PAMQ_DEFAULTS'

    self.read_from_config_file=False
   # if os.path.isfile("PAMQ_DEFAULTS"):
    if os.path.isfile(self.pamq_defaults_file):
      self.read_configure_file()
      self.read_from_config_file=True
     # print 'defaults read from file'

    if self.iprint >8 :
      print ' setting home directory:',self.home_dir
      print 'INITIALIZATION of class object done'

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def set_pam_mpi(self):
  #
  # class method: set up mpi parameter for pam script
  #
    if GlobalVars.iprint > 7:
      print '\nsetting mpi for pam' 
      print 'set_pam_mpi: extracted pam arguments:', GlobalVars.pam_arguments
    # print GlobalVars.pam_arguments

    # extend pam flags with the necessary -mpi arguments
    mpi_procs=GlobalVars.nnodes*GlobalVars.proc_per_node

    if GlobalVars.iprint >= 5:
      print 'set_pam_mpi: calculated # mpi=',mpi_procs
      print 'using GlobalVars.nnodes=',GlobalVars.nnodes
      print 'and GlobalVars.proc_per_node=',GlobalVars.proc_per_node

    GlobalVars.pam_arguments=' -mpi '+str(mpi_procs)+GlobalVars.pam_arguments

    if GlobalVars.iprint >= 7:
      print 'set_pam_mpi: mpi-extended pam arguments:', GlobalVars.pam_arguments

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def create_configure_file(self):
  #
  #      class method: create configure file with important parameters as DEFAULTS
  #
    if GlobalVars.iprint >= 12:
      print 'preparing PAMQ_DEFAULTS file'
    
    error_exit=False
 
   # f = open('PAMQ_DEFAULTS', 'w') 
   # f = open('~/PAMQ_DEFAULTS', 'w') 
   # f = open(GlobalVars.pamq_defaults_file, 'w') 
    f = open(self.pamq_defaults_file, 'w') 

   # ... check if the file is open ...
   # if not file_is_open(f) :
   #   print '\n File ',self.pamq_defaults_file,' is not open ! Error ! \n' 
   #   sys.exit()

    print 'opening file for saving defaults in ~/ dir:',self.pamq_defaults_file

    f.write('###  \n')
    f.write('###  SOME IMPORTANT DEFAULT VARIABLES FOR THE "pamq.py" SCRIPT ###\n')
    f.write('###  \n')

   # print 'len(GlobalVars.queue_system)=',len(GlobalVars.queue_system)

    if (len(GlobalVars.queue_system)>0) :
      f.write('# ... default queue system \n')
      f.write('QUEUE_SYSTEM='+GlobalVars.queue_system+'\n')
      print 'saved QUEUE_SYSTEM='+GlobalVars.queue_system
    else :
      print '\ncreate_configure_file: YOU MUST SPECIFY QUEUE SYSTEM!\n'
      error_exit=True 

    if (len(GlobalVars.queue_name)>0) :
      f.write('# ... default queue name \n')
      f.write('QUEUE_NAME='+GlobalVars.queue_name+'\n')
      print 'saved QUEUE_NAME=',GlobalVars.queue_name

    if (GlobalVars.nnodes>0) :
      f.write('# ... default number of nodes for (parallel) run\n')
      f.write('NODES='+str(GlobalVars.nnodes)+'\n')
      print 'saved NODES=',GlobalVars.nnodes
    else :
      print '\ncreate_configure_file: YOU MUST SPECIFY DEFAULT NUMBER OF NODES !\n'
      error_exit=True 

    if (GlobalVars.proc_per_node>0) :
      f.write('# ... default number of processors (CPU) per 1 node\n')
      f.write('PROCESSORS_PER_NODE='+str(GlobalVars.proc_per_node)+'\n')
      print 'saved PROCESSORS_PER_NODE=',GlobalVars.proc_per_node
    else :
      print '\ncreate_configure_file: YOU MUST SPECIFY DEFAULT NUMBER OF PROCESSORS_PER_NODE !\n'
      error_exit=True 

    if (len(GlobalVars.pam_script)>0) :
      f.write('# ... default location of the DIRAC pam script \n')
      f.write('PAM_SCRIPT='+GlobalVars.pam_script+'\n')
      print 'saved PAM_SCRIPT=',GlobalVars.pam_script
    else :
      print '\ncreate_configure_file: YOU MUST SPECIFY LOCATION OF DIRAC PAM SCRIPT!\n'
      error_exit=True 

    if (GlobalVars.is_serial) :
      f.write('# .... default serial jobs !\n')
      f.write('SERIAL_JOBS_ONLY\n')
      print 'saved SERIAL_JOBS_ONLY'

    f.close()

    if (error_exit) :
      print('script error exit !')
      sys.exit()

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def read_configure_file(self):
  #
  # class method: read parameters (DEFAULTS) from the configure file
  #
    if self.iprint >= 12 :
      print 'in ... read_configure_file'

    #file = open('PAMQ_DEFAULTS','r')
   # file = open('PAMQ_DEFAULTS','r')
    file = open(self.pamq_defaults_file,'r')

    # read line after line
    for line in file:
      #print line
      if (line.find('QUEUE_SYSTEM') >= 0 and not line.find('#') >= 0):
         queue_system=line[line.find('=')+1:len(line)]
         #print 'read QUEUE_SYSTEM=',queue_system
         self.save_queue_system(queue_system) 
    
      if (line.find('QUEUE_NAME') >= 0 and not line.find('#') >= 0):
         queue_name=line[ line.find('=')+1:len(line)]
        # print 'read QUEUE_NAME=',queue_name
         self.queue_name=queue_name
         self.queue_name_specified=True

      if (line.find('NODES') >= 0 and not line.find('#') >= 0):
         nodes=line[ line.find('=')+1:len(line)]
         #print 'read NODES=',int(nodes)
         self.nnodes=nodes

      if (line.find('PROCESSORS_PER_NODE') >= 0 and not line.find('#') >= 0):
         ppn=line[line.find('=')+1:len(line)]
         #print 'read PPN=',int(ppn)
         self.proc_per_node=ppn

      if (line.find('PAM_SCRIPT') >= 0 and not line.find('#') >= 0):
         pam_script=line[line.find('=')+1:len(line)]
        # print 'read PAM_SCRIPT dir=',pam_script
         self.save_pam_script(pam_script,True)

      if (line.find('SERIAL_JOBS_ONLY') >= 0 and not line.find('#') >= 0):
        # print 'serial runs as default !'
         self.is_serial=True

    file.close()

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def save_queue_system(self,qsys_in):
  #
  #  class method: check and save queue_system into proper variable
  #
  #  we are inside class method definition, so we are to use 'self.' specification
  #

 #   print ' hello world '.replace(' ', '')
  #  z=' hello world '.replace(' ', '')
  #  print 'z=',z

   # print 'qsys_in=',qsys_in
   # print 'len(qsys_in)=',len(qsys_in)
   # print 'qsys_in.isalnum=',qsys_in.isalnum()

    # .. check if there is empty space in the string
   # qsys_in.strip() 

   # qsys_in.replace('\n','')
   # print 'after removing extra EOL character, qsys_in=',qsys_in
   # print 'len(qsys_in)=',len(qsys_in)
   # print 'qsys_in.isalnum=',qsys_in.isalnum()

    # check the ending character !
    if qsys_in.endswith('\n'):
     # print 'ending with EOL !' 
      qsys_in=qsys_in[:-1]

 #   print 'after removing extra EOL character, qsys_in=',qsys_in
 #   print 'len(qsys_in)=',len(qsys_in)
 #   print 'qsys_in.isalnum=',qsys_in.isalnum()

    # check if the queue system is in the list ...
    if (qsys_in not in self.implemented_queue_systems) :
      print '\nsave_queue_system: entering queue system, qsys_in=',qsys_in
      print 'self.implemented_queue_systems=',self.implemented_queue_systems
      sys.exit()
    else :
      #print 'save_queue_system: entering queue system VERYFIED!'
      self.queue_system=qsys_in

    if self.iprint >= 12:
      print 'save_queue_system: saved class member queue_system=',self.queue_system

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def save_pam_script(self,pam_script_in,from_file=False):
  #
  #  class method: check and save queue_system into proper variable
  #
  #  we are inside class method definition, so we are to use 'self.' specification
  #

    # check the ending character (if reading from file)!
    if pam_script_in.endswith('\n'):
      #print 'pam_script ending with EOL !' 
      pam_script_in=pam_script_in[:-1]

    # ... verify the existence of the pam script !
    if os.path.isfile(pam_script_in) :
      #print 'pam script is present!'
      self.pam_script=pam_script_in
    else:
      print 'save_pam_script: entering pam_script_in=',pam_script_in
      print 'save_pam_script: pam script does not exist! Error exit !'
      if (from_file) :
        print 'Bad pam_script value in the file PAMQ_DEFAULTS !' 
      sys.exit()

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_pbs_script():
#
#
#   prepare the file for script and submit
#
#
  global GlobalVars

  if GlobalVars.iprint >8 :
     print '\nin prepare_pbs_script...'

  # unite into one local string
  GlobalVars.pam_script_args=GlobalVars.pam_script+GlobalVars.pam_arguments

  if GlobalVars.iprint >8 :
    print ' home directory:',GlobalVars.home_dir
    print ' GlobalVars.pam_script_args=',GlobalVars.pam_script_args

  f=open(GlobalVars.queue_script,'w') # open file with queue commands

  s = ''
  s += '''\
#!/bin/bash
#
#         Queue controls are directly readable from the submitting batch script
#----------------------------------------------------------------------------------------
#
#
#     job name
#
#PBS -N %s
  ''' % GlobalVars.job_name
  s += '''\
#
# number nodes (NOTE: different from number of processors specified by -mpi pam flag !)
#
#PBS -l nodes=%i
''' % GlobalVars.nnodes

#
#   if the queue_name is specified in the input, add it
#
  if GlobalVars.queue_name_specified:
    s += '''\
#
# specify queue name
#
#PBS -q %s
    ''' % GlobalVars.queue_name
  else:
    s += '''\
#
#
# queue name not specified (command incative), using queue system default 
#
##PBS -q %s
    ''' % GlobalVars.queue_name
#
#  if time is specified in the input, add it
#
  if GlobalVars.queue_time_specified:
    s += '''\
#
#
# specify time
#
#PBS -l walltime=%s
    ''' % GlobalVars.queue_time
  else:
    s += '''\
#
# queue time not specified, using system default and inactive
#
##PBS -l walltime=%s
    '''%GlobalVars.queue_time
# ...continue ...
  s += '''\
#
# unite stdout with stderr
#
#PBS -j oe
#
#   queue job standard output - stdout
#
#PBS -o %s
''' % GlobalVars.job_out
  s += '''\
#
#
# ... print out important PBS variables into standard output
#
echo '------------------------------------'
echo '  *** PBS control variables ***  ' 
echo '------------------------------------'
echo '  $ PBS_O_HOME='$PBS_O_HOME
echo '   $PBS_O_HOST='$PBS_O_HOST
echo '   $PBS_SERVER='$PBS_SERVER
echo '  $PBS_O_QUEUE='$PBS_O_QUEUE
echo '$PBS_O_WORKDIR='$PBS_O_WORKDIR
echo '  $PBS_JOBNAME='$PBS_JOBNAME
echo '    $PBS_JOBID='$PBS_JOBID
echo ' $PBS_NODEFILE='$PBS_NODEFILE
echo '------------------------------------'
#
#
echo ''
echo '****  assigned nodes (content of the $PBS_NODEFILE): ****'
cat $PBS_NODEFILE
echo ''
#
# ... always change to the home directory 
#  from where you execute DIRAC pam  
#
cd %s
''' % GlobalVars.home_dir
  s += '''\
#
#   ... execute own DIRAC pam script
#
%s 
''' % GlobalVars.pam_script_args
  s += '''\
#
#
# ... that's whole script, exit
exit
'''
# ... 

# ... write all string to the file
  f.write(s)
  f.close()

 # submit job ...
  if GlobalVars.submit_job:
    if GlobalVars.iprint >= 0 :
      print '\ngoing to submit job with prepared script:',GlobalVars.queue_script
      print '            standard output file will be:',GlobalVars.job_out
      print '                  queue job name will be:',GlobalVars.job_name,'\n'

    os.system("qsub "+GlobalVars.queue_script) # submitting !
    # ... some info ...
    os.system("echo '  -  Current user-status on queueing-system : '")
    os.system("echo 'Job id           Name             User             Time Use S Queue'")
    os.system("echo '---------------- ---------------- ---------------- -------- - -----'")
    os.system("qstat | grep $USER")
  else:
    print '\n not submitting,prepared only the script for the batch system:',GlobalVars.queue_script
    print ''


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_loadleveler_script():
#
# TODO : add !
#
  global GlobalVars

  if GlobalVars.iprint > 7:
    print 'prepare_loadleveler_script...'

  # unite into one local string
  GlobalVars.pam_script_args=GlobalVars.pam_script+GlobalVars.pam_arguments

  if GlobalVars.iprint >8 :
    print ' home directory:',GlobalVars.home_dir
    print ' GlobalVars.pam_script_args=',GlobalVars.pam_script_args

  f=open(GlobalVars.queue_script,'w') # open file with queue commands

  s = ''
  s += '''\
#
# loadl stuff
#
# @ shell            = /bin/sh
# @ network.MPI      = sn_all,not_shared,USnetwork.MPI   
# @ notification     = never
# @ queue      
# 
# 
# @ initialdir       = %s
''' % GlobalVars.home_dir
  s += '''\
#
# @ output           = %s
''' % GlobalVars.job_out
  s += '''\
#
# @ error            = %s
''' % GlobalVars.job_out
  s += '''\
#
# @ node             = %s
''' % str(GlobalVars.nnodes)
  s += '''\
# @ tasks_per_node   = %s
''' % str(GlobalVars.proc_per_node)
  s += '''\
#
# @ wall_clock_limit = %s
''' % GlobalVars.queue_time
  s += '''\
#
# @ environment      = LLSHELL=%s ;ENVIRONMENT=BATCH;
''' % os.getenv("SHELL")
  s += '''\
#
#  change to original directory
#
cd %s
''' % GlobalVars.home_dir
  s += '''\
#
export XLFRTEOPTS="namelist=old"
#
''' 

# ... write all string to the file
  f.write(s)
  f.close()

 # submit job ...
  if GlobalVars.submit_job:
    if GlobalVars.iprint >= 0 :
      print '\ngoing to submit job with prepared script:',GlobalVars.queue_script
      print '            standard output file will be:',GlobalVars.job_out
      print '                  queue job name will be:',GlobalVars.job_name,'\n'

    os.system("llsubmit "+GlobalVars.queue_script) # submitting !
    # ... some info ...
 #   os.system("echo '  -  Current user-status on queueing-system : '")
 #   os.system("echo 'Job id           Name             User             Time Use S Queue'")
 #   os.system("echo '---------------- ---------------- ---------------- -------- - -----'")
 #   os.system("qstat | grep $USER")
  else:
    print '\n not submitting,prepared only the script for the batch system:',GlobalVars.queue_script
    print ''


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_bsub_script():
#
# TODO : add !
#
  if GlobalVars.iprint > 7:
    print 'preparing bsub script...'

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_queue_script():
#
# prepare the script according to queue syustem and submit it
#
  global GlobalVars

  if GlobalVars.iprint > 7:
    print 'Going to prepare script for batch system...'

  if GlobalVars.queue_system=='pbs':
    process_pbs_script()	
  elif GlobalVars.queue_system=='loadleveler': 
    process_loadleveler_script()	
  elif GlobalVars.queue_system=='bsub': 
    process_bsub_script()	
  elif GlobalVars.queue_system=='xgrid': 
    print
#	process_xgrid_script()	
  else:
   print '\nprocess_queue_script: NOT defined queue_system, GlobalVars.queue_system=',GlobalVars.queue_system
   print 'process_queue_script: list of queue systems:\n',GlobalVars.implemented_queue_systems
   sys.exit()


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_pamq_arguments():
#----------------------------------------------------------------------------------------
#
#   process script arguments 
#
#----------------------------------------------------------------------------------------
  global GlobalVars
  from optparse import OptionParser

  if GlobalVars.read_from_config_file :
    usage = "usage: %prog arg1 arg2...; defaults read from file:"+GlobalVars.pamq_defaults_file
  else :
    usage = "usage: %prog arg1 arg2...; no defaults get from config file"
    
  parser = OptionParser(usage)

  #  fill GlobalVars.pam_script
  pam_script_add="location of pam script, default:"+GlobalVars.pam_script
  parser.add_option("-p", "--pam", type="string", nargs=1,
                  dest="dest_pam",default=GlobalVars.pam_script,
                  help=pam_script_add, metavar="STRING")

  #  fill GlobalVars.pam_arguments
  dest_pam_arguments=" pam arguments (like '"' -noarch [*.mol] [*.inp]'"')"
  parser.add_option("--pam_args", type="string", dest="dest_pam_args",
                  help=dest_pam_arguments, metavar="'"'arg1 arg2 ... '"' ")

  # ... get in the print level 
  iprint_label='script verbose print level, default:'+str(GlobalVars.iprint)
  parser.add_option("-i", "--iprint", type="int", nargs=1,  
                  action="store", dest="dest_iprint", default=GlobalVars.iprint, 
                  help=iprint_label, metavar="INT_NUM")

  # ... do not submit flag
  bool_submit=str(GlobalVars.submit_job)
  nosub_label='do not submit job, prepare only the queue script for batch, default-submit:'+bool_submit
  parser.add_option("--nosub", action="store_false", default=GlobalVars.submit_job,
                     dest="dest_nosub", help=nosub_label) 

  # ... specify as serial (no mpi)
  bool_serial=str(GlobalVars.is_serial)
  serial_label='submit job as serial (no mpi), default:'+bool_serial
  parser.add_option("-s","--serial", action="store_true",default=GlobalVars.is_serial,
                     dest="dest_is_serial", help=serial_label) 

  # ... get number of nodes
  if (int(GlobalVars.nnodes) <=0) :
    nodes_label='number of reserved nodes, default:'+str(GlobalVars.nnodes)+'=> YOU MUST SET IT!'
  else :
    nodes_label='number of reserved nodes, default:'+str(GlobalVars.nnodes)
  parser.add_option("-n", "--nodes", type="int",  nargs=1,
                  action="store",dest="dest_nnodes", default=GlobalVars.nnodes, 
                  help=nodes_label, metavar="INT_NUM")

  # ... GlobalVars.queue_system
  if (len(GlobalVars.queue_system)>0) :
    queue_system_label='specify the queue system, default :'+"'"+GlobalVars.queue_system+"'"
  else :
    queue_system_label='determine the queue system, no default=> YOU MUST SPECIFY IT !'
  parser.add_option( "--qsys", type="string", nargs=1,
                  action="store",dest="dest_queue_system", default=GlobalVars.queue_system, 
                  help=queue_system_label,metavar="STRING")

  # ... GlobalVars.job_name
  queue_job_name_label='specify the name for the queue run, default:'+"'"+GlobalVars.job_name+"'"
  parser.add_option( "--name", type="string", nargs=1,
                  action="store",dest="dest_queue_job_name",  
                  help=queue_job_name_label,metavar="STRING")

  # ... GlobalVars.queue_name
  if (GlobalVars.queue_name_specified) :
    queue_label='specify the queue, default :'+GlobalVars.queue_name
    parser.add_option("-q", "--queue", type="string", nargs=1,
                  action="store", dest="dest_queue_name",default=GlobalVars.queue_name,
                  help=queue_label,metavar="STRING")
  else :
    queue_label='specify the queue (no default)'
    parser.add_option("-q", "--queue", type="string", nargs=1,
                  action="store", dest="dest_queue_name", 
                  help=queue_label,metavar="STRING")

  # ... GlobalVars.proc_per_node ... "--ppn" not accepting !
  if (int(GlobalVars.proc_per_node) <=0):
    queue_ppn_label='specify number of processor per node, default:'+str(GlobalVars.proc_per_node)+'<= YOU MUST SET IT !'
  else:
    queue_ppn_label='specify number of processor per node, default:'+str(GlobalVars.proc_per_node)

  parser.add_option("--pnode", type="int", nargs=1, 
                   action="store", dest="dest_ppn", default=GlobalVars.proc_per_node,
                 #action="store", dest="dest_ppn", 
                  help=queue_ppn_label, metavar="INT")

  # ... GlobalVars.queue_time
 # queue_time_label='specify the time in the queue, default:'+GlobalVars.queue_time
  queue_time_label='specify the time in the queue (no default)'
  parser.add_option("-t", "--time", type="string", nargs=1,
                 # action="store", dest="dest_queue_time", default=GlobalVars.queue_time, \
                  action="store", dest="dest_queue_time", 
                  help=queue_time_label,metavar="STRING")

  # ... prepare file with selected default variables (configuration file) ...
  config_file_label='save read parameters (like qsys, queue... ) as defaults into PAMQ_DEFAULTS file,default - do not save'
  parser.add_option("--config", action="store_true", default=GlobalVars.prepare_config_file,
                     dest="dest_prep_config", help=config_file_label) 

  ###################################################################################################
  #
  #                                    ... finally assign
  #
  ###################################################################################################
  (options, args) = parser.parse_args()

# .. now assign individual arguments 

  if GlobalVars.iprint > 5:
    print '\n processing arguments:  len(sys.argv)=', len(sys.argv) 

  # ... handle queue job name
  if (options.dest_queue_job_name is not None) :
    GlobalVars.job_name=options.dest_queue_job_name
    #print '* assigned GlobalVars.job_name=',GlobalVars.job_name
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.job_name=',GlobalVars.job_name

  # ... handle iprint
  if (options.dest_iprint is not None) :
    GlobalVars.iprint=int(options.dest_iprint)
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.iprint=',GlobalVars.iprint

  # number of nodes
  if (options.dest_nnodes is not None) :
    GlobalVars.nnodes=int(options.dest_nnodes)
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.nnodes=',GlobalVars.nnodes

  # about submitting job
  if (options.dest_nosub is not None) :
    GlobalVars.submit_job=options.dest_nosub
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.submit_job=',GlobalVars.submit_job

  #  bout no mpi stuff for pam 
  if (options.dest_is_serial is not None) :
    GlobalVars.is_serial=options.dest_is_serial
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.is_serial=',GlobalVars.is_serial

  # ... save the pam directory
  if (options.dest_pam is not None) :
   # GlobalVars.pam_script=options.dest_pam
    GlobalVars.save_pam_script(options.dest_pam)
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.pam_script=',GlobalVars.pam_script

  # ... deal with pam arguments
  if (options.dest_pam_args is not None) :
    GlobalVars.pam_arguments=options.dest_pam_args
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.pam_arguments=',GlobalVars.pam_arguments
  else :
    GlobalVars.submit_job=False
#    print 'warning: no pam arguments in the script - not submitting ! '
#    print 'type: "',sys.argv[0],'-h" for help'
#    parser.error("missing pam arguments !")
    #print options.get_help()

  # ... deal with queue_name
  if (options.dest_queue_name is not None) :
    GlobalVars.queue_name=options.dest_queue_name
    GlobalVars.queue_name_specified=True
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.queue_name=',GlobalVars.queue_name
      print 'assigned GlobalVars.queue_name_specified=',GlobalVars.queue_name_specified

  # ... deal with the queue_system
  if (options.dest_queue_system is not None) :
    GlobalVars.queue_system=options.dest_queue_system
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.queue_system=',GlobalVars.queue_system

  # ... deal with  number of processors per 1 cluster node
 # print '* options.dest_ppn=',options.dest_ppn
  if (options.dest_ppn is not None) :
  #  print 'options.dest_ppn=',options.dest_ppn
    GlobalVars.proc_per_node=int(options.dest_ppn)
    if GlobalVars.iprint >= 7:
      print 'assigned GlobalVars.proc_per_node=',GlobalVars.proc_per_node

  # ... deal with queue_time
  if (options.dest_queue_time is not None) :
    GlobalVars.queue_time=options.dest_queue_time
    GlobalVars.queue_time_specified=True
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.queue_time=',GlobalVars.queue_name
      print 'assigned GlobalVars.queue_time_specified=',GlobalVars.queue_time_specified

  # ... deal with creating variables_defaults file
  if (options.dest_prep_config is not None) :
    GlobalVars.prepare_config_file=options.dest_prep_config
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.prepare_config_file=',GlobalVars.prepare_config_file


  # ... must be at leat one argument - pam_args !
  #if len(args) < 1:
  #  parser.error("incorrect number of arguments")

  if (len(parser.largs)>0) :
    print 'error: leftover arguments found :',parser.largs 
    parser.error("incorrect passing arguments, all must be under flags !")

  if GlobalVars.iprint > 8 :
    print '\n processing arguments' 
    print 'len(sys.argv)=', len(sys.argv)
    for i in range(len(sys.argv)) :
      print i,sys.argv[i]
		
  #------------------------------------------------------------------------
  #      collect all left pamq.py argumnents into  pam_arguments
  #------------------------------------------------------------------------
  for item in parser.largs:
    GlobalVars.pam_arguments= GlobalVars.pam_arguments+' '+item
  GlobalVars.pam_arguments=' '+GlobalVars.pam_arguments

  if GlobalVars.iprint > 7:
    print 'extracted pam arguments:',GlobalVars.pam_arguments

  # set mpi for pam 
  if not GlobalVars.is_serial :
    GlobalVars.set_pam_mpi()

  # if prepare file, do not submit
  if GlobalVars.prepare_config_file:
    GlobalVars.submit_job=False
    GlobalVars.create_configure_file()
    if GlobalVars.iprint > 7:
      print 'preparing configure file-not submitting job'
  
#########################################################################
#
#                        MAIN FUNCTION
#
# MI: I had to remove error catching to be able to localize a bug in the script
#
#########################################################################
def main(*args):
#  try:
    global GlobalVars
    GlobalVars = PamQVariables() 

    # some code here
   # set_defaults()

  #  if len(sys.argv) < 2 :
  #    pamq_usage()
  #    print 'error exit !'
  #    sys.exit(2)

    process_pamq_arguments()

    if (not GlobalVars.prepare_config_file) :
      process_queue_script()

#  except:
    # handle some exceptions
#    print 'caught exception...'
#  else:
#    return 0 # exit errorlessly
 
###########################################################################
###                   process the main function
###########################################################################
if __name__ == '__main__':
    sys.exit(main(*sys.argv))

