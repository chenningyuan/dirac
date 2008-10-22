#! /usr/bin/env python
#
#  script for submitting DIRAC pam job into the queue system
#
# writen by Miro ILIAS,2008 after Luuks pamq
#
#  still in progrss
#

import  sys, os, string

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class PamQVariables:

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def __init__(self):
  #
  #  initialize set of variables ...
  #

    # print level
    self.iprint=1

    self.pam_arguments=' '
    self.submit_job=True

    self.job_name='queue_job.'+str(os.getpid())

    self.job_inp='input_test'

    self.job_out='queue_job.'+str(os.getpid())+'.log_out'

    self.queue_system='pbs'

   #  self.queue_name='default'
   # self.queue_name='workq'
    self.queue_name='express'
    self.queue_name_specified=False

    # time 
    self.queue_time='0:20'
    self.queue_time_specified=False

    # ... is serial job (ie no mpi adding to pam arguments !)
    self.is_serial=False

    self.nnodes=1
    self.proc_per_node=4
    self.pam_script_name="pam"
    self.pam_script_dir="~/bin"
    self.pam_script=self.pam_script_dir + "/" + self.pam_script_name

    self.queue_script="pamq."+str(os.getpid())+".batch_script"

    self.home_dir=os.getenv("PWD")

    if self.iprint >8 :
      print ' setting home directory:',self.home_dir
      print 'INITIALIZATION of class object done'

  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def set_pam_mpi(self):
  #
  # class method: set up mpi parameter for pam script
  #
    if GlobalVars.iprint > 7:
      print 'set_pam_mpi: extracted pam arguments:', GlobalVars.pam_arguments

    # extend pam flags with the necessary -mpi arguments
    mpi_procs=GlobalVars.nnodes * GlobalVars.proc_per_node
    GlobalVars.pam_arguments=' -mpi '+str(mpi_procs)+GlobalVars.pam_arguments

    if GlobalVars.iprint > 7:
      print 'set_pam_mpi: mpi-extended pam arguments:', GlobalVars.pam_arguments

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def set_defaults():
#
# TODO: read default data from some file
#
  global GlobalVars
  if GlobalVars.iprint > 7:
     print 'in set_defaults...'

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def pamq_usage():
#
#
#
  global GlobalVars

  print 'Usage: pamq.py queue-flags pam-flags [mol] [inp]'
  print '             '
  print '  own pamq queue-flags:'
  print '             '
  print '-help                    : print out pam flags'
  print '-nosub                   : do not submit, only prepare the script for batch'
  print '-serial                  : serial job (not mpi), default:',GlobalVars.is_serial
  print '-i <number>              : verbose print level, default: ',GlobalVars.iprint
  print '-qsys <queue_system>     : batch system, default :','"',GlobalVars.queue_system,'"'
  print '-q <queue>               : specify the queue, default :','"',GlobalVars.queue_name,'"'
  print '-t <time>                : specify the time for queue job, default :',GlobalVars.queue_time
  print '-ppn  <number>           : number of processors per 1 node, default:',GlobalVars.proc_per_node
  print '-nodes <number>          : number of nodes, default :',GlobalVars.nnodes
  print '-jobname <job name>      : name of the job batch system, default :','"',GlobalVars.job_name,'"'
  print '-pam <path/pam_script>   : full location of pam, default :',GlobalVars.pam_script
  print '             '
  print 'pam-flags are obtained using pam...'

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_pbs_script():
#
#
#   prepare the file for script and submit
#
#
  global GlobalVars

  if GlobalVars.iprint >8 :
     print 'in prepare_pbs_script...'

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
# Queue controls are directly readable from the submitting batch script
#
# job name
#PBS -N %s
  ''' % GlobalVars.job_name

#### if not serial job
  if GlobalVars.is_serial:
    s += '''\
#
# number nodes DEACTIVATED (serial job!)
#
##PBS -l nodes=%i
    ''' % GlobalVars.nnodes
  else:
    # ... parallel job
    s += '''\
#
# number nodes (NOTE: different from number of processors specified by -mpi pam flag !)
#
#PBS -l nodes=%i
''' % GlobalVars.nnodes

### if queue_name is specified, add it
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
# queue name not specified, using system default, PBS command inactive
#
##PBS -q %s
    ''' % GlobalVars.queue_name
#
# if time is specified in the input, add it
#
  if GlobalVars.queue_time_specified:
    s += '''\
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
#PBS -o %s
''' % GlobalVars.job_out
  s += '''\
#
#
# ... print out important PBS variables into standard output
#
echo '------------------------------------'
echo '   PBS_O_HOME='$PBS_O_HOME
echo '   PBS_O_HOST='$PBS_O_HOST
echo '   PBS_SERVER='$PBS_SERVER
echo '  PBS_O_QUEUE='$PBS_O_QUEUE
echo 'PBS_O_WORKDIR='$PBS_O_WORKDIR
echo '  PBS_JOBNAME='$PBS_JOBNAME
echo '    PBS_JOBID='$PBS_JOBID
echo ' PBS_NODEFILE='$PBS_NODEFILE
echo '------------------------------------'
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
    if GlobalVars.iprint > 7:
      print 'going to submit job with prepared script ',GlobalVars.queue_script
    os.system("qsub "+GlobalVars.queue_script)
  else:
    print 'job not submitting!'


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_loadleveler_script():
#
# TODO : add !
#
  if GlobalVars.iprint > 7:
    print 'prepare_loadleveler_script...'

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
   print 'not defined queue_system=',GlobalVars.queue_system
   exit(-1)


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_pamq_arguments():
#
#             process script arguments 
#
  global GlobalVars

  if GlobalVars.iprint > 7:
    print '\n processing all read arguments' 
    print 'len(sys.argv)=', len(sys.argv)
		
  # ... cycle - get flags
  i=0
  while  i <  len(range(0,len(sys.argv))) :
    # print the current parameter
    if GlobalVars.iprint > 7:
      print 'sys.argv[',i,']=', sys.argv[i]
    #
    #
    #
    if sys.argv[i]=='-nosub':
      GlobalVars.submit_job=False
      if GlobalVars.iprint > 7:
        print 'not submitting ! submit_job=',GlobalVars.submit_job
    #
    #
    #
    if sys.argv[i]=='-serial':
      GlobalVars.is_serial=True
      if GlobalVars.iprint > 7:
        print 'no mpi job ! is_serial=',GlobalVars.is_serial
    #
    #
    #
    elif sys.argv[i]=='-i':	
      i=i+1
      GlobalVars.iprint=int(sys.argv[i])
      if GlobalVars.iprint > 0:
        print 'obtained iprint=',GlobalVars.iprint
        print 'obtained queue_system=',GlobalVars.queue_system
    #
    #
    #
    elif sys.argv[i]=='-qsys':
      i=i+1	
      GlobalVars.queue_system=sys.argv[i]
      if GlobalVars.iprint > 7:
        print 'obtained queue_system=',GlobalVars.queue_system
    #
    #
    #
    elif sys.argv[i]=='-q':
      i=i+1	
      GlobalVars.queue_name=sys.argv[i]
      if GlobalVars.iprint > 7:
        print 'obtained queue_name=',GlobalVars.queue_name
      GlobalVars.queue_name_specified=True
    #
    #
    #
    elif sys.argv[i]=='-t':
      i=i+1	
      GlobalVars.queue_time=sys.argv[i]
      if GlobalVars.iprint > 7:
        print 'obtained time=',GlobalVars.queue_time
      GlobalVars.queue_time_specified=True
    #
    #
    #
    elif sys.argv[i]=='-ppn':
      i=i+1	
      GlobalVars.proc_per_node=int(sys.argv[i])
      if GlobalVars.iprint > 7:
        print 'obtained proc_per_node=',GlobalVars.proc_per_node
    #
    #
    #
    elif sys.argv[i]=='-nodes':
      i=i+1	
      GlobalVars.nnodes=int(sys.argv[i])
      if GlobalVars.iprint > 7:
        print 'obtained nodes=',GlobalVars.nnodes
    #
    #
    #
    elif sys.argv[i]=='-jobname':
      i=i+1	
      GlobalVars.job_name=sys.argv[i]
      if GlobalVars.iprint > 7:
        print 'obtained queue_system=',GlobalVars.job_name
    #
    #
    #
    elif sys.argv[i]=='-pam':
      i=i+1	
      GlobalVars.pam_script=sys.argv[i]
      if GlobalVars.iprint > 7:
        print 'obtained pam script=',GlobalVars.pam_script
    else:
      # store all remaining arguments
      if i>0:
        GlobalVars.pam_arguments=GlobalVars.pam_arguments + ' ' + sys.argv[i]
	if GlobalVars.iprint > 7:
	  print 'pam arg storing:',sys.argv[i]
	  print 'filled pam arguments:',GlobalVars.pam_arguments
			
    i = i + 1
  #
  #
  #
  if GlobalVars.iprint > 7:
    print 'extracted pam arguments:'
    print GlobalVars.pam_arguments

  if not GlobalVars.is_serial:
    GlobalVars.set_pam_mpi()

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
    set_defaults()

    if len(sys.argv) < 2 :
      pamq_usage()
      print 'error exit !'
      sys.exit(2)

    process_pamq_arguments()

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

