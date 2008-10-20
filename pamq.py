#! /usr/bin/env python
#
#  script for submitting DIRAC pam job into the queue system
#
# writen by Miro ILIAS,2008 after Luuks pamq
#
#  still in progrss
#

import  sys, os, string

class PamQVariables:
  def __init__(self):
    self.iprint=12
    self.pam_arguments=' '
    self.submit_job=1
    self.job_name='queue_test'
    self.job_inp='input_test'
    self.job_out='job_output'
#loc=`pwd`
    self.local_dir='`pwd`'
    self.queue_system='pbs'
    self.queue_name='default'
    self.nnodes=1
    self.nodes_per_proc=4
#pam_dir="`pwd`"
    self.pam_script_name="pam"
    self.pam_script_dir="~/bin"
    self.pam_script=self.pam_script_dir + "/" + self.pam_script_name
#  queue_script="pamq_batch_script"+str(os.getpid())
    self.queue_script="pamq_batch_script"

#    os.system('echo "home ..." $HOME')
#    os.system('echo "pwd  ..." `pwd`')
  #home_dir=os.system('echo `pwd`')
 # home_dir=os.path.expanduser('~')
    self.home_dir=os.getenv("PWD")
  #home_dir=string.join(os.system('echo `pwd`'))
 # home_dir=sys.pwd
    if self.iprint >8 :
      print ' setting home directory:',self.home_dir
    print 'INITIALIZATION of class object done'

def set_defaults():
  global GlobalVars
  print 'in set_defaults...'

def pamq_usage():
  global GlobalVars
  print 'Usage: pamq.py queue-flags pam-flags [mol] [inp]'
  print '             '
  print '  own pamq queue-flags:'
  print '             '
  print '-help      : print out pam flags'
  print '-nosub     : do not submit, only prepare the script for batch'
  print '-i <number> : verbose print level, default',GlobalVars.iprint
  print '-qsys <queue_system> : batch system, default:','"',GlobalVars.queue_system,'"'
  print '-ppn  <number> :  number of processors per 1 node, default:',GlobalVars.nodes_per_proc
  print '-nodes <number> :  number of nodes, default:',GlobalVars.nnodes
  print '-jobname <queue job name> : name of the job batch system, default:... '
  print '-pam <path/pam> : location of pam, default:',GlobalVars.pam_script
  print '             '
  print 'pam-flags are obtained using pam...'

def sumbit_job():  
	if iprint > 7:
		print 'submitting job'

def process_pbs_script():
  global GlobalVars

  if GlobalVars.iprint >8 :
    print 'prepare_pbs_script...'

  # unite into one local string
  GlobalVars.pam_script_args=GlobalVars.pam_script+GlobalVars.pam_arguments

  if GlobalVars.iprint >8 :
   print ' home directory:',GlobalVars.home_dir
   print ' GlobalVars.pam_script_args=',GlobalVars.pam_script_args

  f=open(GlobalVars.queue_script, 'w') # open file with queue commands

  s=''
  s += '''\
#!/bin/bash
#
# Queue controls are directly readable from the submitting batch script
#
# job name
#PBS -N %s
''' % GlobalVars.job_name
#  print 'GlobalVars.job_name=',GlobalVars.job_name
  s += '''\
#
# nodes
#PBS -l nodes=%i
''' % GlobalVars.nnodes
  s += '''\
#
# queue
##PBS -q %s
''' % GlobalVars.queue_name
  s += '''\
#
# unite stdout with stderr
#PBS -j oe
#
#
#   queue job standard output - stdout
#PBS -o %s
''' % GlobalVars.job_out
  s += '''\
#
#
# ... print out important PBS variables
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
#
#
#
def process_loadleveler_script():
  if GlobalVars.iprint > 7:
    print 'prepare_loadleveler_script...'

def process_bsub_script():
  if GlobalVars.iprint > 7:
    print 'preparing bsub script...'

def process_queue_script():
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


#
#             process script arguments 
#
def process_pamq_arguments():
  global GlobalVars

  if GlobalVars.iprint > 7:
    print '\n processing arguments' 
    print 'len(sys.argv)=', len(sys.argv)
		
	#range_m1=range(0,len(sys.argv))-1
	#print 'range(0,len(sys.argv))=',range(0,len(sys.argv))

#	for i in range(0,len(sys.argv)):
  i=0
  while  i <  len(range(0,len(sys.argv))) :
    if GlobalVars.iprint > 7:
      print 'sys.argv[',i,']=', sys.argv[i]
    if sys.argv[i]=='-nosub':
      GlobalVars.submit_job=0
      print 'not submitting ! submit_job=',GlobalVars.submit_job
    elif sys.argv[i]=='-i':	
      i=i+1
      GlobalVars.iprint=sys.argv[i]
      print 'iprint=',GlobalVars.iprint
    elif sys.argv[i]=='-qsys':
      i=i+1	
      GlobalVars.queue_system=sys.argv[i]
      print 'obtained queue_system=',GlobalVars.queue_system
    elif sys.argv[i]=='-ppn':
      i=i+1	
      GlobalVars.nodes_per_proc=sys.argv[i]
      print 'obtained nodes_per_proc=',GlobalVars.nodes_per_proc
    elif sys.argv[i]=='-nodes':
      i=i+1	
      GlobalVars.nnodes=sys.argv[i]
      print 'obtained nodes=',GlobalVars.nnodes
    elif sys.argv[i]=='-jobname':
      i=i+1	
      GlobalVars.job_name=sys.argv[i]
      print 'obtained queue_system=',GlobalVars.job_name
    elif sys.argv[i]=='-pam':
      i=i+1	
      GlobalVars.pam_script=sys.argv[i]
      print 'obtained pam script=',GlobalVars.pam_script
    else:
      # store all remaining arguments
      if i>0:
        GlobalVars.pam_arguments=GlobalVars.pam_arguments + ' ' + sys.argv[i]
	if GlobalVars.iprint > 7:
	  print 'pam arg storing:',sys.argv[i]
	  print 'filled pam arguments:',GlobalVars.pam_arguments
			
    i = i + 1

  if GlobalVars.iprint > 7:

    print 'extracted pam arguments:'
    print GlobalVars.pam_arguments

  # unite into one local string
  #  pam_script_args=pam_script+pam_arguments

#########################################################################
#                        MAIN FUNCTION
#########################################################################
def main(*args):
#  try:
    global GlobalVars
    GlobalVars = PamQVariables() 
   # global GlobalVars

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
 
####
###    process the main function
####
if __name__ == '__main__':
    sys.exit(main(*sys.argv))

