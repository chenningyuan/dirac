#! /usr/bin/env python
#
#  script for submitting DIRAC pam job into the queue system
#
# writen by Miro ILIAS,2008 after Luuks pamq
#
#  still in progrss
#

import  sys, os, string


global iprint, job_name, job_inp, job_out, loc, queue_system, queue_name, nnodes, \
  nodes_per_proc, pam_script , queue_script, pam_arguments, submit_job, home_dir, \
  pam_script_args

def set_defaults():
  global iprint, job_name, job_inp, job_out, loc, queue_system, queue_name, nnodes, \
 		nodes_per_proc, pam_script, queue_script, pam_arguments, submit_job, home_dir

  iprint=12

  pam_arguments=' '

  submit_job=1

  job_name='queue_test'
  job_inp='input_test'

  job_out='job_output'
#loc=`pwd`
  local_dir='`pwd`'

  queue_system='pbs'

  queue_name='default'

  nnodes=2
  nodes_per_proc=4

#pam_dir="`pwd`"
  pam_script_name="pam"
  pam_script_dir="~/bin"
  pam_script=pam_script_dir + "/" + pam_script_name

  queue_script="pamq_batch_script"

#  os.system('`echo $HOME`')
  os.system('echo "home ..." $HOME')
  os.system('echo "pwd  ..." `pwd`')
  #home_dir=os.system('echo `pwd`')
 # home_dir=os.path.expanduser('~')
  home_dir=os.getenv("PWD")
  #home_dir=string.join(os.system('echo `pwd`'))
 # home_dir=sys.pwd

  if iprint >8 :
   print ' setting home directory:',home_dir

#class queue_system:
       # initialization
#	def __init__():

#def set_defaults():
#	iprint=5

def pamq_usage():
	print 'Usage: pamq.py queue-flags pam-flags [mol] [inp]'
	print '             '
	print '  own pamq queue-flags:'
	print '             '
	print '-nosub     : do not submit, only prepare the script for batch'
	print '-i <number> : verbose print level, default',iprint
	print '-qsys <queue_system> : batch system, default:','"',queue_system,'"'
	print '-jobname <queue job name> : name of the job batch system, default:... '
	print '-pamdir <path/pam> : location of pam, default:',pam_script
	print '             '
	print 'pam-flags are obtained using pam...'

def sumbit_job():  
	if iprint > 7:
		print 'submitting job'

def process_pbs_script():
  global queue_script

  if iprint >8 :
    print 'prepare_pbs_script...'

  # unite into one local string
  pam_script_args=pam_script+pam_arguments
  if iprint >8 :
   print ' home directory:',home_dir

  f=open(queue_script, 'w')
  s = ''
  s += '''
#!/bin/bash
#
# Queue controls are directly readable from the submitting batch script
#
# job name
#PBS -N %s
  ''' % job_name
  s += '''
#
# nodes
#PBS -l nodes=%i
  ''' % nnodes
  s += '''
#
# queue
##PBS -q %s
  ''' % queue_name
  s += '''
#
# unite stdout with stderr
#PBS -j oe
#
#
#   queue job standard output - stdout
#PBS -o %s
  ''' % job_out
  s += '''
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
#   ... always change to the home directory
#
cd %s
  ''' % home_dir
  s += '''
#
#   ... execute own DIRAC pam script
#
%s 
  ''' % pam_script_args
  s += '''
#
#
#
exit
  '''

  f.write(s)
  f.close()

 # submit job ...
  if submit_job:
    if iprint > 7:
      print 'going to submit job with script ',queue_script
    #os.system('qsub queue_script'),queue_script 
    os.system("qsub "+queue_script)
    #os.system('ls')

#
#
#
def process_loadleveler_script():
  if iprint > 7:
    print 'prepare_loadleveler_script...'


def process_bsub_script():
	if iprint > 7:
		print 'preparing bsub script...'

def process_queue_script():
	if iprint > 7:
		print 'Going to prepare script for batch system...'

	if queue_system=='pbs':
		process_pbs_script()	
 	elif queue_system=='loadleveler': 
		process_loadleveler_script()	
 	elif queue_system=='bsub': 
		process_bsub_script()	
 	elif queue_system=='xgrid': 
#		process_xgrid_script()	
	else:
		print 'not defined queue_system=',queue_system
		exit(0)


#
# process script arguments 
#
def process_pamq_arguments():

  global iprint, pam_arguments, queue_system, job_name

  if iprint > 7:
    print '\n processing arguments' 
    print 'len(sys.argv)=', len(sys.argv)
		

  i=0
	#range_m1=range(0,len(sys.argv))-1
	#print 'range(0,len(sys.argv))=',range(0,len(sys.argv))

#	for i in range(0,len(sys.argv)):
  i=0
  while  i <  len(range(0,len(sys.argv))) :
    if iprint > 7:
      print 'sys.argv[',i,']=', sys.argv[i]
    if sys.argv[i]=='-nosub':
      print 'not submitting !'
    elif sys.argv[i]=='-i':	
      i=i+1
      iprint=sys.argv[i]	
      print 'iprint=',iprint
    elif sys.argv[i]=='-qsys':
      i=i+1	
      queue_system=sys.argv[i]	
      print 'obtained queue_system=',queue_system
    elif sys.argv[i]=='-jobname':
      i=i+1	
      job_name=sys.argv[i]
      print 'obtained queue_system=',job_name
    else:
      # store all remaining arguments
      if i>0:
        pam_arguments=pam_arguments + ' ' + sys.argv[i]
	if iprint > 7:
	  print 'pam arg storing:',sys.argv[i]
	  print 'filled pam arguments:',pam_arguments
			
    i = i + 1

  if iprint > 7:
    print 'remaining pam arguments:'
    print pam_arguments

  # unite into one local string
	#  pam_script_args=pam_script+pam_arguments

	

#print ' ... python script'
#
#                   MAIN SCRIPT
#
#

set_defaults()

if len(sys.argv) < 2 :
  pamq_usage()
  exit

process_pamq_arguments()

process_queue_script()

