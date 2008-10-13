#! /usr/bin/env python
#
#  script for submitting DIRAC pam job into the queue system
#
# writen by Miro ILIAS,2008 after Luuks pamq
#
#  still in progrss
#

import re, sys, os, string


iprint = 12

job_name='queue test'
job_inp='input_test'

job_out='job_output'
#loc=`pwd`
loc='`pwd`'

queue_system='pbs'

queue_name='default'

nnodes=2
nodes_per_proc=4

#pam_dir="`pwd`"
pam_script_name="pam"
pam_script_dir="~/bin"
pam_script=pam_script_dir + "/" + pam_script_name

queue_script="pamq_batch_script"


#class queue_system:
       # initialization
#	def __init__():

def set_defaults():
	iprint=5

def pamq_usage():
	print 'Usage: pamq.py pamq-flags pam-flags [mol] [inp]'
	print '             '
	print '  special pamq flags:'
	print '             '
	print '-nosub     : do not submit, only prepare the script for batch'
	print '-i <number> : verbose print level, default',iprint
	print '-qsys <queue_system> : batch system, default:','"',queue_system,'"'
	print '-pamdir <path/pam> : location of pam, default:',pam_script
	print '             '

def sumbit_job():  
	if iprint > 7:
		print 'submitting job'

def process_pbs_script():
	print 'prepare_pbs_script...'
	f=open(queue_script, 'w')
	s = ''
	s += '''
#!/bin/bash
#
# Use controls since these are directly readable fromthe batch script
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
echo '   PBS_O_HOME='\$PBS_O_HOME
echo '   PBS_O_HOST='\$PBS_O_HOST
echo '   PBS_SERVER='\$PBS_SERVER
echo '  PBS_O_QUEUE='\$PBS_O_QUEUE
echo 'PBS_O_WORKDIR='\$PBS_O_WORKDIR
echo '  PBS_JOBNAME='\$PBS_JOBNAME
echo '    PBS_JOBID='\$PBS_JOBID
echo ' PBS_NODEFILE='\$PBS_NODEFILE
echo '------------------------------------'
#
#   ... change to original directory
#
cd $loc
	'''

	f.write(s)
	f.close()

def process_loadleveler_script():
	print 'prepare_loadleveler_script...'


def process_bsub_script():
	print 'prepare_loadleveler_script...'

def process_queue_script():
	if iprint > 7:
		print 'Going to prepare script for batch system...'
	if queue_system=='pbs':
		process_pbs_script()	
 	elif queue_system=='loadleveler': 
		process_loadleveler_script()	
 	elif queue_system=='bsub': 
		process_bsub_script()	
	else:
		print 'not defined queue_system=',queue_system
		exit(0)


def process_pamq_arguments():

	if iprint > 7:
	        print '\n processing arguments' 
		print 'len(sys.argv)=', len(sys.argv)
		for i in range(0,len(sys.argv)):
			print 'sys.argv[',i,']=', sys.argv[i]

#print ' ... python script'
#
#                   MAIN SCRIPT
#
#
if len(sys.argv) < 2 :
	pamq_usage()
	exit

set_defaults()

process_pamq_arguments()

process_queue_script()

#raw_input("press Enter>")
