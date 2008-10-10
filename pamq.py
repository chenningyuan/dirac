#! /usr/bin/env python
#
#  script for submitting DIRAC pam job into the queue system
#
#
#
#

import re, sys, os, string


iprint = 12
qsys='pbs'

#class queue_system:
       # initialization
#	def __init__():

def set_defaults():
	iprint=5

def pamq_usage():
	print 'Usage: pamq.py pamq-flags pam-flags [mol] [inp]'
	print '             '
	print '  own pamq flags:'
	print '-nosub     : do not submit, only prepare the script for batch'
	print '-i <number> : verbose print level, default',iprint
	print '-qsys <queue_system> : batch system, default:','"',qsys,'"'

def sumbit_job():  
	if iprint > 7:
		print 'submitting job'

def prepare_pbs_script():
	print 'prepare_pbs_script...'

def prepare_loadleveler_script():
	print 'prepare_loadleveler_script...'

def prepare_bsub_script():
	print 'prepare_loadleveler_script...'

def process_pamq_arguments():

	if iprint > 7:
	        print '\n processing arguments' 
		print 'len(sys.argv)=', len(sys.argv)
		for i in range(0,len(sys.argv)):
			print 'sys.argv[',i,']=', sys.argv[i]


#print ' ... python script'

if len(sys.argv) < 2 :
	pamq_usage()
	exit

set_defaults()
process_pamq_arguments()


#raw_input("press Enter>")
