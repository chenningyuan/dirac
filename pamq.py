#! /usr/bin/env python
#
#  script for submitting DIRAC pam job into the queue system
#
# writen by Miro ILIAS,2008 after Luuks pamq
#
#  still in progrss
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

    # ... is pure serial job (no mpi added to pam)
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
      print 'setting mpi for pam' 
      print 'set_pam_mpi: extracted pam arguments:', GlobalVars.pam_arguments
    # print GlobalVars.pam_arguments

    # extend pam flags with the necessary -mpi arguments
    mpi_procs=GlobalVars.nnodes * GlobalVars.proc_per_node
    GlobalVars.pam_arguments=' -mpi '+str(mpi_procs)+GlobalVars.pam_arguments

    if GlobalVars.iprint > 7:
      print 'set_pam_mpi: mpi-extended pam arguments:', GlobalVars.pam_arguments

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
#
# queue name not specifies, using system default and inactive
#
##PBS -q %s
    ''' % GlobalVars.queue_name
#
# if time is specified in the input, add it
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
#----------------------------------------------------------------------------------------
#
#   process script arguments 
#
#----------------------------------------------------------------------------------------
  global GlobalVars
  from optparse import OptionParser

  usage = "usage: %prog [pamq_options] pam_arg1 pam_arg2..."
  parser = OptionParser(usage)

  #  fill GlobalVars.pam_script
  pam_script_add="location of pam script, default:"+GlobalVars.pam_script
  parser.add_option("-p", "--pam", type="string", nargs=1,dest="dest_pam",default=GlobalVars.pam_script,
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

  #  ... get number of nodes
  nodes_label='number of reserved nodes, default:'+str(GlobalVars.nnodes)
  parser.add_option("-n", "--nodes", type="int",  nargs=1,
                  action="store",dest="dest_nnodes", default=GlobalVars.nnodes, 
                  help=nodes_label, metavar="INT_NUM")

  # ... GlobalVars.queue_system
  queue_system_label='specify the queue system, default :'+"'"+GlobalVars.queue_system+"'"
  parser.add_option( "--qsys", type="str", nargs=1,
                  action="store",dest="dest_queue_system", default=GlobalVars.queue_system, 
                  help=queue_system_label,metavar="STRING")

  # ... GlobalVars.queue_name
  queue_label='specify the queue, default :'+GlobalVars.queue_name
  parser.add_option("-q", "--queue", type="string", nargs=1,
                  action="store", dest="dest_queue_name", default=GlobalVars.queue_name, \
                  help=queue_label,metavar="STRING")

  # ... GlobalVars.queue_time
  queue_time_label='specify the time in the queue, default:'+GlobalVars.queue_time
  parser.add_option("-t", "--time", type="string", nargs=1,
                 # action="store", dest="dest_queue_time", default=GlobalVars.queue_time, \
                  action="store", dest="dest_queue_time", 
                  help=queue_time_label,metavar="STRING")

  # ... GlobalVars.proc_per_node
  queue_ppn_label='specify number of processor per node, default:'+str(GlobalVars.proc_per_node)
  parser.add_option("--ppn", type="int", nargs=1,
                 # action="store", dest="dest_ppn", default=GlobalVars.proc_per_node,
                  action="store", dest="dest_ppn", 
                  help=queue_ppn_label, metavar="STRING")

  ###################################################################################################
  #
  #                                    ... finally assign
  #
  ###################################################################################################
  (options, args) = parser.parse_args()

# .. now assign individual arguments 

  if GlobalVars.iprint > 5:
    print '\n processing arguments:  len(sys.argv)=', len(sys.argv) 

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

  if (options.dest_pam is not None) :
    GlobalVars.pam_script=options.dest_pam
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.pam_script=',GlobalVars.pam_script

  # ... deal with pam arguments
  if (options.dest_pam_args is not None) :
    GlobalVars.pam_arguments=options.dest_pam_args
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.pam_arguments=',GlobalVars.pam_script
  else :
    GlobalVars.submit_job=False
    print 'warning: no pam arguments in the script - not submitting ! '

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

  # ... deal with processors per node
  if (options.dest_ppn is not None) :
    GlobalVars.proc_per_node=options.dest_ppn
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.proc_per_node=',GlobalVars.proc_per_node

  # ... deal with queue_time
  if (options.dest_queue_time is not None) :
    GlobalVars.queue_time=options.dest_queue_time
    GlobalVars.queue_time_specified=True
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.queue_time=',GlobalVars.queue_name
      print 'assigned GlobalVars.queue_time_specified=',GlobalVars.queue_time_specified

#  if len(args) < 3:
#    parser.error("incorrect number of arguments")

  if (len(parser.largs)>0) :
    print 'warning: leftover arguments:',parser.largs  

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

