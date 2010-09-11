#!/usr/bin/python 
##! /usr/bin/env python        
#
                                      
############################################################################################
#
#  Hoping for someone (a student ?) to turn the bash-based configure script 
# into the python equivalent.
#
############################################################################################

import  sys, os, string, subprocess, optparse, distutils.util

# Search in the system for these compilers
F_compilers=['ifort','gfortran','pgf90','g95','lfc']
C_compilers=['icc','gcc','pgcc']

def guess_system():
# guess the operating system
  print "guessing oper.system",os.uname()
  print "is it Windows ? --> ",is_windows()

def is_windows():
  # returns true if this machine is Windows operating system
  platform=distutils.util.get_platform()
  if platform.find("win") == -1 :
    return False
  else:
    return True



#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class ConfigureVariables:
  # ... class containing all important variables, methods for submitting 
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def __init__(self):
  #
  #           ... initialize the global set of variables ...
  #
    self.Bash_Configure_Arguments=""

    #  pure python mode or bash shell mode
    self.PythonMode=False
    
	# installation directory
    self.InstallDirectory=""

	#Fortran compiler
    self.Fortran_Compiler=""
	
	#C-compiler
    self.C_Compiler=""
	
    # print level default
    self.iprint=1

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_configure_arguments():
#-----------------------------------------------------
#
#-----------------------------------------------------

  print 'processing configure arguments....'

  from optparse import OptionParser
  usage="usage: %prog arg1 arg2..."
  parser = OptionParser(usage)

  #... get in the print level 
  iprint_label='script verbose print level, default:'+str(GlobalVars.iprint)
  parser.add_option("-i", "--iprint", type="int", nargs=1,  
        action="store", dest="dest_iprint", default=GlobalVars.iprint, 
        help=iprint_label, metavar="INT_NUM")
  
  #  fill configure_arguments
  label_configure_arguments=" configure arguments (like '"' -nopt ....'"')"
  parser.add_option("--conf_args", type="string", dest="dest_configure_args",
                  help=label_configure_arguments, metavar="'"'arg1 arg2 ... '"' ")

  # is this pure python script configuration ?
  use_python_label='use python only based script, otherwise loading bash-based configure ' 
  parser.add_option("--py", action="store_true", default=GlobalVars.PythonMode,
                     dest="dest_py_config", help=use_python_label)
 
    # assign aruments...
  (options, args) = parser.parse_args()

  # extract read arguments
  print '\n processing read arguments:  len(sys.argv)=', len(sys.argv)

  # ... handle iprint
  if (options.dest_iprint is not None) :
    GlobalVars.iprint=int(options.dest_iprint)
    if GlobalVars.iprint > 7:
      print 'assigned GlobalVars.iprint=',GlobalVars.iprint
    
  GlobalVars.PythonMode=options.dest_py_config
  print 'assigned GlobalVars.PythonMode=',GlobalVars.PythonMode
  if (GlobalVars.PythonMode) :
    print 'Pure python mode of the configure script ! Program it !'
  else :
    print 'Using the old bash written trunk/configure script'
    
  if (options.dest_configure_args is not None) :
    GlobalVars.Bash_Configure_Arguments=options.dest_configure_args
    print 'Read old bash-configure arguments:',GlobalVars.Bash_Configure_Arguments
  else :
    print 'No bash-configure arguments...'
  
 # sys.exit() 

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def call_bash_configure():
#-----------------------------------------------------
#   calls the external bash-based configure script
#-----------------------------------------------------
    print 'Doing the old-style bash configure ...\n'

    bash_script="cd ..; ./configure "
    #bash_script_arguments="--noopt"
    if (GlobalVars.Bash_Configure_Arguments is not None) :
      bash_script_arguments=GlobalVars.Bash_Configure_Arguments
      bash_script=bash_script+bash_script_arguments
      print '...old configure stuff:',bash_script

    # old-style configuration...
    subprocess.call(bash_script, executable="bash", shell=True)


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#class PrintLevel:
#  def __init__(self):
#    self.iprint=0
#  def set_iprint(self):
#    self.iprint = iprint_read


#### check where configure.py sits - must be the main trunk
#def check_position():
#  print 'checking position of configure.......'

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def call_python_configure():
#-----------------------------------------------------
#
#      do the job of the old bash configure
#
#  create the trunk/Makefile.config file
#
#-----------------------------------------------------
  print 'write me ! :)'
 
  guess_system()

  # take care of the installation directory ...

  # scratch direcory ...

#####################################################################################################
#
#  Define main function
#
#####################################################################################################
def main(*args):
  global GlobalVars
  GlobalVars = ConfigureVariables() 

  process_configure_arguments()

  if (GlobalVars.PythonMode) :
    call_python_configure()
  else :
    call_bash_configure()  

###########################################################################
###                   process the main function
###########################################################################
if __name__ == '__main__':
    sys.exit(main(*sys.argv))
