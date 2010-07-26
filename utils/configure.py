#! /usr/bin/env python        
#
                                      
############################################################################################
#
#  Hoping for someone (a student ?) to turn the bash-based configure script 
# into the python equivalent.
#
############################################################################################

import  sys, os, string, subprocess, optparse

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class ConfigureVariables:
  # ... class containing all important variables, methods for submitting 
  #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  def __init__(self):
  #
  #              ... initialize the global set of variables ...
  #
    self.Bash_Configure_Arguments=""

    self.PythonMode=False


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def process_arguments():
#-----------------------------------------------------
#
#-----------------------------------------------------
  from optparse import OptionParser
  usage="usage: %prog arg1 arg2..."
  parser = OptionParser(usage)

  #  fill configure_arguments
  label_configure_arguments=" configure arguments (like '"' -nopt ....'"')"
  parser.add_option("--conf_args", type="string", dest="dest_configure_args",
                  help=label_configure_arguments, metavar="'"'arg1 arg2 ... '"' ")

  # is this pure python script configuration ?
  use_python_label='use python only based script, otherwise loading bash-based configure ' 
  parser.add_option("--config", action="store_true", default=GlobalVars.PythonMode,
                     dest="dest_prep_config", help=use_python_label)
 
  # assign ...
  (options, args) = parser.parse_args()

  # extract read arguments


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def call_bash_configure():
#-----------------------------------------------------
#   calls the external bash-based configure script
#-----------------------------------------------------
    bash_script="cd ..; ./configure "
    bash_script_arguments="--noopt"
    bash_script=bash_script_file+bash_script_arguments
    subprocess.call(bash_script, executable="bash", shell=True)


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
class PrintLevel:
  def __init__(self):
    self.iprint=0
  def set_iprint(self):
    self.iprint = iprint_read

#### check where configure.py sits - must be the main trunk
def check_position():
  print 'checking position of configure.......'

def process_configure_arguments():
  print 'processing configure arguments....'

def main(*args):
  global GlobalVars
  GlobalVars = ConfigureVariables() 
  print 'write me ! :)'

if __name__ == '__main__':
    call_bash_configure()  
    sys.exit(main(*sys.argv))

