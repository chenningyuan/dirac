#! /usr/bin/env python        
#
                             
import  sys, os, string, subprocess, optparse
                                      
############################################################################################
#
# Hoping for someone (a student ?) to turn the bash-based configure into the python equivalent.
#
############################################################################################

def bash_configure():
    subprocess.call("./configure", executable="bash", shell=True)


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
  print 'write me ! :)'


if __name__ == '__main__':
    bash_configure()  
    sys.exit(main(*sys.argv))

