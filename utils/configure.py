#! /usr/bin/env python        
#
                             
import  sys, os, string, optparse
                                      
############################################################################################
#
# Hoping for someone (a student ?) to turn the bash-based configure into the python equivalent.
#
############################################################################################

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
    sys.exit(main(*sys.argv))

