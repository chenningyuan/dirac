#
#
# common functionalities used throughout the DIRAC python scripts 
#
#
import os

def is_windows():
  # returns true if this machine is Windows operating system
  platform=distutils.util.get_platform()
  if platform.find("win") == -1 :
    return False
  else:
    return True

def get_username():
  # get the user name
  return  'XY'

def get_hostname():
  # get the hostname
  return  'AB'



