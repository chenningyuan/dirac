#
#
# common functionalities used throughout the DIRAC python scripts 
#
#

def is_windows():
  # returns true if this machine is Windows operating system
  platform=distutils.util.get_platform()
  if platform.find("win") == -1 :
    return False
  else:
    return True


