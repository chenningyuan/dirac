#########################################################################
#
# common functionalities used throughout the DIRAC python scripts 
#
# Please DIRAC user add more common routines/classes here. 
#
#########################################################################
#
# to utilize them in your python scripts use this command once (Unix):
#
#  sys.path.append("$DIRAC_dir/utils")
#
# and in each routine or at the beginning of the file (as globally visible):
#
#  from common import get_username,get_hostname
#  from common import *
#
#########################################################################
#
# Written by Miro ILIAS, October 2010, Banska Bystrica, Slovakia
#            .......  
#
#########################################################################
import os

def is_windows():
  import distutils
  # returns true if this machine is Windows operating system
  platform=distutils.util.get_platform()
  if platform.find("win") == -1 :
    return False
  else:
    return True

def get_username():
  import getpass
  if is_windows():
    return win32api.GetUserName()
  else:
    return getpass.getuser()

def get_ip():
#########################################################################################
#
# module for getting the lan ip address of the computer
#
# from  http://stackoverflow.com/questions/166506/finding-local-ip-addresses-in-python
#########################################################################################
  import os
  import socket
  if os.name != "nt":
    import fcntl
    import struct
    def get_interface_ip(ifname):
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        return socket.inet_ntoa(fcntl.ioctl(
                        s.fileno(),
                        0x8915,  # SIOCGIFADDR
                        struct.pack('256s', ifname[:15])
                )[20:24])
    ip = socket.gethostbyname(socket.gethostname())
    if ip.startswith("127.") and os.name != "nt":
        interfaces = ["eth0","eth1","eth2","wlan0","wlan1","wifi0","ath0","ath1","ppp0"]
        for ifname in interfaces:
                try:
                        ip = get_interface_ip(ifname)
                        break;
                except IOError:
                        pass
    return ip

def get_hostname():
###  get the hostname with ip-adress ###
  import socket
  this_ip = get_ip()
  this_host_id =  socket.gethostname()+" (IP="+this_ip+")"
  return this_host_id

def backup_file(some_file,upto):
  import os.path, shutil
#
# backup-s file by copiyng it  'upto' version
#
  #print 'backing up the file ',some_file
  upto_default = 7
  if (upto is None):
    upto = upto_default
  #for vers in range(upto-1,0,-1):
  for vers in range(upto-1,-1,-1):
    backup_file_name=some_file+"."+str(vers)
    backup_file_name_p1=some_file+"."+str(vers+1)
    #print backup_file_name 
    if (os.path.isfile(backup_file_name)):
      shutil.copy(backup_file_name,backup_file_name_p1)
      print 'backing up file:',backup_file_name,'-->',backup_file_name_p1
    # copy the existing out into out.0
  if (os.path.isfile(some_file)):
    shutil.copy(some_file,backup_file_name)
    print 'backing up file:',some_file,'-->',backup_file_name

