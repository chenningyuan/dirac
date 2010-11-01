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
# and to use the external module type in each routine or at the beginning of the file (as globally visible):
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

#no global 'import XXX' here, these are within individual fuctions
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

class extract_path_file:
#
# used for backup procedure
#
  def __init__(self):
    self.pure_file_name=''
    self.pure_file_path=''
    self.extracted=False
  def get_path_file(self,full_file_name):
# extracts file name and it's path (directory) by splitting string
    from string import rfind
    #print 'extract. path, file..',full_file_name
    if is_windows():
      string_to_find='\\' # Windows
    else: 
      string_to_find='/' # Unix
    position=full_file_name.rfind(string_to_find)
    if (position+1 < len(full_file_name)) :
      self.pure_file_name=full_file_name[position+1:len(full_file_name)]
      self.pure_file_path=full_file_name[0:position]
      self.extracted=True
      #print 'was extracted'
      #print pure_file_path,pure_file_name
      #return sepure_file_path,pure_file_name
    else: 
      print 'can not extract' 
      self.extracted=False

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def backup_file(some_file,upto):
#############################################################################
#
# Important routine doing backup of the file "some_file", version up to
# numebr "upto"
#
#
#############################################################################
  import os.path, shutil
  def do_copy_print(backup_file_name,backup_file_name_p1):
      shutil.copy(backup_file_name,backup_file_name_p1)
      #         "file.n"  -->  "file.n+1"
      pure_file_path='';pure_file_name=''
      pure_file_path_p1='';pure_file_name_p1=''
      result=extract_path_file()
      result.get_path_file(backup_file_name)
      result_p1=extract_path_file()
      result_p1.get_path_file(backup_file_name_p1)
      #print result.extracted, result.pure_file_path
      #print result_p1.extracted, result_p1.pure_file_path 
      if (result.extracted is True and result_p1.extracted is True  and \
          result.pure_file_path == result_p1.pure_file_path ) :
        print 'backing up file:', result.pure_file_name,'-->',result_p1.pure_file_name,\
        'in directory:',result.pure_file_path
      else:
        print 'backing up file:',backup_file_name,'-->',backup_file_name_p1

#
# backup-s file by copying it up to 'upto' version
#
  upto_default = 7
  if (upto is None):
    upto = upto_default
  for vers in range(upto-1,-1,-1):
    backup_file_name=some_file+"."+str(vers)
    backup_file_name_p1=some_file+"."+str(vers+1)
    if (os.path.isfile(backup_file_name)):
      do_copy_print(backup_file_name,backup_file_name_p1)
  if (os.path.isfile(some_file)):
    do_copy_print(some_file,backup_file_name)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def make_old_file(path):
#############################################################
#  Create (if not already there) a file at path,
#  set the date ("time stamp") in the past.
#############################################################
  import os, datetime, time
  if not os.path.isfile(path):
    f = open(path,'w')
    f.close
  today = datetime.datetime.now()
  pastday = today - datetime.timedelta(days=11)
  atime = int(time.mktime(pastday.timetuple()))
  mtime = atime
  times = (atime,mtime)
  os.utime(path,times)

def get_svn_revision():
# returns svn revision
  import os, subprocess
  command="svn info"
  process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
  #svn_revision_raw = os.popen("svn info | grep 'Revision' | awk '{ print $2 }'", "r").read()
  output = process.stdout.read().splitlines(False)
  #print output
  svn_revision=output[4].replace("Revision:","")
  svn_revision=svn_revision.replace(" ","")
  #print svn_revision
  return svn_revision

def get_svn_url():
# returns svn url, TODO: extend for Windows - avoid Linux commands
  import os, subprocess
  command="svn info"
  process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
  output = process.stdout.read().splitlines(False)
  #print output
  #svn_url = os.popen("echo `svn info | grep 'URL' | awk '{ print $2 }'`", "r").read()
  svn_url=output[2].replace("Repository Root: ","")
  return svn_url

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def create_link(source, dest_dir, dest_file):
###################################################################
#
#  create symbolic link (Unix), before make a comprehensive check
#
###################################################################
  import os
  if not os.path.exists(source):
     print "create_link error ! source does not exist:"+source
     sys.exit(source+' does not exist ! No linking possible, error exit !')

  if not os.path.exists(dest_dir):
     print "create_link error ! destination dir does not exist :"+dest_dir

  dest = dest_dir+"/"+dest_file

  # always remove old stuff (file or link) before
  if os.path.isfile(dest) or os.path.islink(dest):
      #print 'removing old regular file:'+dest
      os.remove(dest)

  if not is_windows():
    try:
       os.symlink(source,dest)
    except (OSError, ValueError):
       print "*** linking error !!! "
       print "source : ",source
       print "destination : ",dest
       sys.exit("Linking error ! Check source and destination !")


