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
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


def is_windows():
    """
    returns true if this machine is Windows operating system
    """
    import distutils
    platform = distutils.util.get_platform()
    return platform.find("win") == -1

def process_time_string_core(time_string):
    """
    core subroutine :
    process entering time_string of the form "#d#h#m#s" 
    or as pure integer number; 
    returns the resulting time in seconds (integer)
    """
    secs=0
    mins=0
    hours=0
    days=0
    # check if the time string is pure integer number - then return  the "virgin" value :)
    if time_string.isdigit(): return (int(time_string))

    # otherwise continue to cut the string into substrings according to the presence of time units: s,m,h,d
    d_pos=time_string.find("d"); 
    h_pos=time_string.find("h"); 
    m_pos=time_string.find("m"); 
    s_pos=time_string.find("s"); 
    #print 'd,h,m,s positions:',d_pos,is_d, h_pos,is_h, m_pos,is_m, s_pos, is_s

    positions=[s_pos,m_pos,h_pos,d_pos]; j=0
    num=[0,0,0,0]
    #print 'positions: ', positions
    for i in positions:
        if i >-1 :
            num_str=''
            #print 'active pos:',i
            i=i-1
            while (i>=0 and time_string[i] >= '0' and time_string[i] <= '9' ):
                #print  'number=',time_string[i]
                num_str+=str(time_string[i])
                i=i-1
            # revert string to get proper order of numbers
            num[j]=int(num_str[::-1])
            j=j+1
    # take care of time units
    k=0
    if s_pos > -1:
        secs=num[k]; k=k+1 
    if m_pos > -1:
        mins=num[k]; k=k+1 
    if h_pos > -1:
        hours=num[k]; k=k+1 
    if d_pos > -1:
        days=num[k]; k=k+1 

    return secs + (60*mins) + (60*60*hours) + (60*60*60*days)

def process_time_string(time_string):
    """
    process entering time_string of the form "#d#h#m#s" 
    or as pure integer number; 
    returns the resulting time in seconds (integer)
    """
    import sys
    try:
        return process_time_string_core(time_string)
    except:
        for stream in [sys.stdout,sys.stderr]:
            print >> stream,'error, wrong entering time string:',time_string
            print >> stream,'proper examples: 1h2m30s  2m40s 1m40s  1h20s 3h ...'
        sys.exit(13)

def get_username():
    """
    returns name of the user
    """
    import getpass
    return getpass.getuser()

def get_ip():
    """
    module for getting the LAN IP-address of the computer
    from  http://stackoverflow.com/questions/166506/finding-local-ip-addresses-in-python
    """
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
    """
    get the server's short hostname 
    NOTE: short hostnames are used in batch-system  "machinefiles" 
    only short master's hostname can be compared with slaves hostnames
    """
    import socket
    return socket.gethostname()

def get_full_hostname():
    """
    get the server's full hostname using its IP-address
    """
    import socket
    host=socket.gethostbyaddr(get_ip()) 
    return host[0]

def get_hostname_ip():
    """
    get the hostname's ip-adress; use own function get_ip
    """
    return get_ip()


def children_pid(pid):
    """
    get the list of the children pids of the entering pid (linux only)
    windows solution: http://mail.python.org/pipermail/python-win32/2005-March/002991.html
    http://web.archiveorange.com/archive/v/MFLkXe0JXTy7BRtHxfbo    
    see also:
    http://www.damonkohler.com/2008/03/how-to-check-if-process-is-still.html
    http://stackoverflow.com/questions/4669754/python-kill-all-processes-owned-by-user
    """
    import subprocess as SP
    proc = SP.Popen('ps -o pid,ppid ax | grep "%d"' % pid, shell=True, stdout=SP.PIPE)
    pidppid = [x.split() for x in proc.communicate()[0].split("\n") if x]
    return list(int(p) for p, pp in pidppid if int(pp) == pid)

def get_pid_label(input_pid):
    """
    get process identification based on its input_pid
    TODO: adapt for Windows
    """
    import subprocess
    # Linux solution
    command="ps -l -F -p "+str(input_pid)
    return subprocess.Popen(command,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE).communicate()[0]

def get_all_children(input_pid):
    """
    get PIDs of all childrens processes of the process "input_pid"
    (http://www.linuxquestions.org/questions/linux-general-1/trying-to-kill-a-process-and-all-of-its-children-and-grandchildren-how-680123/)
    TODO: adapt for Windows
    http://mail.python.org/pipermail/python-win32/2005-February/002943.html
    """
    import os, sys
    # get the pid, pgid, ppid of our current processes:
    command = "ps eo pid,pgid,ppid"
    psraw = os.popen(command).readlines()
    psList = []
    childrenList = []
    for ps in psraw[1:]: # 1: gets rid of header
        psList.append(map(int,ps.split()))
    pgid = 0

    # find the pgid (parent) of the spawned subprocess:
    for ps in psList:
        if int(input_pid) in ps:
            pgid = ps[1]
            break
    if pgid == 0:
        print >>sys.stderr, "Something screwed up trying to find pids. fudge."
        return

    # get a list of all pids in the pgid except the group owner:
    for ps in psList:
        if pgid in ps and pgid != ps[0]:  # check [0] so we don't kill ourselves
            childrenList.append(ps[0])

    # return the list
    return childrenList

class disk_space:
    """
    extracts disk space data - capacity, available, used -  in GB
    TODO: Adapt for Windows
    """
    def __init__(self,path):
        """
        fill (initialize) wanted variables, no other info needed sofar
        TODO: check the path
        """
        import os
        gb_factor = float(1024*1024*1024)
        self.disk_info = os.statvfs(path)
        self.capacity = ( self.disk_info.f_bsize * self.disk_info.f_blocks ) / gb_factor
        self.available = ( self.disk_info.f_bsize * self.disk_info.f_bavail ) /  gb_factor
        self.used = ( self.disk_info.f_bsize * (self.disk_info.f_blocks - self.disk_info.f_bavail)) / gb_factor

class extract_path_file:
    """
    general class which can separate the pure file name from its path
    used (not only) for the backup procedure
    """
    def __init__(self):
        self.pure_file_name = ''
        self.pure_file_path = ''
        self.short_file_path= '' 
        self.extracted = False
    def get_path_file(self,full_file_name):
        """
        extracts file name and it's path (directory) by splitting string
        """
        from string import rfind
        if is_windows():
            string_to_find='\\' # Windows
        else: 
            string_to_find='/' # Unix
        position = full_file_name.rfind(string_to_find)
        if (position+1 < len(full_file_name)) :
            self.pure_file_name=full_file_name[position+1:len(full_file_name)]
            self.pure_file_path=full_file_name[0:position]
            self.extracted=True
        else: 
            print 'can not extract' 
            self.extracted=False

        if (len(self.pure_file_path) > 20):
           max_len=12
           temp_path=self.pure_file_path
           len_tp=len(temp_path)
           #print  temp_path 
           self.short_file_path=temp_path[0:max_len]+'...'+temp_path[len_tp-max_len:len_tp]
 
def backup_file(some_file,upto,iprint):
    """
    Important routine doing backup of the file "some_file", version up to
    number "upto"
    """
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
        if (result.extracted is True and result_p1.extracted is True  and \
            result.pure_file_path == result_p1.pure_file_path ) :
            if (iprint >=3):
                print 'backing up file:', result.pure_file_name,'-->',result_p1.pure_file_name,\
            'in:',result.short_file_path
            pass
        else :
            if (iprint >=3):
                print 'backing up file:',backup_file_name,'-->',backup_file_name_p1
            pass

#   backup-s file by copying it up to 'upto' version
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

    epf=extract_path_file()
    epf.get_path_file(some_file)
    if (iprint >= 0):
        #print 'backed up file ',some_file
        print '  backed up file ',epf.pure_file_name

def make_old_file(path):
    """
    Create (if not already there) a file at path,
    set the date ("time stamp") in the past.
    """
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

class svn_info:
    """
    handle svn repository data
    """
    def __init__(self):
        """
        initialize  svn_revision, svn_url
        """
        import os, subprocess
        command="svn info" # valid for both Unix,Windows
        try:
            process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
            output = process.stdout.read().splitlines(False)
            self.svn_revision=output[4].replace("Revision:","")
            self.svn_revision=self.svn_revision.replace(" ","")
            self.svn_url=output[2].replace("Repository Root: ","")
        except:
#           this is not a subversion working copy
#           (some people are working with git)
            self.svn_revision = 'unknown; not a subversion working copy'
            self.svn_url      = 'unknown; not a subversion working copy'
    def get_svn_revision(self):
        return self.svn_revision
    def get_svn_url(self):
        return self.svn_url

def create_link(source, dest_dir, dest_file):
    """
    create symbolic link (Unix), before make a comprehensive check
    """
    import os
    if not os.path.exists(source):
        print "create_link error ! source does not exist:"+source
        sys.exit(source+' does not exist ! No linking possible, error exit !')
    
    if not os.path.exists(dest_dir):
        print "create_link error ! destination dir does not exist :"+dest_dir
    
    dest = dest_dir+"/"+dest_file
    
    # always remove old stuff (file or link) before
    if os.path.isfile(dest) or os.path.islink(dest):
        os.remove(dest)
    
    if not is_windows():
        try:
            os.symlink(source,dest)
        except (OSError, ValueError):
            print "*** linking error !!! "
            print "source : ",source
            print "destination : ",dest
            sys.exit("Linking error ! Check source and destination !")

def copy_parent(src,dest):
    """
    ensures copying with parent, ie Linux 'cp --parent' functionality
    this function serves for 'Windows'-ing DIRAC
    TODO: write me :)
    """
    import os, shutil
    command = 'cp --parent %s %s' % (src, dest)
    os.system(command)

def vacuum_cleaner(where,what,recursive=False):
    """
    remove all files specified as 'what' from the 'where' directory
    examples of 'what': '*.o','*.r','*.opt','.svn','*~'...
    written as the replacement of: 
             os.system("rm -f  `find %s -name '*.o'`"   % scratch_dir)
    TODO: implement for Windows as pure python function
    """
    import os, shutil
    find_command="`find %s -name '"+what +"'`"
    if (recursive):
        rm_command="rm -rf  "+find_command
    else:
        rm_command="rm -f "+find_command
    #print rm_command
    os.system( rm_command % where )

def find_executable(executable, path=None):
    """Try to find 'executable' in the directories listed in 'path' (a
    string listing directories separated by 'os.pathsep'; defaults to
    os.environ['PATH']).  Returns the complete filename or None if not
    found
    from: http://snippets.dzone.com/posts/show/6313
    """
    import os
    import sys
    if path is None:
        path = os.environ['PATH']
    paths = path.split(os.pathsep)
    extlist = ['']
    if os.name == 'os2':
        (base, ext) = os.path.splitext(executable)
        # executable files on OS/2 can have an arbitrary extension, but
        # .exe is automatically appended if no dot is present in the name
        if not ext:
            executable = executable + ".exe"
    elif sys.platform == 'win32':
        pathext = os.environ['PATHEXT'].lower().split(os.pathsep)
        (base, ext) = os.path.splitext(executable)
        if ext.lower() not in pathext:
            extlist = pathext
    for ext in extlist:
        execname = executable + ext
        if os.path.isfile(execname):
            return execname
        else:
            for p in paths:
                f = os.path.join(p, execname)
                if os.path.isfile(f):
                    return f
    else:
        return None

