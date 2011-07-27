#! /usr/bin/env python
#
import sys
import os
import string
import glob
import shlex
import shutil
def handle_diracrc():
    install_dir = '${CMAKE_SOURCE_DIR}'
    fullargs = []
    rc_file_found=None
    ret_var=None
# hierarchy - from closest to the farthest
    paths = ['.diracrc','diracrc',os.path.join(install_dir,'.diracrc'),
          os.path.join(install_dir,'diracrc'),
          os.path.join(os.path.expanduser('~'),'.diracrc')]
    for path in paths:
        try:
            f = open(path,'r')
            read_args=False
            rc_file_found=path
            for l in f.readlines():
                l = l.strip()
                if len(l) > 0:
                    if l[0] != '#':
                        for a in shlex.split(l):
                            fullargs.append(os.path.expandvars(a))
                            read_args=True
            if read_args:
                break
        except Exception:
            pass
    if rc_file_found:
        print '-- cmake python script found DIRAC configuration file:',rc_file_found
    if fullargs: # save info about the diracrc file
        print '-- DIRAC config file contains these flags for pam:',fullargs
    else:
        if rc_file_found:
            print '-- but this file contains no active flags for pam - fix it !'
        else:
            dest='${CMAKE_SOURCE_DIR}/cmake/diracrc.in'
            src='${CMAKE_SOURCE_DIR}/diracrc'
            print '-- no dirac cofiguration file found - created default template:',src
            shutil.copy(dest,src)
# launch
handle_diracrc()
