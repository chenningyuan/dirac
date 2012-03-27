#! /usr/bin/env python
#
# Purpose: wrapper for launching lower Dirac scripts (pam, wrapper.py, runtest) in Grid environment. 
#
# Inspired by the simplicity of the Ulf's script "wrapper.py".
#
# Written by Miro ILIAS, Matej Bel University,
#        Banska Bystrica, Slovakia, 2012
#
# Contributors:  
#
#
##############################################################################################################################

import os, subprocess, shutil, glob, signal, optparse, sys, shlex, tarfile, tempfile, zipfile, socket

if sys.version < '2.4':
    sys.exit("lauch.py: wrong python version !!!")

#
# some global variables
#
vo=None 

def process_arguments():

    usage = "launch.py --vo=<your virt.organization>"
    parser = optparse.OptionParser(usage)

    group = optparse.OptionGroup(parser, 'Input files - mandatory specifications')
    group.add_option('--vo', type='string', action='store',
                         dest='vo_string',
                         help='virtual organization')
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    # assign VO
    if options.vo_string is not None:
        vo=options.vo_string
        print 'vo=',vo

def main(*args):
    print "Welcome, here python launcher for grid computing."
    process_arguments()


if __name__ == '__main__':
    sys.exit(main(*sys.argv))


