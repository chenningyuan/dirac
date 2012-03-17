#! /usr/bin/env python
#
# Purpose: wrapper for launching lower Dirac scripts (pam, wrapper.py, runtest) in Grid environment. 
#
# Written by Miro ILIAS, Matej Bel University,
#        Banska Bystrica, Slovakia, 2012
#
# Contributors:  
#
#
##############################################################################################################################

import sys

if sys.version < '2.4':
    sys.exit("bad python version !!!")





def main(*args):
    print "hello python launcher in grid !!"



if __name__ == '__main__':
    sys.exit(main(*sys.argv))


