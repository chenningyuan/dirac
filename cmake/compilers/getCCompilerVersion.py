#!/usr/bin/env python

import sys, subprocess, re

#! GNU, Intel, PGI, XL
COMPILER_VENDOR = sys.argv[1]

if COMPILER_VENDOR == 'GNU':
    try:
        outputGCC = subprocess.check_output(['gcc', '--version'])
    except:
        outputGCC = ''
    finally:
        matchGCC = re.search(r'\d\.\d\.\d', outputGCC)
        if matchGCC:
            GCCversion = matchGCC.group()
        else:
            GCCversion = '0.0.0'

    print GCCversion
elif COMPILER_VENDOR == 'Intel':
    #! not implemented yet
    print '0.0.0'
elif COMPILER_VENDOR == 'PGI':
    #! not implemented yet
    print '0.0.0'
elif COMPILER_VENDOR == 'XL':
    #! not implemented yet
    print '0.0.0'
else:
	#! unknown, not supported vendor
    print '0.0.0'
