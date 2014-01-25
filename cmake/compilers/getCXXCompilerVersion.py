#!/usr/bin/env python

import sys, subprocess, re

#! GNU, Intel, PGI, XL
COMPILER_VENDOR = sys.argv[1]

if COMPILER_VENDOR == 'GNU':
    try:
        outputG++ = subprocess.check_output(['g++', '--version'])
    except:
        outputG++ = ''
    finally:
        matchG++ = re.search(r'\d\.\d\.\d', outputG++)
        if matchG++:
            G++version = matchG++.group()
        else:
            G++version = '0.0.0'

    print G++version
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
