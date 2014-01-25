#!/usr/bin/env python

import sys, subprocess, re

#! GNU, Intel, PGI, XL
COMPILER_VENDOR = sys.argv[1]

if COMPILER_VENDOR == 'GNU':
    try:
        outputGFORTRAN = subprocess.check_output(['gfortran', '--version'])
    except:
        outputGFORTRAN = ''
    finally:
        matchGFORTRAN = re.search(r'\d\.\d\.\d', outputGFORTRAN)
        if matchGFORTRAN:
            GFORTRANversion = matchGFORTRAN.group()
        else:
            GFORTRANversion = '0.0.0'

    print GFORTRANversion
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
