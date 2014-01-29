#!/usr/bin/env python

import sys, subprocess, re

#! GNU, Intel, PGI, XL
C_COMPILER_VENDOR = sys.argv[1]
CXX_COMPILER_VENDOR = sys.argv[2]
Fortran_COMPILER_VENDOR = sys.argv[3]

#! looking for C compiler version
if C_COMPILER_VENDOR == 'GNU':
    try:
        outputGCC = subprocess.check_output(['gcc', '--version'])
    except:
        outputGCC = ''
    finally:
        matchGCC = re.search(r'\d+\.\d+\.\d+', outputGCC)
        if matchGCC:
            GCCversion = matchGCC.group()
        else:
            GCCversion = '0.0.0'

    C_COMPILER_VERSION = GCCversion
elif C_COMPILER_VENDOR == 'Intel':
    try:
#! Intel uses stderr
        outputIntelC = subprocess.check_output(['icc', '-V'], stderr=subprocess.STDOUT)
    except:
        outputIntelC = ''
    finally:
        matchIntelC = re.search(r'\d+\.\d+\.\d+\.\d+', outputIntelC)
        if matchIntelC:
            IntelCversion = matchIntelC.group()
        else:
            IntelCversion = '0.0.0.0'
			
    C_COMPILER_VERSION = IntelCversion
elif C_COMPILER_VENDOR == 'PGI':
    try:
        outputPGIC = subprocess.check_output(['pgcc', '-V'])
    except:
        outputPGIC = ''
    finally:
        matchPGIC = re.search(r'\d+\.\d+', outputPGIC)
        if matchPGIC:
            PGICversion = matchPGIC.group()
        else:
            PGICversion = '0.0'
			
    C_COMPILER_VERSION = PGICversion
elif C_COMPILER_VENDOR == 'XL':
    #! not implemented yet
    C_COMPILER_VERSION = '0.0.0'
else:
    #! unknown, not supported vendor
    C_COMPILER_VERSION = '0.0.0'

#! looking for C++ compiler version
if CXX_COMPILER_VENDOR == 'GNU':
    try:
        outputGXX = subprocess.check_output(['g++', '--version'])
    except:
        outputGXX = ''
    finally:
        matchGXX = re.search(r'\d+\.\d+\.\d+', outputGXX)
        if matchGXX:
            GXXversion = matchGXX.group()
        else:
            GXXversion = '0.0.0'

    CXX_COMPILER_VERSION = GXXversion
elif CXX_COMPILER_VENDOR == 'Intel':
    try:
#! Intel uses stderr
        outputIntelCXX = subprocess.check_output(['icpc', '-V'], stderr=subprocess.STDOUT)
    except:
        outputIntelCXX = ''
    finally:
        matchIntelCXX = re.search(r'\d+\.\d+\.\d+\.\d+', outputIntelCXX)
        if matchIntelCXX:
            IntelCXXversion = matchIntelCXX.group()
        else:
            IntelCXXversion = '0.0.0.0'
			
    CXX_COMPILER_VERSION = IntelCXXversion
elif CXX_COMPILER_VENDOR == 'PGI':
    try:
        outputPGICXX = subprocess.check_output(['pgCC', '-V'])
    except:
        outputPGICXX = ''
    finally:
        matchPGICXX = re.search(r'\d+\.\d+', outputPGICXX)
        if matchPGICXX:
            PGICXXversion = matchPGICXX.group()
        else:
            PGICXXversion = '0.0'
			
    CXX_COMPILER_VERSION = PGICXXversion
elif CXX_COMPILER_VENDOR == 'XL':
    #! not implemented yet
    CXX_COMPILER_VERSION = '0.0.0'
else:
    #! unknown, not supported vendor
    CXX_COMPILER_VERSION = '0.0.0'

#! looking for Fortran compiler version
if Fortran_COMPILER_VENDOR == 'GNU':
    try:
        outputGFORTRAN = subprocess.check_output(['gfortran', '--version'])
    except:
        outputGFORTRAN = ''
    finally:
        matchGFORTRAN = re.search(r'\d+\.\d+\.\d+', outputGFORTRAN)
        if matchGFORTRAN:
            GFORTRANversion = matchGFORTRAN.group()
        else:
            GFORTRANversion = '0.0.0'

    Fortran_COMPILER_VERSION = GFORTRANversion
elif Fortran_COMPILER_VENDOR == 'Intel':
    try:
#! Intel uses stderr
        outputIntelFORTRAN = subprocess.check_output(['ifort', '-V'], stderr=subprocess.STDOUT)
    except:
        outputIntelFORTRAN = ''
    finally:
        matchIntelFORTRAN = re.search(r'\d+\.\d+\.\d+\.\d+', outputIntelFORTRAN)
        if matchIntelFORTRAN:
            IntelFORTRANversion = matchIntelFORTRAN.group()
        else:
            IntelFORTRANversion = '0.0.0.0'
			
    Fortran_COMPILER_VERSION = IntelFORTRANversion
elif Fortran_COMPILER_VENDOR == 'PGI':
    try:
        outputPGIFORTRAN = subprocess.check_output(['pgf90', '-V'])
    except:
        outputPGIFORTRAN = ''
    finally:
        matchPGIFORTRAN = re.search(r'\d+\.\d+', outputPGIFORTRAN)
        if matchPGIFORTRAN:
            PGIFORTRANversion = matchPGIFORTRAN.group()
        else:
            PGIFORTRANversion = '0.0'
			
    Fortran_COMPILER_VERSION = PGIFORTRANversion
elif Fortran_COMPILER_VENDOR == 'XL':
    #! not implemented yet
    Fortran_COMPILER_VERSION = '0.0.0'
else:
    #! unknown, not supported vendor
    Fortran_COMPILER_VERSION = '0.0.0'


#! this is output for CMake (must be in this order and with ;)
print C_COMPILER_VERSION + ';' + CXX_COMPILER_VERSION + ';' + Fortran_COMPILER_VERSION
