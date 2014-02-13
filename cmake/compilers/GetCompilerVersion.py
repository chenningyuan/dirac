#!/usr/bin/env python
#
# Script to determine and verify compilers (C, C++ and Fortran) versions.
#
# For each vendor/system call compiler specific subprocess command (like "gcc --version")
# returning desired compiler version.
#
# This script returns (as stdout) for CMake a list of 3 strings - compilers versions - separated by semicommas.
#
# Written by Ivan Hrasko, february 2014
#
#


import sys, subprocess, re, platform

#! GNU, Intel, PGI, XL
C_COMPILER_VENDOR = sys.argv[1]
CXX_COMPILER_VENDOR = sys.argv[2]
Fortran_COMPILER_VENDOR = sys.argv[3]
#! operating system name
CMAKE_SYSTEM_NAME = sys.argv[4]

#! in every call of subprocess.check_output() stderr is merged into stdout
#! by option stderr=subprocess.STDOUT
#! this is really needed probably only by Intel, but it probably also makes no problem
#! when used everywhere

#! looking for C compiler version
if C_COMPILER_VENDOR == 'GNU':
    try:
        outputGCC = subprocess.check_output(['gcc', '--version'], stderr=subprocess.STDOUT)
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
        if CMAKE_SYSTEM_NAME == 'Windows':
            outputIntelC = subprocess.check_output(['icl', '/V'], stderr=subprocess.STDOUT)
        else:
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
        outputPGIC = subprocess.check_output(['pgcc', '-V'], stderr=subprocess.STDOUT)
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
    try:
        outputXLC = subprocess.check_output(['xlc', '-qversion'], stderr=subprocess.STDOUT)
    except:
        outputXLC = ''
    finally:
        matchXLC = re.search(r'\d+\.\d+', outputXLC)
        if matchXLC:
            XLCversion = matchXLC.group()
        else:
            XLCversion = '0.0'
			
    C_COMPILER_VERSION = XLCversion
else:
    #! unknown, not supported vendor
    C_COMPILER_VERSION = '0.0'

#! looking for C++ compiler version
if CXX_COMPILER_VENDOR == 'GNU':
    try:
        outputGXX = subprocess.check_output(['g++', '--version'], stderr=subprocess.STDOUT)
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
        if CMAKE_SYSTEM_NAME == 'Windows':
            outputIntelCXX = subprocess.check_output(['icl', '/V'], stderr=subprocess.STDOUT)
        else:
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
        outputPGICXX = subprocess.check_output(['pgCC', '-V'], stderr=subprocess.STDOUT)
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
    try:
        outputXLCXX = subprocess.check_output(['xlC', '-qversion'], stderr=subprocess.STDOUT)
    except:
        outputXLCXX = ''
    finally:
        matchXLCXX = re.search(r'\d+\.\d+', outputXLCXX)
        if matchXLCXX:
            XLCXXversion = matchXLCXX.group()
        else:
            XLCXXversion = '0.0'
			
    CXX_COMPILER_VERSION = XLCXXversion
else:
    #! unknown, not supported vendor
    CXX_COMPILER_VERSION = '0.0'

#! looking for Fortran compiler version
if Fortran_COMPILER_VENDOR == 'GNU':
    try:
        outputGFORTRAN = subprocess.check_output(['gfortran', '--version'], stderr=subprocess.STDOUT)
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
        if CMAKE_SYSTEM_NAME == 'Windows':
            outputIntelFORTRAN = subprocess.check_output(['ifort', '/V'], stderr=subprocess.STDOUT)
        else:
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
        outputPGIFORTRAN = subprocess.check_output(['pgf90', '-V'], stderr=subprocess.STDOUT)
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
    try:
        outputXLFORTRAN = subprocess.check_output(['xlf90', '-qversion'], stderr=subprocess.STDOUT)
    except:
        outputXLFORTRAN = ''
    finally:
        matchXLFORTRAN = re.search(r'\d+\.\d+', outputXLFORTRAN)
        if matchXLFORTRAN:
            XLFORTRANversion = matchXLFORTRAN.group()
        else:
            XLFORTRANversion = '0.0'
			
    Fortran_COMPILER_VERSION = XLFORTRANversion
else:
    #! unknown, not supported vendor
    Fortran_COMPILER_VERSION = '0.0'


#! this is the standard output for CMake (must be in this order and with ;)
print C_COMPILER_VERSION + ';' + CXX_COMPILER_VERSION + ';' + Fortran_COMPILER_VERSION
