#!/usr/bin/env python
#
# Script to identfy MPI type
# Written by Johann Pototschnig, spring 2021

# calls "<path>/mpiexec --version"
# Type is saved as MPI_TYPE variable.

import os, sys, subprocess, re

#! path to mpiexec
MPI_command = sys.argv[1]

def identify_mpi_type(compiler,version_flag):
    """
    identify MPI type
      commmand - full path to mpiexec
      version_flag - flag to print compiler version
    """
    
    try:
        output_compiler = subprocess.Popen([compiler,version_flag], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output_compiler, stderr = output_compiler.communicate()
    except:
        output_compiler = b""
        stderr=""
        
    if re.search('mpich', output_compiler.decode('utf-8'),re.IGNORECASE):
        mpitype='MPICH'
    elif re.search('hydra', output_compiler.decode('utf-8'),re.IGNORECASE):
        mpitype='MPICH'
    elif re.search('intel', output_compiler.decode('utf-8'),re.IGNORECASE):
        mpitype='INTEL'
    elif re.search('openmpi', output_compiler.decode('utf-8'),re.IGNORECASE):
        mpitype='OPENMPI'
    elif re.search('openrte', output_compiler.decode('utf-8'),re.IGNORECASE):
        mpitype='OPENMPI'
    elif re.search('open-mpi', output_compiler.decode('utf-8'),re.IGNORECASE):
        mpitype='OPENMPI'
    elif re.search('Spectrum', output_compiler.decode('utf-8'),re.IGNORECASE):
        #IBM Spectrum MPI
        mpitype='OPENMPI' 
#   elif re.search('Cray', output_compiler.decode('utf-8'),re.IGNORECASE):
#       #IBM Spectrum MPI
#       mpitype='MPICH'
    else:
        mpitype='None'

    # "--version" can fail (e.g. Debian Open MPI 5 cannot find its help files).
    # Then, look at the real launcher, e.g. /usr/bin/mpirun.openmpi, and for
    # ompi_info, which every Open MPI installs next to it.
    if mpitype == 'None':
        launcher = os.path.realpath(compiler)
        if re.search('openmpi|open-mpi', launcher, re.IGNORECASE) or \
           os.path.isfile(os.path.join(os.path.dirname(launcher), 'ompi_info')):
            mpitype='OPENMPI'

    return mpitype


#! looking for MPI
MPI_TYPE = identify_mpi_type(MPI_command,'--version')
#! this is the standard output for CMake
print(MPI_TYPE)

