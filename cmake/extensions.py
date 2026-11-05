import shutil
import sys


def contains_flag(sys_argv, flag):
    return (any(x for x in sys_argv if x.startswith('--{0}='.format(flag))))


def postprocess_args(sys_argv, arguments):

    # Compilers that are not given explicitly are left unset (None), in which
    # case setup does not pass CMAKE_<LANG>_COMPILER and CMake detects the
    # compilers itself (honouring the FC, CC and CXX environment variables).
    #
    # If --mpi is selected, compilers that are not given explicitly default to
    # the MPI wrappers mpif90, mpicc and mpicxx, but only if these are found
    # in PATH. On systems where MPI is provided through other wrappers
    # (e.g. ftn/cc/CC on Cray), CMake detection is used instead.
    if arguments['--mpi']:
        for flag, wrapper in [('fc', 'mpif90'), ('cc', 'mpicc'), ('cxx', 'mpicxx')]:
            if not contains_flag(sys_argv, flag) and shutil.which(wrapper):
                arguments['--{0}'.format(flag)] = wrapper

    # if one of the compilers contains "mpi" and --mpi is not selected, it is probably a user error
    # in this case stop the configuration
    asking_for_mpi_compiler = False
    for flag in ['fc', 'cc', 'cxx']:
        if contains_flag(sys_argv, flag):
            if 'mpi' in arguments['--{0}'.format(flag)]:
                asking_for_mpi_compiler = True
    if asking_for_mpi_compiler and not arguments['--mpi']:
        sys.stderr.write('ERROR: you ask for an MPI compiler but have not specified --mpi\n')
        sys.exit(1)

    return arguments
