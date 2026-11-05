import os
import shutil
import sys


def contains_flag(sys_argv, flag):
    return (any(x for x in sys_argv if x.startswith('--{0}='.format(flag))))


# (setup flag, environment variable honoured by CMake)
_COMPILER_FLAGS = [('fc', 'FC'), ('cc', 'CC'), ('cxx', 'CXX')]


def _compiler_sets(mpi):
    """
    Candidate compiler sets (Fortran, C, C++), in order of preference.
    The first set with all three compilers found in PATH is used.
    """
    # Cray compiler wrappers: they also provide MPI, and must come before
    # mpif90/mpicc/mpicxx, which on Cray bypass the programming environment
    sets = [('ftn', 'cc', 'CC')]
    if mpi:
        sets.append(('mpif90', 'mpicc', 'mpicxx'))
    sets.append(('gfortran', 'gcc', 'g++'))
    return sets


def postprocess_args(sys_argv, arguments):

    # Compilers are selected as follows, for each language separately:
    #   1. the compiler given explicitly with --fc/--cc/--cxx;
    #   2. if the environment variable FC/CC/CXX is set, nothing is passed
    #      and CMake uses the environment variable;
    #   3. otherwise the compiler from the first set in _compiler_sets()
    #      with all three compilers found in PATH;
    #   4. if no such set exists, nothing is passed and CMake detects the
    #      compiler itself.
    # We do not rely on CMake detection by default since CMake e.g. prefers
    # f95 to ftn and to gfortran (for f95 see CMakeDetermineFortranCompiler.cmake).
    chosen_set = None
    for compilers in _compiler_sets(arguments['--mpi']):
        if all(shutil.which(c) for c in compilers):
            chosen_set = compilers
            break

    for i, (flag, env) in enumerate(_COMPILER_FLAGS):
        if contains_flag(sys_argv, flag):
            continue
        if os.environ.get(env):
            arguments['--{0}'.format(flag)] = None
        elif chosen_set is not None:
            arguments['--{0}'.format(flag)] = chosen_set[i]

    # if one of the compilers contains "mpi" and --mpi is not selected, it is probably a user error
    # in this case stop the configuration
    asking_for_mpi_compiler = False
    for flag, _ in _COMPILER_FLAGS:
        if contains_flag(sys_argv, flag):
            if 'mpi' in arguments['--{0}'.format(flag)]:
                asking_for_mpi_compiler = True
    if asking_for_mpi_compiler and not arguments['--mpi']:
        sys.stderr.write('ERROR: you ask for an MPI compiler but have not specified --mpi\n')
        sys.exit(1)

    return arguments
