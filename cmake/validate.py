def contains_flag(sys_argv, flag):
    return (any(x for x in sys_argv if x.startswith('--{0}='.format(flag))))


def validate_setup_args(sys_argv, arguments):

    # if --mpi is selected and compilers are not selected
    # then compilers default to mpif90, mpicc, and mpicxx
    if arguments['--mpi']:
        if not contains_flag(sys_argv, 'fc') and not contains_flag(sys_argv, 'cc') and not contains_flag(sys_argv, 'cxx'):
           arguments['--fc'] = 'mpif90'
           arguments['--cc'] = 'mpicc'
           arguments['--cxx'] = 'mpicxx'

    return arguments
