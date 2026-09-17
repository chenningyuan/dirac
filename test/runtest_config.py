def configure(options, input_files, extra_args):
    """
    This function is used by runtest to configure runtest
    at runtime for code specific launch command and file naming.
    """

    from os import path
    from sys import platform

    (program, mol, inp) = input_files

    # Three shapes are recognised here:
    #   'a:b'    atomic IAO fragments grouped by a labelled xyz file (a = geometry,
    #            b = labelled file) -> groupfrag_*.py
    #   '[a,b,]' a comma-separated list of molecular fragments      -> molfrag_*.py
    #   'a'      a single molecule, atomic fragments                -> atomfrag_*.py
    if ':' in mol:
        launcher = 'groupfrag_{}.py'.format(program)
        output_prefix = '{0}_{1}'.format(program, mol.split(':')[0])
    elif ',' in mol:
        launcher = 'molfrag_{}.py'.format(program)
        output_prefix = '{0}_{1}'.format(program, mol.split(',')[0][1:])
    else:
        launcher = 'atomfrag_{}.py'.format(program)
        output_prefix = '{0}_{1}'.format(program, mol)

    launcher_full_path = path.normpath(path.join(options.binary_dir, launcher))

    command = []
    command.append('python3 {0}'.format(launcher_full_path))
    command.append('{0} {1}'.format(mol, inp))
    if extra_args is not None:
        command.append(extra_args)

    full_command = ' '.join(command)


    relative_reference_path = 'result'

    return launcher, full_command, output_prefix, relative_reference_path
