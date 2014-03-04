
import os
import sys
from runtest import TestRun

class DiracTestRun(TestRun):

    def __init__(self, _file, argv):
        TestRun.__init__(self, _file, argv)

    def execute(self, inp_files, mol_files, f=None, args='', accepted_errors=[]):
        pam_script = os.path.normpath(os.path.join(self.binary_dir, 'pam'))
        if not os.path.exists(pam_script):
            sys.stderr.write('ERROR: pam script not found in %s\n' % pam_script)
            sys.stderr.write('       run %s --help\n' % sys.argv[0])
            sys.exit(-1)
        launcher = 'python "%s"' % pam_script
        launcher += ' --noarch --nobackup'
        if args != '':
            launcher += ' %s' % args
        for inp in inp_files:
            inp_no_suffix = os.path.splitext(inp)[0]
            for mol in mol_files:
                mol_no_suffix = os.path.splitext(mol)[0]
                out = '%s_%s.out' % (inp_no_suffix, mol_no_suffix)
                sys.stdout.write('\nrunning test: %s %s\n' % (inp_no_suffix, mol_no_suffix))
                TestRun.execute(self,
                                command = launcher + ' --inp=%s --mol=%s' % (inp, mol),
                                accepted_errors = accepted_errors)
                if f != None:
                    f.check(self.work_dir, '%s' % out, 'result/%s' % out)
                    sys.stdout.write('passed\n')
                else:
                    sys.stdout.write('finished (no reference)\n')
