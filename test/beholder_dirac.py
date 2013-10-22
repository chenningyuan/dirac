
import os
import sys
from beholder import TestRun

class DiracTestRun(TestRun):

    def __init__(self, _file, argv):
        TestRun.__init__(self, _file, argv)

    def execute(self, inp_files, mol_files, f=None, args='', accepted_errors=[]):
        launcher = os.path.join(self.binary_dir, 'pam')
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
                                success_message = 'exit           : normal',
                                accepted_errors = accepted_errors)
                if f != None:
                    f.check(self.work_dir, '%s' % out, 'result/%s' % out)
                sys.stdout.write('passed\n')
