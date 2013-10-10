
"""
    Beholder - lightweight regression testing

    Copyright (c) 2013, Radovan Bast (lastname at kth.se)
    Licensed under the GNU Lesser General Public License.

    For documentation see https://github.com/rbast/beholder
"""

import os

from beholder import TestRun

class DiracTestRun(TestRun):

    def __init__(self, _file, argv):
        TestRun.__init__(self, _file, argv)

    def execute(self, inp_files, mol_files, f=None):
        launcher = os.path.join(self.binary_dir, 'pam')
        launcher += ' --noarch --nobackup'
        for inp in inp_files:
            for mol in mol_files:
                output = '%s_%s.out' % (inp, mol)
                TestRun.execute(self,
                                command = launcher + ' --inp=%s.inp --mol=%s.mol' % (inp, mol),
                                success_message = 'exit           : normal')
                if f != None:
                    f.check(self.work_dir, '%s' % output, 'result/%s' % output)
