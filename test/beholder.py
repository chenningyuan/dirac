
"""
    Beholder - lightweight regression testing

    Copyright (c) 2013, Radovan Bast (lastname at kth.se)
    Licensed under the GNU Lesser General Public License.

    For documentation see https://github.com/rbast/beholder
"""

import re
import os
import sys
import argparse
import subprocess
import shlex
import shutil

#-------------------------------------------------------------------------------
class TestRun():

    #---------------------------------------------------------------------------
    def __init__(self, _file, argv):
        self.input_dir = input_dir = os.path.dirname(os.path.realpath(_file))
        self.binary_dir, self.work_dir = self._parse_args(input_dir, argv)
        if self.work_dir != self.input_dir:
            self._safe_copy(self.input_dir, self.work_dir)
        os.chdir(self.work_dir)

    #---------------------------------------------------------------------------
    def execute(self, command, success_message=''):
        try:
            process = subprocess.Popen(shlex.split(command),
                                       stdin=subprocess.PIPE,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if stderr != '':
                sys.stderr.write('ERROR: crash during %s\n' % command)
                sys.exit(-1)
            if success_message != '':
                if success_message not in stdout:
                    sys.stderr.write('ERROR: crash during %s\n' % command)
                    sys.exit(-1)
        except OSError:
            sys.stderr.write('ERROR: could not execute command %s\n' % command)
            sys.exit(-1)

    #---------------------------------------------------------------------------
    def _safe_copy(self, root_src_dir, root_dst_dir, exclude_files=[]):
        for src_dir, dirs, files in os.walk(root_src_dir):
            dst_dir = src_dir.replace(root_src_dir, root_dst_dir)
            if not os.path.exists(dst_dir):
                os.makedirs(dst_dir)
            for f in files:
                if f not in exclude_files:
                    src_file = os.path.join(src_dir, f)
                    dst_file = os.path.join(dst_dir, f)
                    shutil.copy(src_file, dst_file)

    #---------------------------------------------------------------------------
    def _parse_args(self, input_dir, argv):
        parser = argparse.ArgumentParser(description='Beholder - lightweight regression testing.')
        parser.add_argument('--binary-dir',
                           action='store',
                           default=input_dir,
                           help='directory containing the binary/launcher [default: %(default)s]')
        parser.add_argument('--work-dir',
                           action='store',
                           default=input_dir,
                           help='working directory [default: %(default)s]')
        args = parser.parse_args(args=argv[1:])
        return (args.binary_dir, args.work_dir)

#-------------------------------------------------------------------------------
class _SingleFilter():

    def __init__(self, *args, **kwargs):
        self.from_string    = kwargs.get('from_string')
        self.to_string      = kwargs.get('to_string')
        self.tolerance      = kwargs.get('tolerance', 1.0e-5)
        self.ignore_sign    = kwargs.get('ignore_sign', False)
        self.ignore_below   = kwargs.get('ignore_below', 1.0e-30)
        self.start_is_regex = kwargs.get('start_is_regex', False)
        self.end_is_regex   = kwargs.get('end_is_regex', False)

#-------------------------------------------------------------------------------
class Filter():

    #---------------------------------------------------------------------------
    def __init__(self):
        self.filter_list = []

    #---------------------------------------------------------------------------
    def add(self, *args, **kwargs):
        self.filter_list.append(_SingleFilter(*args, **kwargs))

    #---------------------------------------------------------------------------
    def check(self, work_dir, output_file, reference_file):

        os.chdir(work_dir)

        log_diff      = open('%s.diff'      % output_file, 'w')
        log_result    = open('%s.result'    % output_file, 'w')
        log_reference = open('%s.reference' % output_file, 'w')

        output     = open(output_file).readlines()
        output_ref = open(reference_file).readlines()
        for f in self.filter_list:
            line_l     = self._find_from_to(output, output_file, f)
            line_l_ref = self._find_from_to(output_ref, reference_file, f)
            for x in line_l:
                for line in x:
                    log_result.write(output[line])
            for x in line_l_ref:
                for line in x:
                    log_reference.write(output_ref[line])

            """
            in:  output      a list of strings (= a file)
                 output_ref  a list of strings (= a file)
                 line_l      list of integers that correspond to line numbers in python numbering
                 line_l_ref  list of integers that correspond to line numbers in python numbering
            if nr_of_floats > 0 then it will only compare the first nr_of_floats in each line
            """
            nr_of_floats = 0 # FIXME hardcoding

            numeric_const_pattern = r"""
            [-+]? # optional sign
            (?:
                (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
                |
                (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
            )
            # followed by optional exponent part if desired
            (?: [EeDd] [+-]? \d+ ) ?
            """

            p  = re.compile(numeric_const_pattern, re.VERBOSE)
            pd = re.compile(r'[dD]')

            f_l           = []
            f_l_ref       = []
            f_to_line     = []
            f_to_line_ref = []
            for x in line_l:
                for line in x:
                    i = 0
                    for w in output[line].split():
                        i += 1
                        if (nr_of_floats > 0) and (i > nr_of_floats): break
                        # do not consider words like TzB1g
                        # otherwise we would extract 1 later
                        if re.match(r'^[0-9\.eEdD\+\-]*$', w):
                            # apply floating point regex
                            for m in p.findall(w):
                                 # substitute dD by e
                                 m = pd.sub('e', m)
                                 f_l.append(float(m))
                                 f_to_line.append(line)

            for x in line_l_ref:
                for line in x:
                    i = 0
                    for w in output_ref[line].split():
                        i += 1
                        if (nr_of_floats > 0) and (i > nr_of_floats): break
                        if re.match(r'^[0-9\.eEdD\+\-]*$', w):
                            for m in p.findall(w):
                                 m = pd.sub('e', m)
                                 f_l_ref.append(float(m))
                                 f_to_line_ref.append(line)

            if len(f_l) == len(f_l_ref):
                for i in range(len(f_l)):
                    r     = f_l[i]
                    r_ref = f_l_ref[i]

                    if f.ignore_sign:
                        # if ignore sign take absolute values
                        r     = abs(r)
                        r_ref = abs(r_ref)
                    error = r - r_ref

                    rel_error = error
                    if abs(r_ref) > f.ignore_below:
                        # calculate relative error only for significant ('nonzero') numbers
                        rel_error /= r_ref
                    if abs(rel_error) > f.tolerance:
                        log_diff.write('line %i: %s' % (f_to_line[i]+1, output[f_to_line[i]]))
                        log_diff.write('    rel error %7.4e > %7.4e\n\n' % (rel_error, f.tolerance))
            else:
                log_diff.write('compare floats sizes do not match\n')

                log_diff.write('own gave %i floats:\n' % len(f_l))
                last_line_printed = -1
                for i in range(len(f_l)):
                    if (f_to_line[i] != last_line_printed):
                       log_diff.write('    %s' % output[f_to_line[i]])
                       last_line_printed = f_to_line[i]

                log_diff.write('refence gave %i floats:\n' % len(f_l_ref))
                last_line_printed = -1
                for i in range(len(f_l_ref)):
                    if (f_to_line_ref[i] != last_line_printed):
                       log_diff.write('    %s' % output_ref[f_to_line_ref[i]])
                       last_line_printed = f_to_line_ref[i]
        log_diff.close()
        log_result.close()
        log_reference.close()
        if os.path.getsize('%s.diff' % output_file) != 0:
            sys.stderr.write('ERROR: test %s failed\n' % output_file)
            log_diff = open('%s.diff' % output_file, 'r')
            for line in log_diff.readlines():
                sys.stderr.write(line)
            log_diff.close()
            sys.exit(-1)

    #---------------------------------------------------------------------------
    def _find_from_to(self, output, output_file, f):

        if f.start_is_regex:
            ps = re.compile(r'.*%s' % f.from_string)
        if f.end_is_regex:
            pe = re.compile(r'.*%s' % f.to_string)

        line_l = []
        for i in range(len(output)):
            start_line_matches = False
            if f.start_is_regex:
                start_line_matches = ps.match(output[i])
            else:
                start_line_matches = (f.from_string in output[i])
            if start_line_matches:
                for j in range(i, len(output)):
                    f.end_line_matches = False
                    if f.end_is_regex:
                        f.end_line_matches = pe.match(output[j])
                    else:
                        f.end_line_matches = (f.to_string in output[j])
                    if f.end_line_matches:
                        line_l.append(range(i, j))
                        break

        if line_l == []:
            sys.stderr.write('ERROR: filter did not extract anything from file %s\n' % output_file)
            sys.exit(-1)

        return line_l

    #---------------------------------------------------------------------------
    def _find_string_plus_nr_lines(self, output, start_string, n, start_is_regex=False):

        if start_is_regex:
            p = re.compile(r'.*%s' % start_string)

        line_l = []
        for i in range(len(output)):
            start_line_matches = False
            if start_is_regex:
                start_line_matches = p.match(output[i])
            else:
                start_line_matches = (start_string in output[i])
            if start_line_matches:
                line_l.append(range(i, i+n))

        return line_l
