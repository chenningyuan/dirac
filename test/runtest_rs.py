"""
 Functions used for testing real-space data;
 code based on runtest_v1.py
 Gosia Olejniczak, 2024
"""

import re
import os
import sys
import subprocess
import shlex
import shutil
import string
from optparse import OptionParser
import runtest_v1 as runtest

try:
    import h5py
    h5py_available = True
except:
    h5py_available = False

if not h5py_available:
    raise runtest.SubprocessError('ERROR: h5py not available.')


class CompareHDF5:

    def __init__(self, work_dir):

        self.work_dir = work_dir
        self.test_name = os.path.dirname(work_dir)

    def verify(self, **kwargs):
        self.filter_list = []
        recognized_kw = ['typ',
                         'gf_label',
                         'gf_name',
                         'gf_test_key',
                         'ftest',
                         'fref',
                         'ignore_sign',
                         'rel_tolerance',
                         'abs_tolerance']

        unrecoginzed_kw = [kw for kw in kwargs.keys() if kw not in recognized_kw]

        if unrecoginzed_kw != []:
            raise runtest.FilterKeywordError('ERROR: keyword(s) (%s) not recognized\n       available keywords: (%s)\n' % ((', ').join(unrecoginzed_kw),(', ').join(recognized_kw)))

        incompatible_pairs = [('rel_tolerance', 'abs_tolerance')]
        incompatible_kw = [(kw1, kw2) for (kw1, kw2) in incompatible_pairs if kw1 in kwargs.keys() and kw2 in kwargs.keys()]
        if incompatible_kw != []:
            raise runtest.FilterKeywordError('ERROR: incompatible keyword pairs: %s\n' % incompatible_kw)

        # now continue with keywords
        self.typ = kwargs.get('typ', '')
        self.gf_label = kwargs.get('gf_label', '')
        self.gf_name = kwargs.get('gf_name', '')
        self.gf_test_key = kwargs.get('gf_test_key', '')
        self.ftest = kwargs.get('ftest', '')
        self.fref = kwargs.get('fref', '')
        self.ignore_sign = kwargs.get('ignore_sign', False)

        if 'rel_tolerance' in kwargs.keys():
            self.rel_tolerance = kwargs.get('rel_tolerance')
            self.rel_tolerance_is_set = True
        else:
            self.rel_tolerance_is_set = False

        if 'abs_tolerance' in kwargs.keys():
            self.abs_tolerance = kwargs.get('abs_tolerance')
            self.abs_tolerance_is_set = True
        else:
            self.abs_tolerance_is_set = False

        self.compare_with_ref()

    def test_hdf5_single(self, fname):
        data_file = h5py.File(fname, 'r')
        try:
            data = data_file[self.gf_label]
        except:
            message = "Label {} not found in {}\n".format(self.gf_label, fname)
            sys.stdout.write(message)
            raise runtest.TestFailedError(message)
        return data
 
    def test_hdf5_array(self, fname):
        data_file = h5py.File(fname, 'r')
        try:
            data = data_file[self.gf_label]
        except:
            message = "Label {} not found in {}\n".format(self.gf_label, fname)
            sys.stdout.write(message)
            raise runtest.TestFailedError(message)
        return data[self.gf_test_key][()]

    def close_values(self, this_value, ref_value, rtol=None, atol=None):
      if self.abs_tolerance_is_set:
          atol = self.abs_tolerance
      else:
          if atol is None:
              #atol = 0.05*ref_value 
              atol = 1.0e-12

      if self.rel_tolerance_is_set:
          rtol = self.rel_tolerance
      else:
          if rtol is None:
              rtol = 1.0e-6

      if self.ignore_sign:
          this_value = abs(this_value)
          ref_value = abs(ref_value)

      is_close = abs(this_value-ref_value) <= max(rtol * max(abs(this_value), abs(ref_value)), atol)
      return is_close
 
    def compare_with_ref(self):
        
        # gosia 19/5/2023
        # alternatively, there is h5diff tool,
        # but this may not always work, e.g.:
        # https://forum.hdfgroup.org/t/identical-hdf5-files-according-to-h5diff-differ-in-comparison-with-diff-unix-command/3285/5
        # so instead, we compare selected datasets

        if self.typ == 'single':
            test = self.test_hdf5_single(self.ftest)
            ref  = self.test_hdf5_single(self.fref)
            res = tuple(i == j for i, j in zip(test, ref))
            if False in res:
                message = "ERROR: test {} failed: Different single data on hdf5 file: {}\n".format(self.test_name, self.ftest)
                sys.stdout.write(message)
                raise runtest.TestFailedError(message)

        elif self.typ == 'array':
            test = self.test_hdf5_array(self.ftest)
            ref  = self.test_hdf5_array(self.fref)
            res = tuple(self.close_values(i, j) for i, j in zip(test, ref))
            if False in res:
                message = "ERROR: test {} failed: Different array data on hdf5 file: {}\n".format(self.test_name, self.ftest)
                sys.stdout.write(message)
                raise runtest.TestFailedError(message)

        else:
            message = "ERROR: unknown keyword"
            raise runtest.BadFilterError(message)

