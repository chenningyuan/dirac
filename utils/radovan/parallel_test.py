#! /usr/bin/python

# radovan (2009-12-08):
# - i use this script to quickly run many tests, really in parallel
#   (all submitted at the same time)
#   i use it to quickly check whether a fresh compilation on the cluster is ok for production
#   maybe it can be useful for you too (but please change email_recipient)
# - it can be used in the following way (example):
#   rm -f test/run.sh
#   ln -s ../utils/radovan/parallel_test.py test/run.sh
#   cd test
#   for dir in 46* 05*; do
#       cd $dir
#       ./run.sh
#       cd ..
#   done
# - it will submit all selected tests to the queue at the same time
# - for each test that fails i will get an email containing the filtered diff
# - obviously, you need to adapt the script below "script = ..."
# - it has to fail for tests that contain a series of calculations that depend on each other

email_recipient = 'radovan.bast@uit.no'
nr_cores = 8

#-------------------------------------------------------------------------------

import string
import os
import sys

script = '''#!/bin/bash

#PBS -S /bin/bash
#PBS -A nn4654k
#PBS -l nodes=nr_cores:ib
#PBS -l walltime=00:10:00
#PBS -l pvmem=2000mb

export LD_LIBRARY_PATH=/global/apps/intel/mkl/10.2.0.013/lib/em64t:$LD_LIBRARY_PATH

pam=/home/bast/dirac/pam_test

scratch=/global/work/$USER/$PBS_JOBID
mkdir -p $scratch

cd $PBS_O_WORKDIR
pam_command

cd $PBS_O_WORKDIR
../crep result/output result/filter >            result/output.rep
../crep        output result/filter > output.rep
../bdrep                              output.rep result/output.rep > diff_$PBS_JOBID\.txt

if [ -s diff_$PBS_JOBID\.txt ]; then
  cat diff_$PBS_JOBID\.txt | mail -s "test directory output failed" email_recipient
fi

cd $PBS_O_WORKDIR
rm -rf $scratch

exit 0'''

#-------------------------------------------------------------------------------

directory = string.split(os.getcwd(), '/')[-1]

molecule_file = os.environ.get('MOL')
input_file    = os.environ.get('MENU')

molecule_file_no_suffix = os.path.splitext(os.path.basename(molecule_file))[0]
input_file_no_suffix    = os.path.splitext(os.path.basename(input_file))[0]

if molecule_file == input_file:
   output = input_file_no_suffix + '.out'
else:
   output = input_file_no_suffix + '_' + molecule_file_no_suffix + '.out'

pam_command = "$pam -keepwrk -noarch -wrkdir $scratch -mpi %i -mem 1200 -slavemem 1200 -aw 2" % nr_cores

for s in sys.argv[1:]:
    pam_command += ' ' + s

pam_command += ' ' + molecule_file
pam_command += ' ' + input_file

script = string.replace(script, 'pam_command',     pam_command)
script = string.replace(script, 'nr_cores',        str(nr_cores))
script = string.replace(script, 'output',          output)
script = string.replace(script, 'directory',       directory)
script = string.replace(script, 'email_recipient', email_recipient)

run_file = 'run_' + directory
f = open(run_file, 'w')
f.write(script)
f.close()

command = 'chmod 700 ' + run_file + '; qsub ' + run_file + '; rm -f ' + run_file
os.system(command)
