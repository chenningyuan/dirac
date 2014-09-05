#!/usr/bin/env python

"""
This script will parse Testing/Temporary/LastTestsFailed.log
and create a list of failed tests to run to be used
with ctest -I.

Example (sequential tests run on 4 cores):
  $ cd build
  $ ctest -j4 -I `python ../maintenance/get_failing.py`
"""

tests = []
for line in open('Testing/Temporary/LastTestsFailed.log', 'r').readlines():
    tests.append(line.split(':')[0])

first = tests[0]
rest = ','.join(tests[1:])

# ctest expects fror -I option: first, last, stride, test, test, test, ...
print('%s,%s,,%s' % (first, first, rest))
