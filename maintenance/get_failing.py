#!/usr/bin/env python

import sys

if sys.version < '2.6':
    print('Requires python version >= 2.6')
    sys.exit(1)
    
try:
    import re, requests, subprocess
except:
    sys.exit('To continue install Python requests module:\n1. option: pip install requests \n2. option: easy_install requests')

"""
This script will parse local file build/Testing/Temporary/LastTestsFailed.log
or dashboard at https://testboard.org/cdash/viewTest.php for defined build id
and create and run a list of failed tests with or without additional arguments for ctest, 
for example: -j4, -V.

Consider you **want** to run this script from build directory:
Usage (Linux):      python ../maintenance/get_failing.py buildid|local "param1 param2 ..."|none
Usage (Windows):    python ..\maintenance\get_failing.py buildid|local "param1 param2 ..."|none

If you want to use only one parameter without spaces for ctest then using of "" is not necessary.

If you are using CMake/CTest version 3.0 or later you can use to run locally failed tests this command:
    ctest --rerun-failed [param1 param2 ...](see CTest docs)
    
Note:
Project xcfun and its tests are not included.
If you use old build id from dashboard some old tests will not be recognized. 

Example (sequential tests of locally failed tests without arguments for ctest):
  $ cd build
  $ python ../maintenance/get_failing.py local none
  
Example (sequential tests of locally failed tests run on 4 cores with verbose output):
  $ cd build
  $ python ../maintenance/get_failing.py local "-j4 -V"
  
Example (sequential tests of tests obtained from dashboard (buildid 1110) without arguments for ctest):
  $ cd build
  $ python ../maintenance/get_failing.py 1110 none
  
Example (sequential tests of tests obtained from dashboard (buildid 1110) run on 4 cores with verbose output):
  $ cd build
  $ python ../maintenance/get_failing.py 1110 "-j4 -V"
"""

def local():
    tests = []
    
    # read recently failed tests
    for line in open('Testing/Temporary/LastTestsFailed.log', 'r').readlines():
        tests.append(line.split(':')[0])
    
    # add ,
    tests = ','.join(tests[0:])
    
    # ctest expects for -I option: first, last, stride, test, test, test, ...
    print 'command: ' + "ctest -I " + "0,0,1," + tests + " " + params
    try:
        subprocess.call("ctest -I " + "0,0,1," + tests + " " + params, shell=True)
    except:
        sys.exit(1)
    
def dashboard():
    # create request
    payload = {'buildid': buildid, 'onlyfailed' : '1'}
    r = requests.get('https://testboard.org/cdash/viewTest.php', params=payload)
    r.encoding = 'utf-8'

    # names of failed tests from dashboard
    names = []
    for line in r.iter_lines():
        match = re.match(r'^.*href.*name=(.*)&amp;', line)
        if match:
            names.append(match.group(1))

    # order of tests from ctest configuration
    tests = []
    for value in names:
        order = 0
        for line in open('CTestTestfile.cmake', 'r').readlines():
            match = re.match(r'(^add_test)|(^ADD_TEST)', line)
            if match:
                order = order + 1
            match = re.match(r'(^add_test\(' + value + ' )|(^ADD_TEST\(' + value + ' )', line)
            if match:
                tests.append(order)
                break
    
    # sorting
    tests.sort()
    
    # int to string
    pos = 0
    for value in tests:
        tests[pos] = str(tests[pos])
        pos = pos + 1
    
    # add ,
    tests = ','.join(tests[0:])

    # ctest expects for -I option: first, last, stride, test, test, test, ...
    print 'command: ' + "ctest -I " + "0,0,1," + tests + " " + params
    try:
        subprocess.call("ctest -I " + "0,0,1," + tests + " " + params, shell=True)
    except:
        sys.exit(1)

# script logic depends on parameters
if (len(sys.argv) == 3):
    buildid = sys.argv[1]
    params = sys.argv[2]
    
    if (buildid == 'local' and params == 'none'):
        params = ''
        local()
    elif (buildid != 'local' and params == 'none'):
        params = ''
        dashboard()
    elif (buildid == 'local' and params != 'none'):
        local()
    elif (buildid != 'local' and params != 'none'):
        dashboard()
    else:
        sys.exit(1)
else:
    print '\tUsage (Linux):     python ../maintenance/get_failing.py buildid|local "param1 param2 ..."|none'
    print '\tUsage (Windows):   python ..\maintenance\get_failing.py buildid|local "param1 param2 ..."|none'
    sys.exit('\tMore info and examples in source file')

