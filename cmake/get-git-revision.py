#!/usr/bin/env python

import re
import string
import subprocess

p = subprocess.Popen('git show', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()

if stderr:
    print 'unknown'
else: # process first three lines of the standard output
    commit = string.replace(stdout.splitlines()[0], 'commit ', '')
    author = string.replace(re.sub(r'<.*>', '', stdout.splitlines()[1]), 'Author: ', '')
    date   = string.replace(stdout.splitlines()[2], 'Date:   ', '')
    print '%s (%s- %s)' % (commit, author, date)
