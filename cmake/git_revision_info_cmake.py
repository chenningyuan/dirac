#!/usr/bin/env python
import sys
import os
import string
import subprocess
def remove_email(str):
    ''' remove email address from entering string '''
    start=str.find('<')
    end=str.find('>')
    email=str[start:end+1]
    str=str.replace(email,'')
    return str
command='git show' # instead of 'git rev-parse HEAD 2> /dev/null'
p=subprocess.Popen(command,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
stdout, stderr = p.communicate()
if stderr:
    print 'unknown'
else: # process first three lines of the standard output
    all_out=stdout.splitlines()[:3]
    git_rev=all_out[0].split()
    commiter_id=remove_email(all_out[1])
    who_when=commiter_id+' '+all_out[2]
    who_when=who_when.replace('Author:','')
    who_when=who_when.replace('Date:',' on')
    who_when=who_when.replace('  ','')
    git_rev_output=git_rev[1]+' (last commit by'+who_when+')'
    print git_rev_output
