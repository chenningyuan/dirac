#!/usr/bin/env python

# radovan bast <radovan.bast@irsamc.ups-tlse.fr>
# script to convert fixed form fortran files (written in Fortran77) and includes to free form

import sys
import string

if len(sys.argv) == 1:
    print 'usage: %s fixed-form-file > free-form-file' % sys.argv[0]
    sys.exit()

s = open(sys.argv[1], 'r').read().splitlines()

for i in range(len(s)):
    if len(s[i]) > 0:
        if s[i][0] == 'C':
            # replace C-comment by !
            s[i] = '!' + s[i][1:]
        if s[i][0] == 'c':
            # replace c-comment by !
            s[i] = '!' + s[i][1:]
        if s[i][0] == '*':
            # replace *-comment by !
            s[i] = '!' + s[i][1:]
    if i+1 < len(s):
        if len(s[i+1]) > 5:
            # check the presence of a continuation marker which can be any ASCII character (not only '&')
            if s[i+1][0:4]=='    ' and  s[i+1][5] != ' ':
                # add continuation marker to previous line and verify the line length
                if not '&' in s[i][6:] and len(s[i])>0: # unless it is already there
                    if not s[i][0] == '#': # unless this line starts with CPP
                        # replace the continuation marker with blank character
                        s[i+1]=' ' + s[i+1][6:]
                        for j in range(72 - len(s[i])):
                            s[i] += ' '
                        s[i] += '&'
    print s[i]
