#!/usr/bin/env python
import sys, re, os
    
"""A C pre-preprocessor to take out certain sections of source
files, depending on the preprocessor directives. This can be used in the
future to prepare a release. 

For example, to remove the code which would be included only if -DMOD_DFT
is given you could do

cppp.py -UMOD_DFT */*.F

This will completely remove all sections

#ifdef MOD_DFT
...
#endif

It is also possible to use -DMOD_DEF, which will remove the #ifdef but
not the content between the #ifdef/#endif.

Currently the script can handle 

#ifdef XYZ
#if defined (XYZ)
#ifndef XYZ
#if !defined (XYZ)
#endif

The script cannot handle compound statements (#if !defined(A) && !defined(B)),
so it will stop if such a line contains some of the macros we are removing.
Macros with an actual value (-DMOD_DFT=213 would not work).

WARNING: The script will modify the files in-place, so run it on a copy of 
your source tree!
"""


def usage():
    print 'Usage: cppp [-DSYM ..] [-USYM ..] FILES'
    print ' A C pre-preprocessor'
    print ' Modifies FILES by assuming that some macros are always'
    print ' defined (-D), and some macros never defined (-U).'
    print ' The corresponding #if/#endif and enclosed text is removed.'
    print ' The files are modified in-place!'
    print ' This script requires Python 2.3 or newer.'
    print ' By <ulfek@few.vu.nl> April 2008'

def tokenize(line):
    re_split = re.compile('\W')
    return re_split.split(line)

def union(list1,list2):
    l = []
    for i1 in list1:
        if i1 in list2:
            l.append(i1)
    return l

def do_filter(defined, undefined, filename, infile, outfile):
    re_ifdef = re.compile('^#\W*if(def| defined\W*[\(]*)\W*(\w*)')
    re_ifndef = re.compile('^#\W*if(ndef|\W*!\W*defined\W*[\(]*)\W*(\w*)')
    re_endif = re.compile('^#[ \t]*endif.*')
    re_complicated = re.compile('^#\W*if[^/]*[&|].*$')
    line_nr = 1
    cond_stack = [] # Append 'D' for defined and 'U' for undefined, 'I' for other #if things.    
    for line in infile.xreadlines():
        if re_complicated.match(line):
            tokens = tokenize(line)
            if len(union(defined,tokens)) > 0 or len(union(undefined,tokens)) > 0:
                print 'Error: %s:%i' % (filename,line_nr)
                print '  Line too complicated and it contains my tokens:',line.rstrip()
                print '  Fix the cppp script or make a simpler #if.. statement'
                return -1
            else:
                if 'U' not in cond_stack:
                    outfile.write(line)
        else:
            m = re_ifdef.match(line)
            if m:
#                print '---If line:',line.rstrip()
#                print '---macro: "%s"' % m.group(2)
                if m.group(2) in defined:
                    cond_stack.append('D')
                elif m.group(2) in undefined:
                    cond_stack.append('U')
                else:
                    cond_stack.append('I')
                    outfile.write(line)
            else:
                m = re_ifndef.match(line)
                if m:
#                    print '---Ifndef line:',line.rstrip()
#                    print '---macro: "%s"' % m.group(2)
                    if m.group(2) in defined:
                        cond_stack.append('U')
                    elif m.group(2) in undefined:
                        cond_stack.append('D')
                    else:
                        cond_stack.append('I')
                        outfile.write(line)
                else:
                    m = re_endif.match(line)
                    if m:
                        t = cond_stack.pop() # Might want to check for unbalanced #endif here
                        if t == 'I': # write #ifdef only if we also wrote the #if.. 
                            outfile.write(line)
                    else:
                        if 'U' not in cond_stack:
                            outfile.write(line)
        line_nr += 1
    return 0
            

def main():
    defined = []
    undefined = []
    filenames = []

    target = filenames
    for arg in sys.argv[1:]:
        if arg.find('-D') == 0:
            arg = arg[2:]
            target = defined
        elif arg.find('-U') == 0:
            arg = arg[2:]
            target = undefined
        if len(arg) > 0:
            target.append(arg)
            if target in [defined, undefined]:
                if '=' in arg:
                    print 'Error: Giving a symbol a value (%s) is not supported. Quitting!' % arg
                    return -1
            target = filenames

#    print 'Defined:',defined
#    print 'Undefined:',undefined
#    print 'Files',filenames
    for fn in filenames:
        ofname = fn+'.cppptmp'
        f = open(fn,'r')
        of = open(ofname,'w')
        print fn
        res = do_filter(defined,undefined,fn,f,of)
        if res == 0:
            os.rename(ofname,fn)
        else:
            os.remove(ofname)
            

if len(sys.argv) < 2 or sys.argv[1] in ['-h','--help','-help']:
    usage()
else:
    main()
    
