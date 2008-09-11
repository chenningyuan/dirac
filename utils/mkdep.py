#! /usr/bin/env python
import re, sys, os, sets

cppcmd = ' '.join(sys.argv[1:])

use_re     = re.compile('\s*use\s+(\w+)',                    re.I)
module_re  = re.compile('\s*module(?!\s+procedure)\s+(\w+)', re.I)
include_re = re.compile('^#( |line )1 \"([^<"]*)\"',             re.I)

def scan_file(path,cppcmd):
    cpp = os.popen(cppcmd+' '+path)
    use_mod = sets.Set()
    gen_mod = []
    includes = sets.Set()
    for line in cpp.readlines():
        m = use_re.match(line)
        if m:
            use_mod.add(m.group(1).lower())
        else:
            m = module_re.match(line)
            if m and m.group(1):
                gen_mod.append(m.group(1).lower())
            else:
                m = include_re.match(line)
                if m:
                    includes.add(m.group(2))
    return (includes,use_mod,gen_mod)

def objpath(srcpath):
    return srcpath[:srcpath.rfind(".")]+".o"

filedat = {}

for f in sys.stdin.readlines():
    filedat[f] = scan_file(f,cppcmd)

module_def = {}

for key,val in filedat.items():
    for m in val[2]:
        module_def[m] = key

for key,val in filedat.items():
    print objpath(key)+":",
    for incf in val[0]:
        if len(incf) > 0 and incf[0] != '/': # Do not list dependencies on system headers 
            print incf,
    for m in val[1]:
        if module_def.has_key(m):
            if module_def[m] != key:
                print objpath(module_def[m]), # do not print self interaction
        else:
            sys.stderr.write("Warning: cannot find definition of Fortran module %s.\n" % m)
    print
