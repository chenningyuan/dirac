#! /usr/bin/env python
import re, sys, os

class MySet:
    """A subset of sets.Set from python 2.4"""
    def __init__(self):
        self.things = {}
    def add(self,x):
        self.things[x] = 1;
    def __iter__(self):
        """This is to allow if x in set"""
        return self.things.__iter__()

cppcmd = ' '.join(sys.argv[1:])

use_re     = re.compile('\s*use\s+(\w+)',                    re.I)
module_re  = re.compile('\s*module(?!\s+procedure)\s+(\w+)', re.I)
include_re = re.compile('^#( |line )1 \"([^<"]*)\"',             re.I)

def scan_file(path,cppcmd):
    cpp = os.popen(cppcmd+' '+path)
    use_mod = MySet()
    gen_mod = []
    includes = MySet()
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

def objpath(srcpath,dot_suffix):
    return srcpath[:srcpath.rfind(".")]+dot_suffix

filedat = {}

for f in sys.stdin.readlines():
    filedat[f] = scan_file(f,cppcmd)

module_def = {}

for key,val in filedat.items():
    for m in val[2]:
        module_def[m] = key

# put here all suffixes used in the DIRAC make system
dot_suffixes = [".o",".s",".r",".p",".opt"]

for key,val in filedat.items():
  for this_dot_suffix in dot_suffixes:
    print objpath(key,str(this_dot_suffix))+":",
    for incf in val[0]:
        if len(incf) > 0 and incf[0] != '/': # Do not list dependencies on system headers
            print incf,
    for m in val[1]:
        if module_def.has_key(m):
            if module_def[m] != key:
                print objpath(module_def[m],str(this_dot_suffix)), # do not print self interaction
        else:
            sys.stderr.write("Warning: cannot find definition of Fortran module %s.\n" % m)
    print
