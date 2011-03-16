#! /usr/bin/env python

import sys
import subprocess

print "#define WHO_COMPILED       \"foo\""
print "#define WHERE_COMPILED     \"foo\""
print "#define COMPILED_DIRAC_DIR \"foo\""
print "#define ILP_DATA_MODEL     \"foo\""
print "#define FCOMP              \"foo\""
print "#define CCOMP              \"foo\""
print "#define LIBRARIES          \"foo\""
print "#define DIRAC_BLAS         \"foo\""
print "#define DIRAC_LAPACK       \"foo\""

# get revision and remote

revision = 'unknown'
remote   = 'unknown'

# try whether this is a svn working copy
p = subprocess.Popen('svn info',
                     shell=True,
                     stdin=subprocess.PIPE,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
out = p.communicate()[0]
if 'Revision' in out:
    # this is a svn working copy
    for line in out.split('\n'):
        if 'Revision' in line:
            revision = line.split()[-1]
        if 'URL' in line:
            remote = line.split()[-1]
else:
    # try whether this is a git svn clone
    p = subprocess.Popen('git log | head -n 200 | grep git-svn-id | head -n 1',
                         shell=True,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out = p.communicate()[0]
    if 'git-svn-id' in out:
        # this is a git svn clone
        remote   = out.split()[1].split('@')[0]
        revision = out.split()[1].split('@')[1]

print '#define DIRAC_SVN_URL      "%s"' % revision
print '#define DIRAC_SVN_REVISION "%s"' % remote

sys.exit(0)
