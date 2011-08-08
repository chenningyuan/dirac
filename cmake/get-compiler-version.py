#! /usr/bin/env python
#
# Get (as stdout print) the compiler version of Fortran, C, C++ compiler
#

import sys
import os

class compiler_version:
    def __init__(self):
         self.Fortran_compiler=None
         self.Fortran_compiler_version=None
         self.C_compiler=None
         self.C_compiler_version=None
         self.CPP_compiler=None
         self.CPP_compiler_version=None
         self.get_version=None
         self.version_output=None
         self.stdout=[]
         self.stderr=None

    def read_flags(self):
         """ read the flags """
         from optparse import OptionParser, OptionGroup
         usage = '''
           Examples:
            %prog --fc=... --cc=... --cpp=..."
         '''
         parser = OptionParser(usage)
         group = OptionGroup(parser, 'Mandatory options')
         group.add_option('--fc',
                 type='string',
                 action='store',
                 dest='fc',
                 default=None,
                 help='set the Fortran compiler [default: None]',
                 metavar='STRING')
         group.add_option('--cc',
                 type='string',
                 action='store',
                 dest='cc',
                 default=None,
                 help='set the C compiler [default: None]',
                 metavar='STRING')
         group.add_option('--cpp',
                 type='string',
                 action='store',
                 dest='cpp',
                 default=None,
                 help='set the C++ compiler [default: None]',
                 metavar='STRING')
         parser.add_option_group(group)
         (options, args) = parser.parse_args()
 
         if options.fc: 
             self.Fortran_compiler=options.fc
         if options.cc: 
             self.C_compiler=options.cc
         if options.cpp: 
             self.CPP_compiler=options.cpp

    def get_Fortran_compiler_version_command(self):
        """ GNU, G95, Intel, PGI, XL   """
        if self.Fortran_compiler=="GNU":
            self.get_version="gfortran --version"
        elif self.Fortran_compiler=="G95":
            self.get_version="g95 --version"
        elif self.Fortran_compiler=="Intel":
            self.get_version="ifort --version"
        elif self.Fortran_compiler=="PGI":
            self.get_version="pgf90 -V"
        elif self.Fortran_compiler=="XL":
            self.get_version="xlf -V"
        else:
            self.get_version="unknown"

    def get_C_compiler_version_command(self):
        """ GNU, G95, Intel, PGI, XL   """
        if self.C_compiler=="GNU":
            self.get_version="gcc --version"
        elif self.C_compiler=="Intel":
            self.get_version="icc --version"
        elif self.C_compiler=="PGI":
            self.get_version="pgcc -V"
        elif self.C_compiler=="XL":
            self.get_version="xlc -V"
        else:
            self.get_version="unknown"

    def get_CPP_compiler_version_command(self):
        """ GNU, G95, Intel, PGI, XL   """
        if self.CPP_compiler=="GNU":
            self.get_version="g++ --version"
        elif self.CPP_compiler=="Intel":
            self.get_version="iCC --version"
        elif self.Fortran_compiler=="PGI":
            self.get_version="pgCC -V"
        elif self.Fortran_compiler=="XL":
            self.get_version="xlCC -V"
        else:
            self.get_version="unknown"

    def process_compiler_version_command(self):
        import subprocess
        p=subprocess.Popen(self.get_version, shell=True, \
          stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (self.version_output,self.stderr)=p.communicate()
        p.wait()

        # use only first line, not the copyright stuff
        self.stdout = self.version_output.split('\n')[0]
        print >> sys.stdout,self.stdout
    
        if self.stderr: # flush out error stream
            print >> sys.stderr,self.stderr

    def get_compiler_version(self):
         """ Get compiler version of Fortran/C/C++ compilers """
         if self.Fortran_compiler:
             self.get_version=None
             self.get_Fortran_compiler_version_command()
             self.process_compiler_version_command()
         if self.C_compiler:
             self.get_version=None
             self.get_C_compiler_version_command()
             self.process_compiler_version_command()
         if self.CPP_compiler:
             self.get_version=None
             self.get_CPP_compiler_version_command()
             self.process_compiler_version_command()


# "main" part
cv=compiler_version()
cv.read_flags()
cv.get_compiler_version()

