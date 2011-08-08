#! /usr/bin/env python
#
# Get (as stdout print) the compiler version of Fortran, C, C++ compiler
#

import sys
import os

class compiler_version:
    def __init__(self):
         print 'initialzie...'
         self.Fortran_compiler=None
         self.Fortran_compiler_version=None
         self.C_compiler=None
         self.C_compiler_version=None
         self.CPP_compiler=None
         self.CPP_compiler_version=None
         self.get_version=None

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
        print 'fc=',self.Fortran_compiler
        if self.Fortran_compiler=="GNU":
            self.get_version="gfortran --version"
        elif self.Fortran_compiler=="G95":
            self.get_version="g95 --version"

    def get_C_compiler_version_command(self):
        """ GNU, G95, Intel, PGI, XL   """
        print 'cc=',self.C_compiler

    def get_CPP_compiler_version_command(self):
        """ GNU, G95, Intel, PGI, XL   """
        print 'cpp=',self.CPP_compiler

    def process_compiler_version_command(self):
        print 'command=',self.get_version

    def get_compiler_version(self):
         """ Get compiler version """
         if self.Fortran_compiler:
             self.get_Fortran_compiler_version_command()
             self.process_compiler_version_command()
         if self.C_compiler:
             self.get_C_compiler_version_command()
             self.process_compiler_version_command()
         if self.CPP_compiler:
             self.get_CPP_compiler_version_command()
             self.process_compiler_version_command()


# "main" part
cv=compiler_version()
cv.read_flags()
cv.get_compiler_version()

