#! /usr/bin/env python
##  -*- Mode: Python; py-indent-offset: 8 -*-

##------------------------------------------------------
## add_copyright.py
##------------------------------------------------------
##
## Written by
## Christoph Jacob <jacob@few.vu.nl>, Aug 2004
##
##------------------------------------------------------

import sys, os, fortran

major, minor, dum, dum , dum = sys.version_info
if not (major >= 2 and minor >= 3) :
        print "you need Python version >= 2.3 to run this script "
        sys.exit(1)

def process_file (filename, dirac_copyright, options) :

        try:
                f = fortran.ffile (filename, True, False)
        except IOError:
                print "\nfile "+filename+" not found"
                return 1
        except fortran.error.FortranException, ex :
                print "Error processing file"+filename+":"
                print str(ex)
                return 1
                        
        try:
                f.add_copyright (dirac_copyright)

                if options.backup :
                        try:
                                f.backup(options.backup_suffix, options.overwrite_backup)
                                print "saved backup copy of original version"
                        except fortran.error.FortranIOError, ex :
                                print str(ex)
                                print "no backup was saved"

                f.save_strings ()

        except fortran.error.FortranException, ex :
                print "Error processing file"+filename+":"
                print str(ex)
                return 1        

        return 0

## ---- MAIN PROGRAM ----

## parse command line options
from optparse import OptionParser

version="%prog, Version 1.0 \nA tool to add Dirac copyright header to source files \n"+\
        "Written by Christoph Jacob <jacob@few.vu.nl>, Aug 2004\n"

optionparser = OptionParser(usage="%prog [options] filenames\nUse option --help for details", version=version)

optionparser.add_option("--no-backup", action="store_false", dest="backup", default=True,
                        help="Do not save backups of the original file")
optionparser.add_option("-s", "--suffix", action="store", dest="backup_suffix", type="string", default="orig",
                        metavar="SUFFIX", help="suffix for backup copies [default: orig]")
optionparser.add_option("-o", "--overwrite-old-backups", action="store_true", dest="overwrite_backup", default=False,
                        help="Overwrite existing backup files")
optionparser.add_option("-c", "--copyright-file", action="store", dest="copyright_file", default="copyright_fortran",
                        help="The name of the file with the copyright notice to be added")

(options, args) = optionparser.parse_args()

if not os.path.exists(options.copyright_file) :
	optionparser.error("copyright file "+options.copyright_file+" not found")

if len(args) == 0 :
        optionparser.error("you have to give the names of the files to reformat")

## read the copyright file
             
f = file (options.copyright_file, 'r')
dirac_copyright = f.readlines () 
f.close

for filename in args :
        if options.backup_suffix[0] == '.' :
                backup_filename = filename + options.backup_suffix
        else :
                backup_filename = filename + '.' + options.backup_suffix

        if True :
##	if not os.path.exists(backup_filename) :

                print "\nProcessing file "+filename

                if not process_file (filename, dirac_copyright, options) == 0 :
                        failed.append(filename)
                        print "this file is left unchanged"
                else :
                        print "Successfully saved reformatted version of file " + filename


