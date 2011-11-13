#!/usr/bin/env python

import os
import sys

def create_build_dir(build_tag,build_dir,source_dir):
   """Creates a directory whose files point to the cvs source code.

      This function will perform all the steps required to create
      a given subdirectory from the cvs source code in dirac.
      This amounts to finding out which files are in the cvs for
      each module directory, and making links to those files in the
      build directory structure created
   """

   verbose = True

   print " ... looking up source directory: ",source_dir
   os.chdir(source_dir)

   print "     ... processing source files"
   list = []
   list_files(".",list)

   if verbose:
      print "     ... files found:"
      for entry in list:
         print "         ... file: ",entry

   print "\n ... checking if directory",build_dir,"exists",
   if os.access(build_dir,0) == False:
      print " ... it does not, making directory",build_dir
      os.mkdir(build_dir)
   else:
      print " ... it does"

   os.chdir(build_dir)
#  os.mkdir("basis")
#  os.mkdir("basis/MOLFDIR")

   for entry in list:
      src_file = os.path.join(source_dir,entry)
      dst_file = entry
      (dst_dir,filename) = os.path.split(entry)
      full_dst_dir = build_dir+os.sep+dst_dir
      print "\n ... checking if directory",full_dst_dir,"exists",
      if os.access(full_dst_dir,0) == False:
         print " ... it does not, making directory"
         os.makedirs(full_dst_dir)
      else:
         print " ... it does"
      
      print " ... symlinking",src_file," -> ",dst_file,
      if os.access(dst_file,0) == False:
         os.symlink(src_file,dst_file)
         print " ... done"
      else:
         print " ... skipped, linked file already present"

def list_files(dir,list):
   """Goes over a directory tree and lists the files.

      At this moment skipping CVS directories is hardcoded.
   """

   if os.path.isdir(dir):
      direntries = os.listdir(dir)
      for entry in direntries:
         if entry != "CVS":
            if not os.path.isfile(entry):
               newentry = dir+os.sep+entry
               print "new entry", newentry
               list_files(newentry,list)

            if dir == ".":
               name = entry
            else:
               name = dir[2:]+os.sep+entry
            list.append(name) 


def main(argv=None):
 
# default values
   source_dir     = os.path.expandvars("$HOME/devel/cvs/snapshots/latest/Dirac")
   build_basedir  = os.path.expandvars("$HOME/devel/DIRAC_BUILDS")  
   build_name     = "darwin_intel_g95"

   for arg in sys.argv:
      parsearg=arg.split("=")
      if len(parsearg) == 2:
         if parsearg[0] == "tag":
            build_name    = parsearg[1]
         elif parsearg[0] == "bdir":
            build_basedir = parsearg[1]
         elif parsearg[0] == "src":
            source_dir    = parsearg[1]
         else:
            print "\n unknown key=value combination. exiting"
            sys.exit(-1)

   build_dir = build_basedir+os.sep+build_name 

   create_build_dir(build_name,build_dir,source_dir)
   print "\n ... done"


if __name__ == "__main__":
    main()
