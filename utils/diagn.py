#! /usr/bin/env python
#
# Script aimed to replace this old makefile command by wise python equivalent
#================================================================================
#        @echo -e "\tAnalyzing the dirac.x executable "
#        @echo -e "Size of all static allocations: \c "; \
#          size $(TARGET) 2> /dev/null | grep $(TARGET) | awk '{printf $$3/(1024*1024)}'; \
#          echo ' MB' || echo 'can not be determined'
#        @echo -e "Sorted list of all object sizes (patience please):"
# ... miro: some awk's are not recognizing '--non-decimal-data' !
#           todo: try to make it universal for all architectures
#        @readelf -s $(TARGET) 2> /dev/null | awk --non-decimal-data '{printf "%12.3f %s    %s\n",(($$3)+0)/(1024),"kB",$$8}' | sort -n -r -u || \
#        echo 'can not be determined, readelf not present'
#================================================================================
#
# Script depends on linux commands 'size' and 'readelf'.
#
# TODO: better formatting, count all objects/data sizes
#
#
#   Written by Miro Ilias, Banska Bystrica, 2010
#
#================================================================================


import string, re, sys, os, subprocess

#
# return the integer value of a hexadecimal string s, or return zero if string is wrong
#
def hex2dec(s):
  if  re.match("0x[0-9a-fA-F]+$", s) or re.match("[0-9]+$", s):
   return int(s, 16)
  else :  # no proper hex number - return 0
   return 0

dirac_exe="dirac.x"

if not os.path.isfile(dirac_exe):
   print '\nExecutable dirac.x does not exist for analyze! Exit\n'
   sys.exit(-1)

##### Process size ######
command = "size "+dirac_exe+" 2> /dev/null | grep dirac.x || echo 'can not be determined' "
process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
istat=os.waitpid(process.pid, 0)
output = process.stdout.read().strip()

if istat[1]==0 :
  print '\nTotal size of all static allocations in the dirac.x executable:\n ',\
  float(output.split()[2])/(1024*1024),'MB (what makes ',float(output.split()[2])/(1024), 'kB)'
else :
  print 'can not determine'; sys.exit('exit!')

#### Process readelf ####

print 'Sizes of  individual objects:'

command = " readelf -s "+dirac_exe+" 2> /dev/null || echo 'can not be determined, readelf not present'"
process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

output = process.stdout.read().splitlines(False)

### substitute strings-hex numbers by converted decimal numbers
data_as_lists = [ line.strip().split() for line in output ]  ## convert to a list of lists
for line in data_as_lists:
  if (len(line) >= 8) :
     numstr=string.join(line[2:3]); #print 'numstr=',numstr
     # replace hexadecimal value with decimal value (not string, rather float)
     decstr=hex2dec(numstr);       # print 'decstr=',decstr
     # format float value (into MB) and insert it
     float_val="%11.7f" % (float(decstr)/(1024*1024))
     line.insert(3, float_val )
     del line[2:3] # remove list member

#### sort according to the size of the object - position 2:3  ####
data_as_lists.sort(key=lambda data_as_lists:data_as_lists[2:3], reverse=True)

# selected printing
i=1; imax=100
for el in data_as_lists : 
  if (len(el) >= 8) :
    # avoid printing duplicities, limits the size of the output by imax #
    if (i>1 and el[2:3]==el1[2:3] and el[7:8]==el1[7:8] and i<=imax ):
       #obj_size=(str(el[2:3]));
       print_string=str(el[2:3])+" MB  "+str(el[7:8])+"----"+str(el[3:4])
       print "* %s "%print_string
    el1=el
    i=i+1 # increase counter
