#! /usr/bin/env python
#
# Script aimed to replace this makefile command by python equivalent
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


import re, sys, os, subprocess

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
  print '\nTotal size of all static allocations in the dirac.x executable:',float(output.split()[2])/(1024*1024),'MB'
else :
  print 'can not determine'; sys.exit('exit!')

#### Process readelf ####

print 'analyzing individual objects...'

command = " readelf -s "+dirac_exe+" 2> /dev/null || echo 'can not be determined, readelf not present'"
process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
#p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True).stdout

#while 1:
#  line = p.readline()
#  print line

#for line in p.read():
#  print line

#output = process.stdout.read()
output = process.stdout.read().split()
print output[1]
for word in output:
  print word

#object_line=[[]]
#for line in output:
#for line in process.stdout.read().strip():
#for line in process.stdout.read().strip():
 #line = line.strip()
 #object_line.append(line)
#  print line

