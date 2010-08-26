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

#output=[]
#output = process.stdout.read()
output = process.stdout.read().splitlines(False)
#output = process.stdout.readline().split()
#print output

#output = process.stdout
#print output.read()
#print output.splitlines(True)

#for line in output:
#  print line.split()


data_as_lists = [ line.strip().split() for line in output ]     ## convert to a list of lists
#data_as_lists.sort(key=operator.itemgetter(-1))  ## sort on last element
#data_as_lists.sort(key=lambda (3))  ## sort somehow...
#data_as_lists.sort(key='3')  ## sort somehow...
data_as_lists.sort()  ## sort somehow...
for el in data_as_lists : ## print each item in the list
  #print el.split()[1]
  z=el
  print len(el),el
  #print z[2]


#print output[1],output[2],output[3],output[4],output[5],output[6],output[7],output[21]
#nl=(len(output)-21)/8
#print (len(output)-21)/8
#z=21
#for i in range(nl):
#  for j in range():
#    print 'i=',i,'***',output[i+21]

#for word in output:
#  print word

#object_line=[[]]
#for line in output:
#for line in process.stdout.read().strip():
#for line in process.stdout.read().strip():
 #line = line.strip()
 #object_line.append(line)
#  print line

