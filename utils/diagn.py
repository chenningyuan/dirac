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
#
#   Written by Miro Ilias, Banska Bystrica, 2010
#
#================================================================================


import string, re, sys, os, subprocess

#
def hex2dec(s):
# return the integer value of a hexadecimal string s
  if  re.match("0x[0-9a-fA-F]+$", s) or re.match("[0-9]+$", s):
   return int(s, 16)
  else :
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
  print '\nTotal size of all static allocations in the dirac.x executable:',float(output.split()[2])/(1024*1024),'MB'
else :
  print 'can not determine'; sys.exit('exit!')

#### Process readelf ####

print 'analyzing individual objects...'

command = " readelf -s "+dirac_exe+" 2> /dev/null || echo 'can not be determined, readelf not present'"
process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

output = process.stdout.read().splitlines(False)

### substitute strings-hex numbers by converted decimal numbers
for line in output:
   line=line.strip().split()
   numstr=string.join(line[2:3]);#  print 'numstr=',numstr
   decstr=hex2dec(numstr);        #print 'decstr=',decstr

   #print decstr,line[2:3],line
   # replace hexadecimal value with decimal value (string)
   line.insert(3,str(decstr))
   #print line


data_as_lists = [ line.strip().split() for line in output ]     ## convert to a list of lists

#data_as_lists = [ for line in output:  line.strip().split(); print line[2:3] ]     ## convert to a list of lists


#data_as_lists.sort(key=operator.itemgetter(-1))  ## sort on last element
#data_as_lists.sort(key=lambda (3))  ## sort somehow...
#data_as_lists.sort(key='3')  ## sort somehow...
##data_as_lists.sort(reverse=True)  ## sort somehow...
##data_as_lists.sort(key=str,reverse=True)  ## sort somehow...

### works....
#data_as_lists.sort(key=lambda data_as_lists:data_as_lists[3:3], reverse=True)  ## sort somehow...
#data_as_lists.sort(key=lambda data_as_lists:data_as_lists[1:3], reverse=True)  ## sort somehow...
#data_as_lists.sort(key=lambda data_as_lists:data_as_lists[-5: :3], reverse=True)  ## sort somehow...
#data_as_lists.sort(key=lambda data_as_lists:data_as_lists[:3 :-3 ], reverse=True)  ## sort somehow...
#data_as_lists.sort(key=lambda data_as_lists:data_as_lists[:-3 :-2 ], reverse=True)  ## sort somehow...
#data_as_lists.sort(key=lambda data_as_lists:data_as_lists[:+5 :-2 ], reverse=False)  ## sort somehow...
#data_as_lists.sort(key=lambda data_as_lists:data_as_lists[3:8 ], reverse=False)  ## sort somehow...
data_as_lists.sort(key=lambda data_as_lists:data_as_lists[2:3], reverse=False)  ## sort somehow...


#data_as_lists.sort(key=lambda data_as_lists:data_as_lists[2::3], reverse=False)  ## sort somehow...

#data_as_lists[:-2]

for el in data_as_lists : ## print each item in the list
  #print el.split()[1]
  # z=el
  #print len(el),el
  #print len(el),z[-2: :-1]
  #print el[-1: :-5]
  if (len(el) >= 8) :
    #print el[2:3], el[7:8]
    #if e[2:3].isdigit():
    #  print e[2:3]
    print el[2:3], el[7:8]
     #print el


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

