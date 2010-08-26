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


import re, sys, os

dirac_exe="dirac.x"

if not os.path.isfile(dirac_exe):
   print 'executable dirac.x does not exist! Exit'
   sys.exit(-1)

print 'analyzing static allocation in the dirac.x executable'
#cmd = 'size '+dirac_exe+' 2> /dev/null | grep dirac.x | awk "{printf $$3/(1024*1024)}"; echo " MB" || echo "can not be determined"'
cmd = 'size '+dirac_exe+' 2> /dev/null | grep dirac.x || echo "can not be determined"'
status = os.system(cmd)
#status = os.popen(cmd)
#print status


