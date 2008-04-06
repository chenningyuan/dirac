/*
 * $Id$
 *
 * tracker_util.c : functions to retrieve information from memory usage (addresses, amount of free memory)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <system_flavor.h>

#if (SYS_FLAVOR == LINUX)
#define LINUX_IFACE_MEMORY
#endif

// using linux's interface to memory
#if defined (LINUX_IFACE_MEMORY)
#include <sys/sysinfo.h>
#endif /* LINUX_IFACE_MEMORY */

#define BUFFERSIZE (512)

void get_memory_address_(void *variable, long int *address);
void get_memory_info_(char *fstring, int *flength);

// variable: is the real 
// address_string: is a character string with the memory address of variable
//                 question is whether or not this should always be long int in all archs.. we need to be able to hold up to 64 bit ints

void get_memory_address_(void *variable, long int *address)
{
   char bufvar[BUFFERSIZE]; 

   snprintf (bufvar,(BUFFERSIZE-1),"%li",variable);

   *address = atol(bufvar); 
}


// reads memory information. note, sysinfo is linux-specific 

void get_memory_info_(char *fstring, int *flength)
{
   unsigned long total_bytes = 1, total_kbytes = 1;
   size_t i, len;
   int status;

#if defined (LINUX_IFACE_MEMORY)
   struct sysinfo info;
   
   status = sysinfo(&info);
   if (status == 0)
   {
      total_bytes  = info.freeram * info.mem_unit;
      total_kbytes = total_bytes / 1024;
   }
#endif /* LINUX_IFACE_MEMORY */

   snprintf(fstring,(size_t)(*flength),"%li Kb",total_kbytes);
   len = strlen(fstring);
   
   if ( len <= *flength )
   {
      for ( i = len; i < *flength; i++)
         fstring[i] = ' ';
      fstring[i] = '\0';
   }
}
