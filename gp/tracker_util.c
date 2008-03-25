/*
 * $Id$
 *
 * tracker_util.c : functions to retrieve information from memory usage (addresses, amount of free memory)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined (SYS_LINUX)
#include <sys/sysinfo.h>
#endif /* SYS_LINUX */

#define BUFFERSIZE (512)

void get_memory_address_(void *variable, int *address);
void get_memory_info_(char *fstring, int *flength);

// variable: is the real 
// address_string: is a character string with the memory address of variable

void get_memory_address_(void *variable, int *address)
{
   char bufvar[BUFFERSIZE]; 

   snprintf (bufvar,(BUFFERSIZE-1),"%i",variable);

// address_string = malloc(strlen(bufvar)*sizeof(char)); 
// strncpy(address_string,bufvar,strlen(bufvar));

   *address = atoi(bufvar); 
}


// reads memory information. note, sysinfo is linux-specific 

void get_memory_info_(char *fstring, int *flength)
{
   unsigned long total_bytes = 1, total_kbytes = 1;
   size_t i, len;
   int status;

#if defined (SYS_LINUX)
   struct sysinfo info;
   
   status = sysinfo(&info);
   if (status == 0)
   {
      total_bytes  = info.freeram * info.mem_unit;
      total_kbytes = total_bytes / 1024;
   }
#endif /* SYS_LINUX */

   snprintf(fstring,(size_t)(*flength),"%li Kb",total_kbytes);
   len = strlen(fstring);
   
   if ( len <= *flength )
   {
      for ( i = len; i < *flength; i++)
         fstring[i] = ' ';
      fstring[i] = '\0';
   }
}
