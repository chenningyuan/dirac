/*
 * $Id$
 *
 * tracker_util.c : functions to retrieve information from memory usage (addresses, amount of free memory)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "system_flavor.h"

#if (SYS_FLAVOR == LINUX)
#define LINUX_IFACE_MEMORY
#endif

/* using linux's interface to memory */
#if defined (LINUX_IFACE_MEMORY)
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif /* LINUX_IFACE_MEMORY */

#define BUFFERSIZE (512)
#define TYPICAL_FCHAR_LENGTH (132)

void get_memory_address_(void *variable, long int *address);
void get_memory_info_(char *fstring, int *flength);
/* void get_free_mwords_(int *freemw); */
void get_free_mwords_(size_t *freemw);

/* variable: is the real 
   address_string: is a character string with the memory address of variable
                   question is whether or not this should always be long int in all archs.. we need to be able to hold up to 64 bit ints */

void get_memory_address_(void *variable, long int *address)
{
   char bufvar[BUFFERSIZE]; 

   snprintf (bufvar,(BUFFERSIZE-1),"%li",variable);

   *address = atol(bufvar); 
}


/* reads memory information. note, sysinfo is linux-specific */

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
   
   if ( *flength > TYPICAL_FCHAR_LENGTH )
      printf ("warning: possible problem with fortran character array size\n"
              "         value passed to get_memory_info() is %d\n"
              "         typical value in the code is however %d\n",
              *flength,TYPICAL_FCHAR_LENGTH);

   if ( len < *flength )
   {
      for ( i = len; i < *flength; i++)
         fstring[i] = ' ';
      fstring[i] = '\0';
   }
}



void get_free_mwords_(size_t *freemw)
{
   unsigned long total_bytes = 1, total_kbytes = 0;
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
   *freemw = (size_t)(total_kbytes/(1024*8)); 
/*
   printf ("Debug: getFreeMwords() found %d (%ld,%ld) mWords avaiable\n",*freemw,*freemw,total_kbytes);
   catfile_();
*/
#else
   *freemw = (size_t)(total_kbytes);
/*
   printf ("Debug: getFreeMwords() found %d (%ld) mWords avaiable\n",*freemw,*freemw);
*/
#endif /* LINUX_IFACE_MEMORY */
}



void catfile_(void)
{
   (void)system("cat /proc/meminfo"); 
}

void showlimits_(void)
{
   int as_query, data_query, memlock_query;

#if defined (LINUX_IFACE_MEMORY)
   struct rlimit lim_as, lim_data, lim_memlock;

/*
 * from the manpage of getrlimit(2):
 *
 *
       getrlimit() and setrlimit() get and set resource limits respectively.
       Each resource has an associated soft and hard limit, as defined by 
       the rlimit structure (the rlim argument to both getrlimit() and 
       setrlimit()):

            struct rlimit {
                rlim_t rlim_cur;     Soft limit   
                rlim_t rlim_max;     Hard limit (ceiling for rlim_cur)   
            };

       The  soft limit is the value that the kernel enforces for the 
       corresponding resource.  The hard limit acts as a ceiling for 
       the soft limit: an unprivileged process may only set its soft 
       limit to a value in the range from 0 up to the hard limit, and 
       (irreversibly) lower its hard limit.  A privileged process 
       (under Linux: one with the CAP_SYS_RESOURCE capability) may make 
       arbitrary changes to either limit value.

       The value RLIM_INFINITY denotes no limit on a resource (both 
       in the structure returned by getrlimit() and in the structure 
       passed to setrlimit()).

       resource must be one of:
       RLIMIT_AS
              The maximum size of the process's virtual memory (address 
              space) in bytes.  This limit affects calls to brk(2), 
              mmap(2) and mremap(2), which fail with the error ENOMEM upon 
              exceeding this limit.  Also  automatic  stack expansion  will
              fail (and generate a SIGSEGV that kills the process if no 
              alternate stack has been made available via sigaltstack(2)).
              Since the value is a long, on machines with a 32-bit long 
              either this limit is at most 2 GiB, or this resource is unlimited.

       RLIMIT_DATA
              The  maximum size of the process's data segment (initialized 
              data, uninitialized data, and heap).  This limit affects calls 
              to brk() and sbrk(), which fail with the error ENOMEM upon 
              encountering the soft limit of this resource.

       RLIMIT_MEMLOCK
              The  maximum  number  of  bytes of memory that may be locked 
              into RAM.  In effect this limit is rounded down to the nearest
              multiple of the system page size.  This limit affects mlock(2) 
              and mlockall(2) and the mmap(2) MAP_LOCKED operation.  Since 
              Linux 2.6.9 it also affects the shmctl(2) SHM_LOCK operation, 
              where it sets a maximum on the total bytes in shared memory 
              segments (see shmget(2)) that may be locked by the real user 
              ID  of the  calling process.  The shmctl(2) SHM_LOCK locks are 
              accounted for separately from the per-process memory locks 
              established by mlock(2), mlockall(2), and mmap(2) MAP_LOCKED; 
              a process can lock bytes up to this limit in each of these 
              two categories.  In Linux kernels before 2.6.9, this limit 
              controlled the amount of memory that could be locked by a 
              privileged process.  Since Linux 2.6.9, no limits are placed on 
              the amount of memory that a privileged process may lock, and 
              this limit instead governs the amount of memory that an 
              unprivileged process may lock.
*
*/

   if ( (as_query = getrlimit(RLIMIT_AS,&lim_as)) != 0 )
      printf (" Error querying address space limits!\n");

   if ( (data_query = getrlimit(RLIMIT_DATA,&lim_data)) != 0 )
      printf (" Error querying data segment limits!\n"); 

   if ( (memlock_query = getrlimit(RLIMIT_MEMLOCK,&lim_memlock)) != 0 )
      printf (" Error querying memory lock limits!\n"); 

   if ( as_query == 0 && data_query == 0 && memlock_query == 0 )
   {
      printf(" Querying system limits in place for this process:\n"
             "   maximum size of address space (soft,hard): %lli, %lli, %li, %li, \n"
             "   maximum size of data segment  (soft,hard): %lli, %lli, %li, %li, \n"
             "   maximum size of locked memory (soft,hard): %lli, %lli, %li, %li, \n\n",
             lim_as.rlim_cur, lim_as.rlim_max, lim_as.rlim_cur, lim_as.rlim_max,
             lim_data.rlim_cur, lim_data.rlim_max, lim_data.rlim_cur, lim_data.rlim_max,
             lim_memlock.rlim_cur, lim_memlock.rlim_max, lim_memlock.rlim_cur, lim_memlock.rlim_max);
    
      if (lim_memlock.rlim_cur != -1)
      {
         printf (" Trying to reset memory lock...\n");
        
         lim_memlock.rlim_cur = 3000000000; /* about 2.8 gb */

         if ( (memlock_query = setrlimit(RLIMIT_MEMLOCK,&lim_memlock)) != 0 )
            printf (" Error setting memory lock limits!\n");
         else 
         {
            memlock_query = getrlimit(RLIMIT_MEMLOCK,&lim_memlock);
            printf (" Successfully (%d) reset memory to %li\n",memlock_query,lim_memlock.rlim_cur);
         }
      }
   }

#else
   printf (" Querying of system limits not implemented for this OS\n");
#endif
}
