/*
 * $Id$
 *
 * tracing.c : functions to retrieve information from the execution stack and process them into fortran-friendly data
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined (SYS_LINUX)
#include <execinfo.h>
#endif

#define MAXSIZE_SYMBOL (300)
#define MAXDEPTH_TRACE  (30)

enum allowedActions { printAllToStdout = 1, printAllocCaller, printPrettyTrace, prinRawTrace };
enum retrieveSymbolKinds {symbolName = 1, symbolRetAddr};

// public (callable from fortran) functions. 

void get_traceback_info_(enum allowedActions*, char *fstring, int *fsize);

// internal (not meant to be called from fortran) function

void retrieve_symbol(enum retrieveSymbolKinds, char *result, char *backtrace_string);
void print_rawtrace_to_stdout(size_t size, char **strings);
void get_name_from_caller(size_t size, char **strings);
void print_trace_to_string(char *symbol_in, char *symbol_out, int length);
void print_prettytrace_to_stdout (size_t size, char **strings);

// interface to fortran

//
// get_traceback_info: will perform the traceback and get whatever information we want and copy it to the result string or prints to stdout
//
// variables:
//    action  : specifies what to retieve and what to do with it
//    result  : when passing information back to the fortran program, *result is the fortran array holding the data
//    fsize   : when passing information back to the fortran program, *fsize holds the size of the "result" array in fortran
//

void get_traceback_info_(enum allowedActions *action, char *fstring, int *fsize)
{
  void *array[MAXDEPTH_TRACE];
  size_t size = 0;
  char **strings = NULL;
  size_t i=3, j, k = 0;
  char symbol_name[MAXSIZE_SYMBOL];

// only support to linux for now. for other systems the traceback is yet to be implemented
// mac osx 10.5 should have the same interface
#if defined (SYS_LINUX)
// calls to gnu libc routines
   size    = backtrace (array, MAXDEPTH_TRACE);
   strings = backtrace_symbols (array, size);
#endif /* SYS_LINUX */

   if ( strings != NULL ) 
   {
// now we branch off for the different pieces of the trace we want 
// to use
      if ( *action == printAllToStdout ) 
      {
         retrieve_symbol(symbolRetAddr,symbol_name,strings[i]);
//         print_full_trace(size,strings,fstring);
      }
      else if ( *action == printAllocCaller )
      {
         retrieve_symbol(symbolName,symbol_name,strings[i]);
         print_trace_to_string(symbol_name,fstring,*fsize);
      }
      else if ( *action == printPrettyTrace )
      {
         print_prettytrace_to_stdout(size,strings);
      }
      else if ( *action == prinRawTrace )
      {
         print_rawtrace_to_stdout(size,strings);
      }

      free (strings);
   }
}


// internal functions


//
// print_rawtrace_to_stdout: gets the information from the trace and prints it as it is to stdout
//
// size    : depth of the trace
// strings : information from the trace

void print_rawtrace_to_stdout(size_t size, char **strings)
{
  size_t i;

  printf ("\nTraceback information:\n\n\tObtained %zd stack frames.\n", size);

  for (i = 0; i < size; i++)
     printf ("\t\t%100s\n", strings[i]);

   printf ("\n");
}


/* Obtain a backtrace and print it to stdout. */
void print_prettytrace_to_stdout (size_t size, char **strings)
{
   size_t i;
   char symbol_name[MAXSIZE_SYMBOL];
 
   // minus five below is 2 from the glibc, one from this function 
   // and 2 from the fortran side (error handler and the allocation
   
   printf ("\n Traceback from allocation call: %zd frames.\n\n",(size-3));

   for (i = 3; i < size; i++)
   {
      retrieve_symbol(symbolName,symbol_name,strings[i]);
      
      if (strncmp(symbol_name,"",1))
      {
         printf("%3d   function name:   %30s\n",(size-i),symbol_name);
      }
      else
      {
         printf("%3d   function name:   Unknown\n",(size-i));
      }
   }
   printf("\n");

}



// print_trace_to_string:
//
// symbol_in  :  holds a C string with the already processed (no parenthesis or brackets) trace symbol (function name+entry address or return adress)
// symbol_out : 
// length     : integer that 

void print_trace_to_string(char *symbol_in, char *symbol_out, int length)
{
   size_t i,j=0,n;
   
   n = strlen(symbol_in);

   for (i = 0; i < n; i++)
   {
      if (( symbol_in[i] != ' ' )&&( i < length ))
      {
         symbol_out[j] = symbol_in[i];
         j++;
      }
   }
   for (i = j ; i < length; i++)
      symbol_out[i] = ' ';

   symbol_out[length-1] = '\0';
}


void retrieve_symbol(enum retrieveSymbolKinds kind, char *result, char *backtrace_string)
{
  size_t j, k;
  int symbol_name_start, symbol_name_end;
  char element;

  symbol_name_start = -1;
  symbol_name_end   = -1;
  k = 0;

  for ( j = 0 ; ; j++)
     {
        element = *(backtrace_string+j);

        if ( element == '\0' )
        {
           result[k] = element;
           break;
        }

        if ( ((element == '(')&&(kind == symbolName))    ||
             ((element == '[')&&(kind == symbolRetAddr))  )
        {
           symbol_name_start = j+1;
        }
        else if ( ((element == ')')&&(kind == symbolName))   ||
                  ((element == ']')&&(kind == symbolRetAddr)) )
        {
           symbol_name_end   = j-1;
           symbol_name_start = -1;
        }

        if ( (symbol_name_end == -1)&&( j >= symbol_name_start ) )
        {
           result[k] = element;
           k++;
        }
     }
}

