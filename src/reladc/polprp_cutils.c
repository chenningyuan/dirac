#include <stdio.h>

/*
 * these are fields defined in the polprp_cfg Module and
 * they are visible externally. The data are located
 * either in the .DATA or .BSS segment and their start
 * address can therefore be determined and the contents
 * accessed via our C utility. 
 *
 */

extern char ppinfo_code;    // character array made available in polprp_cfg module

int icode_activated_(char * icode){

   int i;
   char * a;

   printf("\n ---entering icode_activated---\n");
   a =  &ppinfo_code;
   printf("ppinfo string: %s\n",a);
   printf("\n ---end---\n");
   return 0;
}
