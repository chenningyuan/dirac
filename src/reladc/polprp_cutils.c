#include <stdio.h>

/*
 * these are fields defined in the polprp_cfg Module and
 * they are visible externally. The data are located
 * either in the .DATA or .BSS segment and their start
 * address can therefore be determined and the contents
 * accessed via our C utility. 
 *
 */

extern int __polprp_cfg_MOD_polprp_printlev;    // integer value
extern int __polprp_cfg_MOD_polprp_inforeqest;  // integer field of length 32

int icode_activated_(int * icode){

   int i;
   int * a;

   a =  &__polprp_cfg_MOD_polprp_inforeqest;
   for(i=0;i<32;i++){
     if(*icode == *(a+i)) return 1;
   }
   return 0;
}
