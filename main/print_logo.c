/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Print logo if master node
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdio.h>

#include "print_logo.head"
 ;


#include "print_logo.head"

{
#include "title.h"
#include "logo.h"
#include "subtitle.h"

   if (mytid != 0) return;

   printf("%s",strtitle);
   printf("%s",strlogo);
   printf("%s",strsubtitle);

   sysinf();
   
}
