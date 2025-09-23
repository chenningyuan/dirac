/* dirac_copyright_start */
/*
 *
 *     Copyright (c) by the authors of DIRAC.
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU Lesser General Public
 *     License version 2.1 as published by the Free Software Foundation.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *     Lesser General Public License for more details.
 *
 *     If a copy of the GNU LGPL v2.1 was not distributed with this
 *     code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
 */
/* dirac_copyright_end */

// File: gendet.c
#include <stdio.h>      // printf, ?
#include <string.h>     // strcat

#if defined (SYS_LINUX) || defined (SYS_AIX) || defined (SYS_DARWIN) || defined (SYS_WINDOWS)

#define cgendet     cgendet_

#define ibtabini    ibtabini_
#define fnbits      fnbits_
#define quit        quit_
#define getnbitscalls getnbitscalls_

#endif

#if defined (INT_STAR8)
#include <stdint.h>
typedef int64_t integer;
typedef unsigned long long int unsigned_integer;
#else
typedef int integer;
typedef unsigned int unsigned_integer;
#endif

/*
 * cgendet generates all n-bits integers with m bits set under
 * the symmetry constraints, occupation constraints, and M_k constraints
 * given by the functions gascip_{symok,gasok,mkok}.
 *
 * The code is written in C as fortran77 does not support recursive
 * functions.
 *
 */

extern void
getnbitscalls ( integer *fReset, long long int *nCN, long long int *nFN );

extern void
ibtabini ( unsigned_integer an[] );

extern void
quit ( char sz[], integer nLen );

extern void
fnbits ( integer *n , unsigned long long int *i, unsigned_integer anBitTable[] );

extern integer
cnbits ( unsigned long long int i, unsigned_integer anBitTable[] );

extern integer
cgasok( unsigned long long int prefix, integer nGAS, unsigned long long int anGASMask[], integer anGASSpc[] );

extern const char
*lli_to_b_format ( unsigned long long int x , integer bits) ;

extern void
cgendet1( unsigned long long int anStr[], integer nOrb, integer mElec,
	  integer *j, unsigned long long int prefix, integer pos,
	  integer nGAS, integer nGASTot[], unsigned long long int anGASMask[], integer anGASSpc[],
	  integer nMaxStr );

unsigned_integer anBitTable[ 256 ];
integer nElec, nOrb;
long long int nCNBitsCalls = 0;
long long int nFNBitsCalls = 0;


extern void
cgendet ( integer *pnElec, integer *pnAsht,
	  integer *pnStr, unsigned long long int anStr[],
	  integer *pnGAS, integer nGASTot[],
          unsigned long long int anGASMask[], integer anGASSpc[],
	  integer *pnMaxStr ) {

  /*
   *
   *  INTEGER   IEL, NASHT, NSTR, NGAS_DC, nGASTOT, nGASSpc(1:2), MAXSTR
   *  INTEGER*8 JSTR(1:MAXSTR), MSKGAS(1:NGAS)
   *  CALL CGENDET(IEL,NASHT,NSTR,JSTR,
   * &             NGAS_DC,nGASTOT,MSKGAS,nGASSpc,MAXSTR)
   *
   */

  const unsigned long long int i8_zero = 0;

  /* 
   * Initialize bit table
   */

  ibtabini( anBitTable );
  nElec = *pnElec;
  nOrb  = *pnAsht;
    
  *pnStr = 0;

     // printf ("\nHi from cgendet. IEL, NASHT, NGAS %i %i %i\n\n",nElec,nOrb,*pnGAS);

  cgendet1 ( anStr,  nOrb,     nElec,     pnStr,    i8_zero,    0,
	       *pnGAS, nGASTot, anGASMask, anGASSpc, *pnMaxStr );

     // printf ("   cgendet: number of strings %i\n", *pnStr);
     // fflush(stdout);

  /*
  integer i;
  for (i = 0; i < *pnStr; i++ )
    printf ("string %i: %lli\n", i, anStr[i]);
  */

}


extern void
cgendet1( unsigned long long int anStr[], integer n, integer m,
	  integer *j, unsigned long long int prefix, integer pos,
	  integer nGAS, integer nGASTot[],
          unsigned long long int anGASMask[], integer anGASSpc[],
	  integer nMaxStr ) {
/*
 * anStr  : vector of strings
 * n      : number of orbitals
 * m      : number of electrons
 * j      : number of strings
 * prefix : 0 on entry from cgendet; changes during recursive calls
 * pos    : 0 on entry from cgendet; increases during recursive calls
 * nGAS   : number of GAS spaces
 * nGASTot[i-1] : total number of orbitals in GAS spaces 1:i
 * anGASMask : mask for GAS spaces
 * anGASSpc  : min and max electrons each GAS space
 * nMaxStr : max number of strings
 */
 
   integer k1, k2;
   unsigned long long int pre;
   const unsigned long long int i8_one = 1;

//   printf ("entry cgendet1. n_o m_e j_str prefix pos %i %i %i %llu %i\n",n,m,*j,prefix,pos);

   if ( !m ) {  // m == 0
     if ( cgasok( prefix, nGAS, anGASMask, anGASSpc ) ) {
       if ( *j < nMaxStr ) {
	   anStr [ *j ] = prefix;
//            printf("GAS OK, bit string: %64s \n", lli_to_b_format(prefix,pos));
	   (*j)++;
       }
       else {
#ifdef INT_STAR8
	 printf ("*** ERROR in CGENDET1 *** Too many strings: %li! Allocate more memory in calling routine!\n", *j );
#else
	 printf ("*** ERROR in CGENDET1 *** Too many strings: %i! Allocate more memory in calling routine!\n", *j );
#endif
	 quit( "*** ERROR in CGENDET1 ***", 25 );
       }
     }
     else {
//       printf("GAS :(, bit string: %64s \n", lli_to_b_format(prefix,pos));
     }
   }
   else {
      
     for (k1 = 0; k1 < nGAS; k1++ ) {
       k2 = nGASTot[k1];
       if ( pos == k2 ) {
          if ( !cgasok( prefix, k1+1, anGASMask, anGASSpc ) ) {
//             printf("pos %i; NOT OK bit string: %64s \n", pos, lli_to_b_format(prefix,pos));
             return;
          }
       }
     }

     /* "1" + (n-1) bit number with (m-1) bits */

     pre = (i8_one<<pos) | prefix;
//     printf("+1 next bit string: %64s \n", lli_to_b_format(pre,pos+1));
     cgendet1( anStr, n - 1, m - 1, j, pre, pos+1,
		   nGAS, nGASTot, anGASMask, anGASSpc, nMaxStr );

     /* "0" + (n-1) bit number with m bits */

     if ( n > m ) {
//       printf("+0 next bit string: %64s \n", lli_to_b_format(prefix,pos+1));
       cgendet1( anStr, n - 1, m, j, prefix, pos+1,
		nGAS, nGASTot, anGASMask, anGASSpc, nMaxStr );
     }
   }

}

extern integer
cgasok( unsigned long long int prefix, integer nGAS, unsigned long long int anGASMask[], integer anGASSpc[] ) {

  long i, j, n;

  // printf("cgasok debug prefix string: %64s \n", lli_to_b_format(prefix,64));
  for (i = 0, j = 0; i < nGAS; i++, j+=2 ) {

    // printf("cgasok debug i, j     : %lli %lli \n",i,j);
    // printf("cgasok debug anGASMask: %64s \n", lli_to_b_format(anGASMask[ i ],64));
    // printf("cgasok debug bit and  : %64s \n", lli_to_b_format(prefix & anGASMask[ i ],64));
    n = cnbits ( prefix & anGASMask[ i ], anBitTable );
    // printf("cgasok debug cnbits   : %lli \n", n);
    // printf("cgasok debug anGASSpc : %lli %lli\n", anGASSpc[ j ], anGASSpc[ j + 1 ] );

    if ( ( n < anGASSpc[ j ] ) || ( n > anGASSpc[ j + 1 ] ) )
      return 0;

  }

  return 1;

}

extern integer
cnbits ( unsigned long long int i, unsigned_integer anBitTable[] ) {

  unsigned char *j;

  nCNBitsCalls++;

  j = (unsigned char *) &i;

  return anBitTable[ *j ] + anBitTable[ *(j+1) ] +
    anBitTable [ *(j+2) ] + anBitTable[ *(j+3) ] +
    anBitTable [ *(j+4) ] + anBitTable[ *(j+5) ] +
    anBitTable [ *(j+6) ] + anBitTable[ *(j+7) ];

}

extern void
fnbits ( integer *n, unsigned long long int *i, unsigned_integer anBitTable[] ) {

  unsigned char *j;

  nFNBitsCalls++;

  j = (unsigned char *) i;

  *n =   anBitTable[ *j ] + anBitTable[ *(j+1) ] +
    anBitTable [ *(j+2) ] + anBitTable[ *(j+3) ] +
    anBitTable [ *(j+4) ] + anBitTable[ *(j+5) ] +
    anBitTable [ *(j+6) ] + anBitTable[ *(j+7) ];

}

extern void
getnbitscalls ( integer *fReset, long long int *nCN, long long int *nFN ) {

  *nCN = nCNBitsCalls;
  *nFN = nFNBitsCalls;

  if ( *fReset ) {
    nCNBitsCalls = 0;
    nFNBitsCalls = 0;
  }

}

const char *lli_to_b_format ( unsigned long long int x , integer bits)
{
    // Oct 2010 Hans Joergen Aa. Jensen
    // based of code by EviTeach on
    // http://stackoverflow.com/questions/111928/is-there-a-printf-converter-to-print-in-binary-format

    static char b[65];
    b[0] = '\0';

    if ( (bits > 64) || (bits <= 0) ) {
       strcat(b, "\n *** ERROR Illegal number of bits! ***\n");
       return b;
    }

    unsigned long long int z, zmax;
    zmax = 1; zmax  <<= bits-1;

    for (z = zmax; z > 0; z >>= 1)
    {
        strcat(b, ((x & z) == z) ? "1" : "0");
    }

    return b;
}
// end of gendet.c
