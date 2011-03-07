/* dirac_copyright_start */
/*
 *     Copyright (c) 2010 by the authors of DIRAC.
 *     All Rights Reserved.
 *
 *     This source code is part of the DIRAC program package.
 *     It is provided under a written license and may be used,
 *     copied, transmitted, or stored only in accordance to the
 *     conditions of that written license.
 *
 *     In particular, no part of the source code or compiled modules may
 *     be distributed outside the research group of the license holder.
 *     This means also that persons (e.g. post-docs) leaving the research
 *     group of the license holder may not take any part of Dirac,
 *     including modified files, with him/her, unless that person has
 *     obtained his/her own license.
 *
 *     For information on how to get a license, as well as the
 *     author list and the complete list of contributors to the
 *     DIRAC program, see: http://dirac.chem.vu.nl
 */
/* dirac_copyright_end */

#include <stdio.h>

#if defined (SYS_IRIX) || defined (SYS_LINUX) || defined (SYS_DEC) || defined (SYS_AIX) || defined (SYS_HPUX) || defined (SYS_DARWIN) || defined (SYS_WINDOWS)

#define cgendet     cgendet_

#define ibtabini    ibtabini_
#define fnbits      fnbits_
#define quit        quit_
#define getnbitscalls getnbitscalls_

#endif

#if defined (INT_STAR8)
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
getnbitscalls ( long *fReset, long long int *nCN, long long int *nFN );

extern void
ibtabini ( long an[] );

extern void
quit ( char sz[], long nLen );

extern void
fnbits ( long *n , long *i, long anBitTable[] );

extern long
cnbits ( long i, long anBitTable[] );

extern long
cgasok( long prefix, long nGAS, long anGASMask[], long anGASSpc[] );

extern void
cgendet1( long anStr[], long nElec, long mElec,
	  long *j, long prefix, long pos,
	  long nGAS,
	  long anGASMask[],
	  long anGASSpc[],
	  long nMaxStr );

long anBitTable[ 256 ];
long long int nCNBitsCalls = 0;
long long int nFNBitsCalls = 0;


extern void
cgendet ( long *pnElec, long *pnAsht,
	  long *pnStr, long anStr[],
	  long *pnGAS,
	  long anGASMask[],
	  long anGASSpc[],
	  long *pnMaxStr ) {

  long i;

  /* 
   * Initialize bit table
   */

  ibtabini( anBitTable );
    
  *pnStr = 0;

  cgendet1 ( anStr, *pnAsht, *pnElec, pnStr, 0, 0,
	     *pnGAS, anGASMask, anGASSpc, *pnMaxStr);

  /* printf ("   cgendet: number of strings %i\n", *pnStr); */

  /*
  for (i = 0; i < *pnStr; i++ )
    printf ("string %i: %i\n", i, anStr[i]);
  */

}


extern void
cgendet1( long anStr[], long n, long m,
	  long *j, long prefix, long pos,
	  long nGAS,
	  long anGASMask[],
	  long anGASSpc[],
	  long nMaxStr ) {
 
   long pre;

   if ( !m ) {
     if ( cgasok( prefix, nGAS, anGASMask, anGASSpc ) ) {
       if ( *j < nMaxStr ) {
	 anStr [ *j ] = prefix;
	 (*j)++;
       }
       else {
	 printf ("*** ERROR in CGENDET1 *** "
		 "Too many strings: %i! "
		 "Allocate more memory!\n", (int)*j );
	 quit( "*** ERROR in CGENDET1 ***", 26 );
       }
     }
   }
   else {

      /* "1" + (n-1) bit number with (m-1) bits */

      pre = (1<<pos) | prefix;
      cgendet1( anStr, n - 1, m - 1, j, pre, pos+1,
		nGAS, anGASMask, anGASSpc, nMaxStr );

      /* "0" + (n-1) bit number with m bits */

      if ( n > m )
         cgendet1( anStr, n - 1, m, j, prefix, pos+1,
		nGAS, anGASMask, anGASSpc, nMaxStr );

   }

}


extern long
cgasok( long prefix, long nGAS, long anGASMask[], long anGASSpc[] ) {

  long i, j, n;

  for (i = 0, j = 0; i < nGAS; i++, j+=2 ) {
    
    n = cnbits ( prefix & anGASMask[ i ], anBitTable );

    if ( ( n < anGASSpc[ j ] ) || ( n > anGASSpc[ j + 1 ] ) )
      return 0;

  }

  return 1;

}

extern long
cnbits ( long i, long anBitTable[] ) {

  unsigned char *j;

  nCNBitsCalls++;

  j = (unsigned char *) &i;

  return anBitTable[ *j ] + anBitTable[ *(j + 1) ] +
    anBitTable [ *(j+2) ] + anBitTable[ *(j + 3) ];

}

extern void
fnbits ( long *n, long *i, long anBitTable[] ) {

  unsigned char *j;

  nFNBitsCalls++;

  j = (unsigned char *) i;

  *n = anBitTable[ *j ] + anBitTable[ *(j + 1) ] +
    anBitTable [ *(j+2) ] + anBitTable[ *(j + 3) ];

}

extern void
getnbitscalls ( long *fReset, long long int *nCN, long long int *nFN ) {

  *nCN = nCNBitsCalls;
  *nFN = nFNBitsCalls;

  if ( *fReset ) {
    nCNBitsCalls = 0;
    nFNBitsCalls = 0;
  }

}
#else
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
getnbitscalls ( int *fReset, long long int *nCN, long long int *nFN );

extern void
ibtabini ( int an[] );

extern void
quit ( char sz[], int nLen );

extern void
fnbits ( int *n , int *i, int anBitTable[] );

extern int
cnbits ( int i, int anBitTable[] );

extern int
cgasok( int prefix, int nGAS, int anGASMask[], int anGASSpc[] );

extern void
cgendet1( int anStr[], int nElec, int mElec,
	  int *j, int prefix, int pos,
	  int nGAS,
	  int anGASMask[],
	  int anGASSpc[],
	  int nMaxStr );

int anBitTable[ 256 ];
long long int nCNBitsCalls = 0;
long long int nFNBitsCalls = 0;


extern void
cgendet ( int *pnElec, int *pnAsht,
	  int *pnStr, int anStr[],
	  int *pnGAS,
	  int anGASMask[],
	  int anGASSpc[],
	  int *pnMaxStr ) {

  int i;

  /* 
   * Initialize bit table
   */

  ibtabini( anBitTable );
    
  *pnStr = 0;

  cgendet1 ( anStr, *pnAsht, *pnElec, pnStr, 0, 0,
	     *pnGAS, anGASMask, anGASSpc, *pnMaxStr);

  /* printf ("   cgendet: number of strings %i\n", *pnStr); */

  /*
  for (i = 0; i < *pnStr; i++ )
    printf ("string %i: %i\n", i, anStr[i]);
  */

}


extern void
cgendet1( int anStr[], int n, int m,
	  int *j, int prefix, int pos,
	  int nGAS,
	  int anGASMask[],
	  int anGASSpc[],
	  int nMaxStr ) {
 
   int pre;

   if ( !m ) {
     if ( cgasok( prefix, nGAS, anGASMask, anGASSpc ) ) {
       if ( *j < nMaxStr ) {
	 anStr [ *j ] = prefix;
	 (*j)++;
       }
       else {
	 printf ("*** ERROR in CGENDET1 *** "
		 "Too many strings: %i! "
		 "Allocate more memory!\n", *j );
	 quit( "*** ERROR in CGENDET1 ***", 26 );
       }
     }
   }
   else {

      /* "1" + (n-1) bit number with (m-1) bits */

      pre = (1<<pos) | prefix;
      cgendet1( anStr, n - 1, m - 1, j, pre, pos+1,
		nGAS, anGASMask, anGASSpc, nMaxStr );

      /* "0" + (n-1) bit number with m bits */

      if ( n > m )
         cgendet1( anStr, n - 1, m, j, prefix, pos+1,
		nGAS, anGASMask, anGASSpc, nMaxStr );

   }

}


extern int
cgasok( int prefix, int nGAS, int anGASMask[], int anGASSpc[] ) {

  int i, j, n;

  for (i = 0, j = 0; i < nGAS; i++, j+=2 ) {
    
    n = cnbits ( prefix & anGASMask[ i ], anBitTable );

    if ( ( n < anGASSpc[ j ] ) || ( n > anGASSpc[ j + 1 ] ) )
      return 0;

  }

  return 1;

}

extern int
cnbits ( int i, int anBitTable[] ) {

  unsigned char *j;

  nCNBitsCalls++;

  j = (unsigned char *) &i;

  return anBitTable[ *j ] + anBitTable[ *(j + 1) ] +
    anBitTable [ *(j+2) ] + anBitTable[ *(j + 3) ];

}

extern void
fnbits ( int *n, int *i, int anBitTable[] ) {

  unsigned char *j;

  nFNBitsCalls++;

  j = (unsigned char *) i;

  *n = anBitTable[ *j ] + anBitTable[ *(j + 1) ] +
    anBitTable [ *(j+2) ] + anBitTable[ *(j + 3) ];

}

extern void
getnbitscalls ( int *fReset, long long int *nCN, long long int *nFN ) {

  *nCN = nCNBitsCalls;
  *nFN = nFNBitsCalls;

  if ( *fReset ) {
    nCNBitsCalls = 0;
    nFNBitsCalls = 0;
  }

}
#endif
