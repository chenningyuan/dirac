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

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  if (argc > 1)
    {
      printf("Generate a C header file for Fortran->C calling\n");
      return EXIT_SUCCESS;
    }
  printf("#ifndef FORTRAN_INTERFACE_H\n");
  printf("#define FORTRAN_INTERFACE_H\n");

  printf("typedef ");
  if (sizeof(short int) == 8)
    printf("short int");
  if (sizeof(int) == 8)
    printf("int");
  else if (sizeof(long int) == 8)
    printf("long int");
  else if (sizeof(long long int) == 8)
    printf("long long int");
  else
    return EXIT_FAILURE;
  printf(" fortran_int8_t;\n");

  printf("typedef ");
  if (sizeof(short int) == 4)
    printf("short int");
  if (sizeof(int) == 4)
    printf("int");
  else if (sizeof(long int) == 4)
    printf("long int");
  else if (sizeof(long long int) == 4)
    printf("long long int");
  else
    return EXIT_FAILURE;
  printf(" fortran_int4_t;\n");

#ifdef INT_STAR8
  printf("typedef fortran_int8_t fortran_int_t;\n");
#else
  printf("typedef fortran_int4_t fortran_int_t;\n");
#endif
  printf("#endif\n");
  return EXIT_SUCCESS;
}
