/* dirac_copyright_start */
/*
 *     Copyright (c) 2013 by the authors of DIRAC.
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
 *     DIRAC program, see: http://www.diracprogram.org
 */
/* dirac_copyright_end */

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

  if (sizeof(short int) == 8)
  {
    printf("#define FORTRAN_INT8_CONVERSION  \"%%i\" \n");
    printf("#define FORTRAN_STR_TO_INT8(a) atoi(a) \n");
    printf("#define FORTRAN_STRx_TO_INT8(a) (short int)strtol(a,NULL,16) \n");
    printf("typedef short int");
  }
  if (sizeof(int) == 8)
  {
    printf("#define FORTRAN_INT8_CONVERSION  \"%%i\" \n");
    printf("#define FORTRAN_STR_TO_INT8(a) atoi(a) \n");
    printf("#define FORTRAN_STRx_TO_INT8(a) (int)strtol(a,NULL,16) \n");
    printf("typedef int");
  }
  else if (sizeof(long int) == 8)
  {
    printf("#define FORTRAN_INT8_CONVERSION  \"%%li\" \n");
    printf("#define FORTRAN_STR_TO_INT8(a) atol(a) \n");
    printf("#define FORTRAN_STRx_TO_INT8(a) strtol(a,NULL,16) \n");
    printf("typedef long int");
  }
  else if (sizeof(long long int) == 8)
  {
    printf("#define FORTRAN_INT8_CONVERSION  \"%%lli\" \n");
    printf("#define FORTRAN_STR_TO_INT8(a) atoll(a) \n");
    printf("#define FORTRAN_STRx_TO_INT8(a) strtoll(a,NULL,16) \n");
    printf("typedef long long int");
  }
  else
    return EXIT_FAILURE;
  printf(" fortran_int8_t;\n");

  if (sizeof(short int) == 4)
  {
    printf("#define FORTRAN_INT4_CONVERSION  \"%%i\" \n");
    printf("#define FORTRAN_STR_TO_INT4(a) atoi(a) \n");
    printf("#define FORTRAN_STRx_TO_INT4(a) (short int)strtol(a,NULL,16) \n");
    printf("typedef short int");
  }
  if (sizeof(int) == 4)
  {
    printf("#define FORTRAN_INT4_CONVERSION  \"%%i\" \n");
    printf("#define FORTRAN_STR_TO_INT4(a) atoi(a) \n");
    printf("#define FORTRAN_STRx_TO_INT4(a) (int)strtol(a,NULL,16) \n");
    printf("typedef int");
  }
  else if (sizeof(long int) == 4)
  {
    printf("#define FORTRAN_INT4_CONVERSION  \"%%li\" \n");
    printf("#define FORTRAN_STR_TO_INT4(a) atol(a) \n");
    printf("#define FORTRAN_STRx_TO_INT4(a) strtol(a,NULL,16) \n");
    printf("typedef long int");
  }
  else if (sizeof(long long int) == 4)
  {
    printf("#define FORTRAN_INT4_CONVERSION  \"%%lli\" \n");
    printf("#define FORTRAN_STR_TO_INT4(a) atoll(a) \n");
    printf("#define FORTRAN_STRx_TO_INT4(a) strtoll(a,NULL,16) \n");
    printf("typedef long long int");
  }
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
