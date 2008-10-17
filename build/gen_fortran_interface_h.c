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
