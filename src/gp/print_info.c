#include <stdio.h>
#include "git_info.h"
#include "build_info.h"

/*  Print the DIRAC configuration/build info generated at build time by the executable.
 *
 * Printing of long string is available in C, not in Fortran, therefore whole C routine must be employed.
 *
 * */
void print_info(void) 
{
    printf("\n ---  Git version information  ---\n");
    printf("                Git branch : %s\n",GIT_BRANCH);
    printf("      Last git commit hash : %s\n",GIT_COMMIT_HASH);
    printf("    Last git commit author : %s\n",GIT_COMMIT_AUTHOR);
    printf("      Last git commit date : %s\n",GIT_COMMIT_DATE);

    printf("\n---  Configuration and build info  ---\n");
    printf("              Who compiled : %s\n",USER_NAME);
    printf("        Compiled on server : %s\n",HOST_NAME);
    printf("          Operating system : %s\n",SYSTEM);
    printf("            Python version : %s\n",PYTHON_VERSION);
    printf("             CMake version : %s\n",CMAKE_VERSION);
    printf("           CMake generator : %s\n",CMAKE_GENERATOR);
    printf("          CMake build type : %s\n",CMAKE_BUILD_TYPE);
    printf("       MPI parallelization : %s\n",ENABLE_MPI);
    printf("           64-bit integers : %s\n",ENABLE_64BIT_INTEGERS);
    printf("          Fortran compiler : %s\n",FORTRAN_COMPILER);
    printf("  Fortran compiler version : %s\n",FORTRAN_COMPILER_VERSION);
    printf("    Fortran compiler flags : %s\n",FORTRAN_COMPILER_FLAGS);
    printf("                C compiler : %s\n",C_COMPILER);
    printf("        C compiler version : %s\n",C_COMPILER_VERSION);
    printf("          C compiler flags : %s\n",C_COMPILER_FLAGS);
    printf("              C++ compiler : %s\n",CXX_COMPILER);
    printf("      C++ compiler version : %s\n",CXX_COMPILER_VERSION);
    printf("        C++ compiler flags : %s\n",CXX_COMPILER_FLAGS);
    printf("            Static linking : %s\n",STATIC_LINKING);
    printf("      Builtin BLAS library : %s\n",ENABLE_BUILTIN_BLAS);
    printf("    Builtin LAPACK library : %s\n",ENABLE_BUILTIN_LAPACK);
    printf("    Mathematical libraries : %s\n",MATH_LIBS);
    printf("        Explicit libraries : %s\n",EXPLICIT_LIBS);
    printf("       Modules definitions : %s\n",DEFINITIONS);
    printf("       Configuration time :  %s\n",CONFIGURATION_TIME);
    printf("-----------------------------------------------------------------\n");
}

