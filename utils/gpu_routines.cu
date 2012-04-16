/*

 Space for C routines of GPGPU libraries CUDA & CULA.

 Written by Miro Ilias & Daniel Kuzma in order to extend Fortran language handling
 of GPU mathematical libraries.

*/

#include "stdio.h"

// we need CUDA Toolkit routines...
#include "cublas.h"
#include "cuda.h"
#include "cuda_runtime_api.h"


/* C routine is necessary to obtain info about CUDA library (not possible from Fortran routine) */
void cudaabouts_c_(void)
{
 int cudaDriverVersion, cudaRuntimeVersion;
 cudaError_t stat;
 stat=cudaDriverGetVersion(&cudaDriverVersion); // get CUDA driver version
 stat=cudaRuntimeGetVersion(&cudaRuntimeVersion); // get runtimeversion
 printf("--- CUDA GPU info (in C routine) ---\n");
 printf(" CUDA driver version  :%d \n",cudaDriverVersion );
 printf(" CUDA runtime version :%d \n",cudaRuntimeVersion );
 return;
}

/* flush out some info about C compiler */
void c_compiler_info_(void)
{
// http://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
  printf("--- C compiler info (in C routine) ---\n");
#if __GNUC__
  printf(" GNU gcc version: %s \n",__VERSION__);
// printf("file: %s\n",__BASE_FILE__);
#endif

#if __ECC
  printf(" Intel version %s",__ECC);
#endif

#if __ICC
  printf (" Intel __ICC %i \n",__ICC);
#endif

#if __INTEL_COMPILER
  printf(" Intel compiler macro %i\n",__INTEL_COMPILER);
#endif

#if __INTEL_COMPILER_BUILD_DATE
 // printf("INTEL_COMPILER_BUILD_DATE:%s \n",__INTEL_COMPILER_BUILD_DATE);
#endif

}


/* remove empty space from routine name */
char *trim_(char *buf)
{
  char *i=(char*)buf, *j=(char*)buf;  
  do {  
    if (*i != ' ')  
      *(j++) = *i;  
  } while (*(i++));  
  return buf;
}

/*   
 Space for testing GPU routines written in C
*/
void sgemm_cuda_c_test(int *n, int *verbosity){

printf("preslo");
}


void c_gpu_math_tests_(int *n, char *routine, int *verbosity)
{
printf("--- Testing of C based routines --- \n");
printf("entering routine name=%s, size of the problem n=%i, verbosity level=%i \n", trim_(routine), *n, *verbosity);
//printf("trimmed entering routine name=%s\n",trim(routine));
// DOPLNTE SEM VASE C-ovske RUTINY, PAN KUZMA
//printf(" DOPLNTE SEM VASE C-ovske RUTINY, PAN KUZMA !!!!\n");
char *rutina;
rutina = trim_(routine);
switch (rutina){
     case "sgemm_cuda_c": sgemm_cuda_c_test(n,verbosity);
                          break;
}
/*switch(routine){
 case 'sgemm_cuda_c': sgemm_cuda_c_test(n,verbosity);
                    break;
}         */
//if (&routine=='sgemm_cuda_c') sgemm_cuda_c_test(n,verbosity);
return;

}

