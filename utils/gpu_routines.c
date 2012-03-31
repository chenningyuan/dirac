#include "stdio.h"

// we need CUDA Toolkit routines...
#include "cublas.h"
#include "cuda.h"
#include "cuda_runtime_api.h"

void c_compiler_abouts(void);

/* C routine is necessary to obtain info about CUDA library (not possible from Fortran routine) */
void cudaabouts_c_(void)
{
int cudaDriverVersion, cudaRuntimeVersion;
cudaError_t stat;

//first,  print some info about C compiler 
  c_compiler_abouts();

  stat=cudaDriverGetVersion(&cudaDriverVersion); // get CUDA driver version
  stat=cudaRuntimeGetVersion(&cudaRuntimeVersion); // get runtimeversion
  printf("--- Info about CUDA (C routine) ----\n");
  printf("CUDA driver version...%d \n",cudaDriverVersion );
  printf("CUDA runtime version...%d \n",cudaRuntimeVersion );

  return;
}

/* flush out some info about C compiler */
void c_compiler_abouts(void)
{
// http://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
  printf("C compiler info (for C routines)...\n");
#if __GNUC__
  printf("GNU gcc version: %s \n",__VERSION__);
// printf("file: %s\n",__BASE_FILE__);
#endif

#if __ECC
  printf("Intel version %s",__ECC);
#endif

#if __ICC
  printf ("Intel __ICC %i \n",__ICC);
#endif

#if __INTEL_COMPILER
  printf("Intel compiler macro %i\n",__INTEL_COMPILER);
#endif

#if __INTEL_COMPILER_BUILD_DATE
 // printf("INTEL_COMPILER_BUILD_DATE:%s \n",__INTEL_COMPILER_BUILD_DATE);
#endif

}




