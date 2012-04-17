/*

 Space for C routines of GPGPU libraries CUDA & CULA.

 Written by Miro Ilias & Daniel Kuzma in order to extend Fortran language handling
 of GPU mathematical libraries.

*/

#include "stdio.h"
#include <sys/time.h>
#include <string.h>

// we need CUDA Toolkit routines...
#include "cublas.h"
#include "cuda.h"
#include "cuda_runtime_api.h"

/* cula */
#include <cula.h>
#include <culablas.h>
#include <culapack.h>

/* cblas */
#if defined HAVE_MKL_BLAS
/* see /opt/intel/mkl/include/mkl_cblas.h */
#include "mkl_cblas.h"
#pragma message "Using Intel MKL <mkl_cblas.h> interface"
#else
/* see /usr/include/cblas.h */
#include "cblas.h"
#pragma message "Using GNU <cblas.h> interface"
#endif

/* lapack */
#if defined HAVE_MKL_LAPACK
/* see /opt/intel/mkl/include/mkl_lapack(e).h  */
//#include <mkl_lapack.h>
#include <mkl_lapacke.h>
#pragma message "Using Intel MKL <mkl_lapacke.h> interface"
#else
#include <clapack.h>
#pragma message "Using GNU <clapack.h> interface"
#endif


int NRHS = 1; //for cula sgesv


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

double getTime()
{
 struct timeval tx;
 double t;
 gettimeofday(&tx, NULL);
 t = (tx.tv_sec) * 1000.0;
 t += (tx.tv_usec) / 1000.0;
 return t;
}


void sgemm_cuda_c_test(int *n, int *verbosity){
        
        cublasStatus status;
        float* h_A;
        float* h_B;
        float* h_C;
        float* h_C_ref;
        float* d_A = 0;
        float* d_B = 0;
        float* d_C = 0;
        float alpha = 1.0f;
        float beta = 0.0f;
        int i;
        float error_norm;
        float ref_norm;
        float diff;
        double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
        double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu
        int n2,n3;
        
        n3 = *n;
        n2 = n3*n3;
        //printf("\n\n\n N2:%i N:%i N*:%i\n\n\n",n2,n3,*n);
        /* Initialize CUBLAS */
        printf("\nCUBLAS SGEMM test running...\n\n");

        status = cublasInit();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! CUBLAS initialization error\n");
                exit(EXIT_FAILURE);
        } else { fprintf(stdout,"GPU cublas initialized.\n"); }


        /* Allocate host memory for the matrices */
        h_A = (float*)malloc(n2 * sizeof(h_A[0]));
        if (h_A == 0) {
                fprintf (stderr, "!!!! host memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated \n"); }


        h_B = (float*)malloc(n2 * sizeof(h_B[0]));
        if (h_B == 0) {
                fprintf (stderr, "!!!! host memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"B mtx allocated \n"); }

        h_C = (float*)malloc(n2 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C mtx allocated \n"); }

        /* Fill the matrices with test data */
        for (i = 0; i < n2; i++) {
                h_A[i] = rand() / (float)RAND_MAX;
                h_B[i] = rand() / (float)RAND_MAX;
                h_C[i] = 2;//rand() / (float)RAND_MAX;
        }
        if(*verbosity>0)fprintf(stdout,"A,B,C matrixes filled \n");

        /* Allocate device memory for the matrices A,B and C */
        status = cublasAlloc(n2, sizeof(d_A[0]), (void**)&d_A);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated on GPU \n"); }

        status = cublasAlloc(n2, sizeof(d_B[0]), (void**)&d_B);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B mtx allocated on GPU \n"); }

        status = cublasAlloc(n2, sizeof(d_C[0]), (void**)&d_C);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C mtx allocated on GPU \n"); }
        /* Initialize the device matrices with the host matrices */
        gpu_time_1 = getTime();
        status = cublasSetVector(n2, sizeof(h_A[0]), h_A, 1, d_A, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1;
        if(*verbosity==2)printf("Filling time of A mtx on GPU: %f sekund\n",gpu_time/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx filled on GPU \n"); }
        gpu_time_1 = getTime();
        status = cublasSetVector(n2, sizeof(h_B[0]), h_B, 1, d_B, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity==2)printf("Filling time of B mtx on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B mtx filled on GPU \n"); }
        gpu_time_1 = getTime();
        status = cublasSetVector(n2, sizeof(h_C[0]), h_C, 1, d_C, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity==2)printf("Filling time of C mtx on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C mtx filled on GPU \n"); }
        
        printf("CPU cblas_sgemm running....");
        fflush(stdout);
    
        cpu_time_1 = getTime();

/*
 void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc);
*/
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n3, n3,n3, alpha, h_A, n3, h_B, n3, beta, h_C, n3);

        cpu_time_2 = getTime();

        printf("...done\n");
        
        h_C_ref = h_C;

        /* Clear last error */
        cublasGetError();

        /* Performs operation using cublas */
       // matrix-matrix multiplication  na GPU ...
        printf("GPU cublasSgemm running....");
        fflush(stdout);
        gpu_time_1 = getTime();
        cublasSgemm('n', 'n', n3, n3, n3, alpha, d_A, n3, d_B, n3, beta, d_C, n3);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        
        status = cublasGetError();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! kernel execution error.\n");
                exit(EXIT_FAILURE);
        } else  {fprintf(stdout,"...done.\n"); }
        if(*verbosity==2)printf("Time of computation on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        /* Allocate host memory for reading back the result from the device memory */
        h_C = (float*)malloc(n2 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"New C mtx ready for GPU data. \n");}

        /* Read the result back */
        gpu_time_1 = getTime();
        status = cublasGetVector(n2, sizeof(h_C[0]), d_C, 1, h_C, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity==2)printf("Obtaining time from GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (read C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C mtx obtained from GPU device. \n");}
        //printf("...done\n");
        printf("\n\n\n");
        if(*verbosity>0)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",h_C_ref[0],h_C[0]);
        printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",h_C_ref[n3-1],h_C[n3-1]);
        
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n\n\n",gpu_time/1000);
        
        /* Memory clean up */
        free(h_A);
        free(h_B);
        free(h_C);
        //free(h_C_ref);
        status = cublasFree(d_A);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (A)\n");
                exit(EXIT_FAILURE);
        }
        status = cublasFree(d_B);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (B)\n");
                exit(EXIT_FAILURE);
        }
        status = cublasFree(d_C);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (C)\n");
                exit(EXIT_FAILURE);
        }

        /* Shutdown */
        status = cublasShutdown();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! shutdown error (A)\n");
                exit(EXIT_FAILURE);
        }

}

void sgemv_cuda_c_test(int *n, int *verbosity){

cublasStatus status;
        float* h_A;
        float* h_B;
        float* h_C;
        float* h_C_ref;
        float* d_A = 0;
        float* d_B = 0;
        float* d_C = 0;
        float alpha = 1.0f;
        float beta = 0.0f;
        int i,j;
        float error_norm;
        float ref_norm;
        float diff;
        double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
        double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu

        int n2,n3;
        
        n3 = *n;
        n2 = n3*n3;
        //printf("\n\n\n N2:%i N:%i N*:%i\n\n\n",n2,n3,*n);
        /* Initialize CUBLAS */
        printf("\nCUBLAS SGEMV test running...\n\n");

        status = cublasInit();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! CUBLAS initialization error\n");
                exit(EXIT_FAILURE);
        } else { fprintf(stdout,"GPU cublas initialized.\n"); }


        /* Allocate host memory for the matrices */
        h_A = (float*)malloc(n2 * sizeof(h_A[0]));
        if (h_A == 0) {
                fprintf (stderr, "!!!! host memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated \n"); }


        h_B = (float*)malloc(n3 * sizeof(h_B[0]));
        if (h_B == 0) {
                fprintf (stderr, "!!!! host memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"B vector allocated \n"); }

        h_C = (float*)malloc(n3 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C vector allocated \n"); }

        /* Fill the matrices with test data */
        for (i = 0; i < n2; i++) {
                h_A[i] =  rand() / (float)RAND_MAX;
                
        }
        for (j = 0; j < n3; j++) {
                
                h_B[j] = rand() / (float)RAND_MAX;
                h_C[j] = 2;//rand() / (float)RAND_MAX;
        }    
        if(*verbosity>0)fprintf(stdout,"A matrixes filled and B,C vector filled\n");

        /* Allocate device memory for the matrices A and vector B and C */
        status = cublasAlloc(n2, sizeof(d_A[0]), (void**)&d_A);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated on GPU \n"); }

        status = cublasAlloc(n3, sizeof(d_B[0]), (void**)&d_B);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B vector allocated on GPU \n"); }

        status = cublasAlloc(n3, sizeof(d_C[0]), (void**)&d_C);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C vector allocated on GPU \n"); }
        /* Initialize the device matrices with the host matrices */
        gpu_time_1 = getTime();
        status = cublasSetVector(n2, sizeof(h_A[0]), h_A, 1, d_A, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1;
        if(*verbosity==2)printf("Filling time of A mtx on GPU: %f sekund\n",gpu_time/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx filled on GPU \n"); }
        gpu_time_1 = getTime();
        status = cublasSetVector(n3, sizeof(h_B[0]), h_B, 1, d_B, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity==2)printf("Filling time of B vector on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B vector filled on GPU \n"); }
        gpu_time_1 = getTime();
        status = cublasSetVector(n3, sizeof(h_C[0]), h_C, 1, d_C, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity==2)printf("Filling time of C vector on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (write C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C vector filled on GPU \n"); }
        
        printf("CPU cblas_sgemv running....");
        fflush(stdout);
    
        cpu_time_1 = getTime();

/*
void cblas_sgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY);
*/
        cblas_sgemv(CblasRowMajor, CblasNoTrans, n3,n3, alpha, h_A, n3, h_B, 1, beta, h_C, 1);

        cpu_time_2 = getTime();

        printf("...done\n");
        
        h_C_ref = h_C;

        /* Clear last error */
        cublasGetError();

        /* Performs operation using cublas */
       // matrix-matrix multiplication  na GPU ...
        printf("GPU cublasSgemv running....");
        fflush(stdout);
        gpu_time_1 = getTime();
        cublasSgemv('t', n3, n3, alpha, d_A, n3, d_B, 1, beta, d_C, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        
        status = cublasGetError();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! kernel execution error.\n");
                exit(EXIT_FAILURE);
        } else  {fprintf(stdout,"...done.\n"); }
        if(*verbosity==2)printf("Time of computation on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        /* Allocate host memory for reading back the result from the device memory */
        h_C = (float*)malloc(n2 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"New C mtx ready for GPU data. \n");}

        /* Read the result back */
        gpu_time_1 = getTime();
        status = cublasGetVector(n3, sizeof(h_C[0]), d_C, 1, h_C, 1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1 + gpu_time;
        if(*verbosity==2)printf("Obtaining time from GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! device access error (read C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C mtx obtained from GPU device. \n");}
        //printf("...done\n");
        printf("\n\n\n");
        if(*verbosity>0)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",h_C_ref[0],h_C[0]);
        printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",h_C_ref[n3-1],h_C[n3-1]);
        
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n\n\n",gpu_time/1000);
        
        /* Memory clean up */
        free(h_A);
        free(h_B);
        free(h_C);
        //free(h_C_ref);
        status = cublasFree(d_A);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (A)\n");
                exit(EXIT_FAILURE);
        }
        status = cublasFree(d_B);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (B)\n");
                exit(EXIT_FAILURE);
        }
        status = cublasFree(d_C);
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! memory free error (C)\n");
                exit(EXIT_FAILURE);
        }

        /* Shutdown */
        status = cublasShutdown();
        if (status != CUBLAS_STATUS_SUCCESS) {
                fprintf (stderr, "!!!! shutdown error (A)\n");
                exit(EXIT_FAILURE);
        }
}


void sgemm_cula_c_test(int *n, int *verbosity){
        
        culaStatus status;
        culaFloat* d_A = 0;
        culaFloat* d_B = 0;
        culaFloat* d_C = 0;
        culaInt* IPIV = NULL;
        float* h_C_ref;
        float* h_A = 0;
        float* h_B = 0;
        float* h_C = 0;
        float alpha = 1.0f;
        float beta = 0.0f;
        int i;
        float error_norm;
        float ref_norm;
        float diff;
        double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
        double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu

        int n2,n3;
        
        n3 = *n;
        n2 = n3*n3;
        //printf("\n\n\n N2:%i N:%i N*:%i\n\n\n",n2,n3,*n);
        /* Initialize CUBLAS */
        printf("\nCULA SGEMM test running...\n\n");

        status = culaInitialize();
        if (status != culaNoError) {
                fprintf (stderr, "!!!! CULA initialization error\n");
                exit(EXIT_FAILURE);
        } else { fprintf(stdout,"GPU cula initialized.\n"); }


        /* Allocate host memory for the matrices */
        h_A = (float*)malloc(n2 * sizeof(h_A[0]));
        if (h_A == 0) {
                fprintf (stderr, "!!!! host memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated \n"); }


        h_B = (float*)malloc(n2 * sizeof(h_B[0]));
        if (h_B == 0) {
                fprintf (stderr, "!!!! host memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"B mtx allocated \n"); }

        h_C = (float*)malloc(n2 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C mtx allocated \n"); }

        /* Allocate device memory for the matrices A,B and C */
        d_A = (culaFloat*)malloc(n2*sizeof(culaFloat));
        if (d_A == 0) {
                fprintf (stderr, "!!!! device memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated on graphic card\n"); }
        
        d_B = (culaFloat*)malloc(n2*sizeof(culaFloat));
        if (d_B == 0) {
                fprintf (stderr, "!!!! device memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B mtx allocated on graphic card\n"); }
        d_C = (culaFloat*)malloc(n2*sizeof(culaFloat));
        if (d_C == 0) {
                fprintf (stderr, "!!!! device memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C mtx allocated on graphic card\n"); }


        /* Fill the matrices with test data */
        for (i = 0; i < n2; i++) {
               d_A[i] = h_A[i] = rand() / (float)RAND_MAX;
                //d_A[i] = rand() / (float)RAND_MAX;
               d_B[i] = h_B[i] = rand() / (float)RAND_MAX;
                //d_B[i] = rand() / (float)RAND_MAX;
                h_C[i] = 2;//rand() / (float)RAND_MAX;
                d_C[i] = 21;//rand() / (float)RAND_MAX;
        }
        if(*verbosity>0)fprintf(stdout,"A,B,C matrixes filled \n");

        
       
        printf("CPU cblas_sgemm running....");
        fflush(stdout);
    
        cpu_time_1 = getTime();
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n3, n3,n3, alpha, h_A, n3, h_B, n3, beta, h_C, n3);
        cpu_time_2 = getTime();

        printf("...done\n");
        
        h_C_ref = h_C;

        /* Clear last error */
        cublasGetError();

        /* Performs operation using cublas */
       // matrix-matrix multiplication  na GPU ...
        printf("GPU culaSgemm running....");
        fflush(stdout);
        gpu_time_1 = getTime();
        status = culaSgemm('N', 'N', n3, n3, n3, alpha, d_A, n3, d_B, n3, beta, d_C, n3);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1;
        //checkStatus(status);
        //status = cublasGetError();
        if (status != culaNoError) {
                fprintf (stderr, "!!!! kernel execution error.\n");
                exit(EXIT_FAILURE);
        } else  {fprintf(stdout,"...done.\n"); }
        if(*verbosity==2)printf("Time of computation on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        /* Allocate host memory for reading back the result from the device memory */
        
        printf("\n\n\n");
        if(*verbosity>0)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",h_C_ref[0],d_C[0]);
        printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",h_C_ref[n2-1],d_C[n2-1]);
        
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n\n\n",gpu_time/1000);
        
        /* Memory clean up */
        free(d_A);
        free(d_B);
        free(d_C);
        free(h_A);
        free(h_B);
        free(h_C);
        //free(h_C_ref);

        /* Shutdown */
        culaShutdown();


}

void sgemv_cula_c_test(int *n, int *verbosity){
        
        culaStatus status;
        culaFloat* d_A = 0;
        culaFloat* d_B = 0;
        culaFloat* d_C = 0;
        culaInt* IPIV = NULL;
        float* h_C_ref;
        float* h_A = 0;
        float* h_B = 0;
        float* h_C = 0;
        float alpha = 1.0f;
        float beta = 0.0f;
        int i;
        float error_norm;
        float ref_norm;
        float diff;
        double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
        double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu

        int n2,n3;
        
        n3 = *n;
        n2 = n3*n3;
        //printf("\n\n\n N2:%i N:%i N*:%i\n\n\n",n2,n3,*n);
        /* Initialize CUBLAS */
        printf("\nCULA SGEMV test running...\n\n");

        status = culaInitialize();
        if (status != culaNoError) {
                fprintf (stderr, "!!!! CULA initialization error\n");
                exit(EXIT_FAILURE);
        } else { fprintf(stdout,"GPU cula initialized.\n"); }


        /* Allocate host memory for the matrices */
        h_A = (float*)malloc(n2 * sizeof(h_A[0]));
        if (h_A == 0) {
                fprintf (stderr, "!!!! host memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated \n"); }


        h_B = (float*)malloc(n3 * sizeof(h_B[0]));
        if (h_B == 0) {
                fprintf (stderr, "!!!! host memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"B mtx allocated \n"); }

        h_C = (float*)malloc(n3 * sizeof(h_C[0]));
        if (h_C == 0) {
                fprintf (stderr, "!!!! host memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        }  else {if(*verbosity>0)fprintf(stdout,"C mtx allocated \n"); }

        /* Allocate device memory for the matrices A,B and C */
        d_A = (culaFloat*)malloc(n2*sizeof(culaFloat));
        if (d_A == 0) {
                fprintf (stderr, "!!!! device memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated on graphic card\n"); }
        
        d_B = (culaFloat*)malloc(n3*sizeof(culaFloat));
        if (d_B == 0) {
                fprintf (stderr, "!!!! device memory allocation error (B)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B mtx allocated on graphic card\n"); }
        d_C = (culaFloat*)malloc(n3*sizeof(culaFloat));
        if (d_C == 0) {
                fprintf (stderr, "!!!! device memory allocation error (C)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"C mtx allocated on graphic card\n"); }


        /* Fill the matrices with test data */
        for (i = 0; i < n2; i++) {
               d_A[i] = h_A[i] = rand() / (float)RAND_MAX;
        }
        for (i = 0; i < n3; i++) {
               d_B[i] = h_B[i] = rand() / (float)RAND_MAX;
                h_C[i] = 2;//rand() / (float)RAND_MAX;
                d_C[i] = 21;//rand() / (float)RAND_MAX;
        }
        if(*verbosity>0)fprintf(stdout,"A,B,C matrixes filled \n");
       
        printf("CPU cblas_sgemv running....");
        fflush(stdout);
    
        cpu_time_1 = getTime();
        cblas_sgemv(CblasRowMajor, CblasNoTrans, n3,n3, alpha, h_A, n3, h_B, 1, beta, h_C, 1);
        cpu_time_2 = getTime();

        printf("...done\n");
        

        /* Clear last error */
        cublasGetError();

        /* Performs operation using cublas */
       // matrix-matrix multiplication  na GPU ...
        printf("GPU culaSgemm running....");
        fflush(stdout);
        gpu_time_1 = getTime();
        status = culaSgemv('t',n3,n3,alpha,d_A,n3,d_B,1,beta,d_C,1);
        gpu_time_2 = getTime();
        gpu_time = gpu_time_2 - gpu_time_1;
        //checkStatus(status);
        //status = cublasGetError();
        if (status != culaNoError) {
                fprintf (stderr, "!!!! kernel execution error.\n");
                exit(EXIT_FAILURE);
        } else  {fprintf(stdout,"...done.\n"); }
        if(*verbosity==2)printf("Time of computation on GPU: %f sekund\n",(gpu_time_2 - gpu_time_1)/1000);
        /* Allocate host memory for reading back the result from the device memory */
        
        printf("\n\n\n");
        if(*verbosity>0)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",h_C[0],d_C[0]);
        printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",h_C[n3-1],d_C[n3-1]);
        
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n\n\n",gpu_time/1000);
        
        /* Memory clean up */
        free(d_A);
        free(d_B);
        free(d_C);
        free(h_A);
        free(h_B);
        free(h_C);
        //free(h_C_ref);  

        /* Shutdown */
        culaShutdown();


}

void sgesv_cula_c_test(int *n, int *verbosity){
    int i, info;
    int n3 = *n;

    double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
    double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu

    culaStatus status;
    
    culaFloat* A = NULL;
    culaFloat* B = NULL;
    culaFloat* X = NULL;
    culaInt* IPIV = NULL;

    culaFloat one = 1.0f;
    culaFloat thresh = 1e-6f;
    culaFloat diff;

    void checkStatus(culaStatus status)
{
    char buf[80];
    if(!status) return;
    culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
    printf("%s\n", buf);
    culaShutdown();
    exit(EXIT_FAILURE);
}

    void verify_result(void)
    {
    printf("Verifying Result\n");
    for(i = 0; i < n3; ++i)
    {
        diff = X[i] - B[i];
        if(diff < 0.0f)
            diff = -diff;
        if(diff > thresh)
            printf("Result check failed:  i=%d  X[i]=%f  B[i]=%f", i, X[i], B[i]);
    }
    } 

    void set_variables(void)
    {
    // Set A to the identity matrix
    memset(A, 0, n3*n3*sizeof(culaFloat)); // first zero all numbers
    for(i = 0; i < n3; ++i)
        A[i*n3+i] = one;
    // Set B to a random matrix (see note at top)
    for(i = 0; i < n3; ++i)
        B[i] = (culaFloat)rand();
    // copy all stuff into X, which will get roots
    memcpy(X, B, n3*sizeof(culaFloat));
    memset(IPIV, 0, n3*sizeof(culaInt));
    }

    printf("-------------------\n");
    printf("   GPU CULA_SGESV  \n");
    printf("-------------------\n");

    printf("Allocating Matrix A, vectors B, X and IPIV, size N=%d \n",n3);
    A = (culaFloat*)malloc(n3*n3*sizeof(culaFloat));
    B = (culaFloat*)malloc(n3*sizeof(culaFloat));
    X = (culaFloat*)malloc(n3*sizeof(culaFloat));
    IPIV = (culaInt*)malloc(n3*sizeof(culaInt));
    if(!A || !B || !IPIV)
        exit(EXIT_FAILURE);

    printf("Initializing CULA\n");
    status = culaInitialize();
    checkStatus(status);

    set_variables();

    printf("Calling culaSgesv...");
    gpu_time_1 = getTime();
    status = culaSgesv(n3, NRHS, A, n3, IPIV, X, n3);
    gpu_time_2 = getTime();
    checkStatus(status);
    printf("done\n");

    verify_result();

    printf("Shutting down CULA\n\n");
    culaShutdown();

    printf("----------------------\n");
    printf("   CPU CLAPACK_SGESV  \n");
    printf("----------------------\n");

    set_variables();

    printf("Starting clapack_sgesv ....");

    cpu_time_1 = getTime();
#if defined HAVE_MKL_LAPACK
/*
lapack_int LAPACKE_sgesv( int matrix_order, lapack_int n, lapack_int nrhs,
                          float* a, lapack_int lda, lapack_int* ipiv, float* b,
                          lapack_int ldb );
*/

    info = LAPACKE_sgesv(CblasRowMajor, n3, NRHS, A,  n3, IPIV, B, n3);
#else // we have GNU clapack.h
/*
GNU version:
 int clapack_sgesv(const enum CBLAS_ORDER Order, const int N, const int NRHS,
                   float *A, const int lda, int *ipiv, float *B, const int ldb)
*/

    info = clapack_sgesv(CblasRowMajor, n3, NRHS, A,  n3, IPIV, B, n3);
#endif
    cpu_time_2 = getTime();
    if (info != 0) { fprintf(stderr, "clapack_sgesv failure with error %d\n", info);}
    else { fprintf(stdout,"clapack_sgesv successfully done.\n"); }

    verify_result();
    printf("\n\n\n");
        if(*verbosity>0)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",B[0],X[0]);
        printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",B[n3-1],X[n3-1]);
        
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n\n\n",(gpu_time_2-gpu_time_1)/1000);

    free(A);
    free(B);
    free(X);
    free(IPIV);
}

void sgegrf_cula_c_test(int *n, int *verbosity){
        double gpu_time_1,gpu_time_2,gpu_time; // zaciatocny a konecny cas na gpu
        double cpu_time_1,cpu_time_2; // zaciatocny a konecny cas na cpu
        
        int i;
        
         culaFloat thresh = 1e-6f;
         culaFloat diff;
        
         srand ( time(NULL) );
        
         void checkStatus(culaStatus status)
{
    char buf[80];
    if(!status) return;
    culaGetErrorInfoString(status, culaGetErrorInfo(), buf, sizeof(buf));
    printf("%s\n", buf);
    culaShutdown();
    exit(EXIT_FAILURE);
}
        
        culaStatus status;
        int n2;
        int info = 0;
        int n3 = *n;
        culaFloat* A = NULL;
        culaFloat* TAU_r = NULL;
        float* B = 0;
        float* TAU_r1 = 0;
        float* WORK_r = 0;
        
        n2 = n3*n3;
        
        
        void verify_result(void)
    {
    printf("Verifying Result\n");
    for(i = 0; i < n2; ++i)
    {
        diff = B[i] - A[i];
        if(diff < 0.0f)
            diff = -diff;
        if(diff > thresh)
        
           if(*verbosity>0) printf("Result check failed:  i=%d  B[i]=%f  A[i]=%f\n", i, B[i], A[i]);
    }
    for(i = 0; i < n3; ++i)
    {
        diff = TAU_r1[i] - TAU_r[i];
        if(diff < 0.0f)
            diff = -diff;
        if(diff > thresh)
        
           if(*verbosity>0) printf("Result check failed:  i=%d  TAU_r1[i]=%f  TAU_r[i]=%f\n", i, TAU_r1[i], TAU_r[i]);
    }
    } 
        
        B = (float*)malloc(n2 * sizeof(B[0]));
        if (B == 0) {
                fprintf (stderr, "!!!! host memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"B mtx allocated \n"); }
        
        TAU_r1 = (float*)malloc(n3 * sizeof(TAU_r1[0]));
        if (TAU_r1 == 0) {
                fprintf (stderr, "!!!! host memory allocation error (TAU_r1)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"TAU_r1 allocated \n"); }
        
        WORK_r = (float*)malloc(n3 * sizeof(WORK_r[0]));
        if (WORK_r == 0) {
                fprintf (stderr, "!!!! host memory allocation error (WORK_r)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"WORK_r allocated \n"); }
        
        A = (culaFloat*)malloc(n2*sizeof(culaFloat));
        TAU_r = (culaFloat*)malloc(n3*sizeof(culaFloat));
        
        if (A == NULL) {fprintf (stderr, "!!!! device memory allocation error (A)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"A mtx allocated on graphic card\n"); }
        
        if (TAU_r == NULL) {fprintf (stderr, "!!!! device memory allocation error (TAU_r)\n");
                exit(EXIT_FAILURE);
        } else {if(*verbosity>0)fprintf(stdout,"TAU_r allocated on graphic card\n"); }
        
        for (i=0;i<n2;i++)
        {
          A[i] = rand() / (float)RAND_MAX;
          B[i] = A[i];
        }
        printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",B[n2-1],A[n2-1]);
        
        printf("-------------------\n");
        printf("   GPU CULA_SGEQRF  \n");
        printf("-------------------\n");
        
        printf("Initializing CULA\n");
        status = culaInitialize();
        checkStatus(status);
        gpu_time_1 = getTime();
        status = culaSgeqrf(n3, n3, A, n3, TAU_r);
        gpu_time_2 = getTime();
        checkStatus(status);
        
        printf("Shutting down CULA\n\n");
        culaShutdown();
        
        printf("----------------------\n");
        printf("   CPU LAPACKE_sgeqrt (MKL)  \n");
        printf("----------------------\n");
      
        printf("Starting clapack_sgeqrf ....");
        cpu_time_1 = getTime();
#if defined HAVE_MKL_LAPACK
        info=LAPACKE_sgeqrf(LAPACK_COL_MAJOR,n3, n3, B, n3, TAU_r1);
#else
        printf("clapack_sgeqrf not present in GNU clapack.h ! returning\n");
        return;
#endif
        cpu_time_2 = getTime();
        printf("\nComputing ... done\n");
        
        verify_result();
        printf(" ... done\n");
        printf("\n\n\n");
        printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",B[n2-1],A[n2-1]);
        if(*verbosity>0)printf("Results : CPU - First element of matrix - %f / GPU - First element of matrix - %f\n",TAU_r1[0],TAU_r[0]);
        printf("Results : CPU - Last element of matrix - %f / GPU - Last element of matrix - %f\n\n\n",TAU_r1[n3-1],TAU_r[n3-1]);
        printf("\n\n\n");
        printf("Time of computation on CPU: %f sekund\n",(cpu_time_2-cpu_time_1)/1000);
        printf("Overall time of computation and transfer on GPU: %f sekund\n\n\n",(gpu_time_2-gpu_time_1)/1000);
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
/*switch (rutina){
     case "sgemm_cuda_c": sgemm_cuda_c_test(n,verbosity);
                          break;
} */
/*switch(routine){
 case 'sgemm_cuda_c': sgemm_cuda_c_test(n,verbosity);
                    break;
}         */
if (strcmp(rutina,"sgemm_cuda_c")==0) sgemm_cuda_c_test(n, verbosity);
if (strcmp(rutina,"sgemv_cuda_c")==0) sgemv_cuda_c_test(n, verbosity);
if (strcmp(rutina,"dgemm_cuda_c")==0) sgemm_cuda_c_test(n, verbosity);
if (strcmp(rutina,"sgemm_cula_c")==0) sgemm_cula_c_test(n, verbosity);
if (strcmp(rutina,"sgemv_cula_c")==0) sgemv_cula_c_test(n, verbosity);
if (strcmp(rutina,"sgesv_cula_c")==0) sgesv_cula_c_test(n, verbosity);
if (strcmp(rutina,"sgeqrf_cula_c")==0) sgegrf_cula_c_test(n, verbosity);

//if (strcmp(rutina,"dgemm_cuda_c")==0) printf("\n\n\ndgemm_cuda_c\n\n\n");

//printf("\n\n\n%s\n\n\n",rutina);

return;

}

