
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "linsolve.h"

void dgesv_(int *N, int *NRHS, double *A, int *LDA,
	    int *IPIV, double *B, int *LDB, int *INFO);

int linsolve(int N, int NRHS, double *A, int LDA,
	     double *rhs,  int LDrhs)
{
  int *ipiv = malloc(N*sizeof(*ipiv));
  int info;
  if (!ipiv)
    return -255;
  dgesv_(&N, &NRHS, A, &LDA, ipiv, rhs, &LDrhs, &info);
  free(ipiv);
  return info;
}

void dgelsd_(int *M, int *N, int *NRHS, 
	     double *A, int *LDA,
	     double *B, int *LDB, 
	     double *S, double *RCOND, int *rank,
	     double *work, int *lwork,
	     int *iwork,
	     int *INFO);

int linsolve_least_squares(int M, int N, int NRHS, double *A, int LDA,
			   double *rhs,  int LDrhs, double rcond)
{
  double *S, *work,dummy,wdum;
  int info,lwork,liwork,idum,*iwork,rank;
  /* Query memory requirements */
  lwork = -1;
  dgelsd_(&M,&N,&NRHS,A,&LDA,rhs,&LDrhs,&dummy,&dummy,&idum,&wdum,
	 &lwork, &liwork, &info);
  assert(info == 0 && "dgelsd workspace query failed");
  lwork = (int)wdum;
  S = malloc(sizeof(*S)*(M > N ? M : N)); /* Singular values */
  assert(S && "malloc failed");
  work = malloc(sizeof(*work)*lwork);
  assert(work && "malloc failed");
  iwork = malloc(sizeof(*iwork)*liwork);
  assert(iwork && "malloc failed");
  dgelsd_(&M,&N,&NRHS,A,&LDA,rhs,&LDrhs,
	  S,&rcond,&rank,work,
	  &lwork, &liwork, &info);
  free(S);
  free(iwork);
  return info;
}

#ifdef TEST_LINSOLVE

int main(int argc, char *argv[])
{
  int i,j,N, NRHS, res;
  double *A, *rhs;
  if (argc != 3 || sscanf(argv[1],"%i",&N) != 1 || 
      sscanf(argv[2],"%i",&NRHS) != 1)
    {
      printf("Usage: linsolve N NRHS\n");
      printf("       Read the NxN matrix A from stdin, then\n");
      printf("       read the NxNRHS right hand side matrix B\n");
      printf("       Solve AX = B and print X.\n");
      return EXIT_FAILURE;
    }
  printf("Reading %ix%i matrix A..\n",N,N);
  A = malloc(sizeof(*A)*N*N);
  assert(A && "malloc failed");
  for (i=0;i<N;i++)
    for (j=0;j<N;j++)
      {
	if (scanf("%lf",&A[i+j*N]) != 1)
	  {
	    fprintf(stderr,"Error reading matrix element (%i,%i).\n",i+1,j+1);
	    return EXIT_FAILURE;
	  }
      }
  printf("Reading %ix%i matrix B..\n",N,NRHS);
  rhs = malloc(sizeof(*rhs)*N*NRHS);
  assert(rhs && "malloc failed");
  for (i=0;i<N;i++)
    for (j=0;j<NRHS;j++)
      {
	if (scanf("%lf",&rhs[i+j*N]) != 1)
	  {
	    fprintf(stderr,"Error reading matrix element (%i,%i).\n",i+1,j+1);
	    return EXIT_FAILURE;
	  }
      }
  printf("Solving AX=B ..\n");
  res = linsolve_least_squares(4,4,3,A,4,rhs,4,-1);
  if (res == 0)
    {
      printf("X:\n");
      for (i=0;i<N;i++)
	{
	  for (j=0;j<NRHS;j++)
	    printf("%lf ",rhs[i+N*j]);
	  printf("\n");
	}
    }
  else
    {
      printf("Error, code: %i\n",res);
    }
  return res;
}

#endif


