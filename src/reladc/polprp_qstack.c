#include<stdio.h>
#include<stdlib.h>
#include<string.h>

/*
 *     Local variables managing the heap space for the vector stack.
 *     Visible only in this source file.
 *     Values are kept during program execution.
 *     All data are on the heap.
 */

static unsigned int initialized = 0;
static unsigned int numofstacks = 0;
static int verbose = 0;
static double *** qsbase = NULL;    /* stack distributor */
static int ** qslen = NULL;  /* length of each vector on each stack */
static int * qsdep = NULL; /* depth of stack #i */



/*
 *               function qstack_init
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
int qstack_init(int *nost, int *vbflag){

  int i;
  int nos;

  if(initialized){
    printf("*** Qstack: Stack already initialized!\n");
    printf("*** Qstack: No action taken.\n");
    return 0;
  }
 
  nos = *nost;

  if(nos<1 || nos>100){
    printf("*** Qstack: Illegal value for number of stacks!\n");
    printf("*** Qstack: Should be between 1 and 100.\n");
    return -1;
  }

/* 
 * now generate requested number of stacks;
 * initialize depths, length array and main distributor
 * array
 */

  qsdep = (int *) malloc(nos*sizeof(int));
  memset(qsdep, 0, nos*sizeof(int));

  qslen = (int **)malloc(nos*sizeof(int *));
  if(qslen == NULL){
    printf("*** Qstack: Could not generate %d stacks.\n",nos);
    return -1;
  }
  for(i=0;i<nos;i++){
    qslen[i] = NULL;
  }

  qsbase = (double ***)malloc(nos*sizeof(double **));
  if(qsbase==NULL){
    printf("*** Qstack: Could not generate %d stacks.\n",nos);
    return -1;
  }
  for(i=0;i<nos;i++){
    qsbase[i] = NULL;
  }
  initialized = 1;
  numofstacks = (unsigned int)nos;
  verbose = (*vbflag)?1:0;

  printf("*** Qstack: Created %d stacks.\n",nos);
  printf("*** Qstack: Verbose mode is ");
  printf((verbose)?"on\n":"off\n");
  return nos;
}



/*
 *               function qstack_push
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
int qstack_push(int *stack, int *len, double *arr){ 

  int i,stackn;
  int curdep;
  int * newlen;
  double * curvec;
  double ** newvp;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }
  stackn = *stack - 1;  // C indexing starts at zero!

  if(*len < 1){
    printf("*** Qstack runtime problem: Illegal data length!\n");
    printf("*** No action taken.\n");
    return -1;
  }

/* first store vector item on the heap */

  curvec = (double *)malloc(*len * sizeof(double));
  memcpy(curvec,arr,*len * sizeof(double));

/* insert this address now to the corresponding stack and
 * adjust stack depths and update corresponding length
 * array.
 * create intermediate structure for copying existing pointers
 * and link it into main array at the end.
 */

  curdep = qsdep[stackn];
  curdep++;
  newvp = (double **)malloc(curdep*sizeof(double *));
  newlen = (int *)malloc(curdep*sizeof(int));

  for(i=0;i<curdep-1;i++) {
    newvp[i]  = qsbase[stackn][i];   
    newlen[i] = qslen[stackn][i];
  }
  newvp[curdep-1] = curvec;
  newlen[curdep-1] = *len;
  free(qsbase[stackn]);   /* NULL pointer save! */
  free(qslen[stackn]);   /* NULL pointer save! */
  qsbase[stackn] = newvp;
  qslen[stackn]  = newlen;

  qsdep[stackn] = curdep;
  if(verbose) printf("PUSH: Stack #%d has now %d vectors stored.\n",*stack,curdep);

  return *len;
}


/*
 *               function qstack_popf
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
int qstack_popf(int *stack, double *retvec){
  int i;
  int stackn;
  int curlen;
  int curdep;
  double ** newvp;
  int * newlen;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return 0 */

  if(qsdep[stackn] == 0) return 0;

/* fetch >>first<< vector from stack (FIFO) and copy to retvec */

  curlen = qslen[stackn][0];
  memcpy(retvec,*qsbase[stackn],curlen*sizeof(double));
/*
  for(i=0;i<curlen;i++) {
    retvec[i] = *(*qsbase[stackn] + i);
  }
*/

/* shrink stack from below and check special case for 1 vector
   on stack */

  curdep = qsdep[stackn];
  curdep--;
  free(qsbase[stackn][0]);  /* free lowest vector on selected stack */
  if(curdep == 0){    /*there was only one element, no loop required */
    free(qsbase[stackn]);
    free(qslen[stackn]);
    qsbase[stackn]  = NULL;
    qslen[stackn] = NULL;
  } else {
    newvp = (double **)malloc(curdep*sizeof(double *));
    newlen = (int *)malloc(curdep*sizeof(int));
    for(i=0;i<curdep;i++) {
      newvp[i]  = qsbase[stackn][i+1];
      newlen[i] = qslen[stackn][i+1];
    }
    free(qsbase[stackn]);
    free(qslen[stackn]);
    qsbase[stackn] = newvp;
    qslen[stackn] = newlen;
  }
  qsdep[stackn] = curdep;
  if(verbose) printf("POPF: Stack # %d has now %d vectors stored.\n",*stack,curdep);
  return curlen;
}



/*
 *               function qstack_popl
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
int qstack_popl(int *stack, double *retvec){
  int i;
  int stackn;
  int curlen;
  int curdep;
  double ** newvp;
  int * newlen;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return 0 */

  if(qsdep[stackn] == 0) return 0;

/* fetch >>last<< vector from stack (LIFO) and copy to retvec */

  curdep = qsdep[stackn];
  curlen = qslen[stackn][curdep-1];
  memcpy(retvec,qsbase[stackn][curdep-1],curlen*sizeof(double));
/*
  for(i=0;i<curlen;i++) {
    retvec[i] = *(qsbase[stackn][curdep-1] + i);
  }
*/

/* shrink stack from above and check special case for 1 vector
   on stack */

  curdep--;
  free(qsbase[stackn][curdep]);  /* free topmost vector on selected stack */
  if(curdep == 0){    /*there was only one element, no loop required */
    free(qsbase[stackn]);
    free(qslen[stackn]);
    qsbase[stackn]  = NULL;
    qslen[stackn] = NULL;
  } else {
    newvp = (double **)malloc(curdep*sizeof(double *));
    newlen = (int *)malloc(curdep*sizeof(int));
    for(i=0;i<curdep;i++) {
      newvp[i]  = qsbase[stackn][i];
      newlen[i] = qslen[stackn][i];
    }
    free(qsbase[stackn]);
    free(qslen[stackn]);
    qsbase[stackn] = newvp;
    qslen[stackn] = newlen;
  }
  qsdep[stackn] = curdep;
  if(verbose) printf("POPL: Stack # %d has now %d vectors stored.\n",*stack,curdep);
  return curlen;
}




/*
 *               function qstack_peekf
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
int qstack_peekf(int *stack, double *retvec){
  int stackn;
  int curlen;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return and do not touch retvec array */

  if(qsdep[stackn] == 0) return 0;

/* copy >>first<< vector from stack (FIFO) to retvec, do not adjust stack */

  curlen = qslen[stackn][0];
  memcpy(retvec,qsbase[stackn][0],curlen*sizeof(double));

  return curlen;
}




/*
 *               function qstack_peekl
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
int qstack_peekl(int *stack, double *retvec){
  int stackn;
  int curlen;
  int curdep;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return and do not touch retvec array */

  if(qsdep[stackn] == 0) return 0;

/* copy >>last<< vector from stack (LIFO) to retvec, do not adjust stack */

  curdep = qsdep[stackn];
  curlen = qslen[stackn][curdep-1];
  memcpy(retvec,qsbase[stackn][curdep-1],curlen*sizeof(double));

  return curlen;
}







/*
 *               function qstack_compl
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
int qstack_compl(int *stack, double *retvec){
  int stackn;
  int curlen;
  int curdep;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stack, return -1 and do not touch retvec array */

  if(qsdep[stackn] == 0) return -1;

/* compare >>last<< vector from stack (LIFO) to retvec, do not adjust stack */

  curdep = qsdep[stackn];
  curlen = qslen[stackn][curdep - 1];
  return memcmp(retvec,qsbase[stackn][curdep - 1],curlen*sizeof(double));
}






/*
 *               function qstack_meminf
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
double qstack_meminf(void){

  double usedmem=0.0;
  unsigned long int umem=0;
  int i,stackn;

  if(!initialized){
    printf("*** Qstack runtime problem: Stack not initialized!\n");
    printf("*** No action taken.\n");
    return 0.0;
  }

// contribution from qsbase

  for(stackn=0;stackn < numofstacks;stackn++){
    for(i=0;i<qsdep[stackn];i++){
      umem+=qslen[stackn][i]*sizeof(double);
    }
    umem+=qsdep[stackn]*sizeof(double *);
  }
  umem+=numofstacks*sizeof(double **);

// contribution from qslen

  for(stackn=0;stackn < numofstacks;stackn++){
    umem+=qsdep[stackn]*sizeof(int);
  }
  umem+=numofstacks*sizeof(int *);

// contribution from qsdep

  umem+=numofstacks*sizeof(int);

  usedmem=(double)(umem)/(1024.0 * 1024.0);  // mem in MB as a double
  printf("*** Qstack: total allocated memory: %lu bytes = (%.2f Megabytes).\n",umem,usedmem);

  return usedmem;
}



/*
 *               function qstack_drop
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
int qstack_drop(int * stack){

  int stackn,i,curdep;
  double ** curstack;

  if(!initialized) {
    printf("*** Qstack runtime problem: Please initialize stack\n");
    printf("*** before use! No action taken.\n");
    return -1;
  }

  if(*stack < 1 || *stack > numofstacks){
    printf("*** Qstack runtime problem: Stack number out of range!\n");
    printf("*** No action taken.\n");
    return -1;
  }

  stackn = *stack - 1;  // C indexing starts at zero!

/* if there are no elements on stackline, return 0 */

  if(qsdep[stackn] == 0) return 0;

/* clear stack line # stackn */

  curstack = qsbase[stackn];
  for(i=0;i<qsdep[stackn];i++){
    free(curstack[i]);
  }
  free(qsbase[stackn]);
  free(qslen[stackn]);
  curdep = qsdep[stackn];
  qsdep[stackn] = 0;

  qsbase[stackn] = NULL;
  qslen[stackn] = NULL;

  if(verbose) printf("*** Qstack: dropped %d elements from stack %d\n",curdep,stackn);

  return curdep;
}




/*
 *               function qstack_shutdown
 *&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 */
void qstack_shutdown(void){

  int i,stackn;
  double ** curstack;

  if(!initialized){
    printf("*** Qstack runtime problem: Stack not initialized!\n");
    printf("*** Qstack: No action taken.\n");
    return;
  }
 
  printf("*** Qstack: Shutting down stack system. Clearing heap space.\n");

  for(stackn=0;stackn < numofstacks;stackn++){
    curstack = qsbase[stackn];       //clear now stack # stackn
    for(i=0;i<qsdep[stackn];i++){
      free(curstack[i]);
    }
    free(qsbase[stackn]);
    free(qslen[stackn]);
  }
  free(qsbase);
  free(qslen);
  free(qsdep);

  qsbase = NULL;
  qslen = NULL;
  qsdep = NULL;

  numofstacks = 0;
  initialized = 0;
  verbose = 0;

  return;
}
