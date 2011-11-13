C ** COMMON block for Mulliken population analysis
C     POPLAB - labels for population analysis
C     ILABDF = 1: base label definitions on primitive AO-labels
C     ILABDF = 2: base label definitions on primitive SO-labels
C     LABDEF    : user defined labels
C     ADDSML    : collect small components for a given center into one 
C                 label (AO-def)
C                 collect small components for given center combination
C                 into one label (SO-def)
C     NPOPLAB   : number of labels for population analysis
C     IPOPLAB(MAXLAB) - pointer from AO/SO-label to POPLAB
C     NSPI  - number fo spinors to analyze
C     ISPI  - offset for spinors to analyze
C     IPRPOP - print level
C
#include <dcblab.h>
      CHARACTER POPLAB*12,VECPOP*72
      COMMON/CBCPOP/POPLAB(MAXLAB),VECPOP(2)
      LOGICAL ADDSML,DONETP,LABDEF
      COMMON/CBLPOP/ADDSML,DONETP,LABDEF
      COMMON/CBIPOP/IPRPOP,IPOPLAB(MAXLAB),NPOPLAB,ILABDF
