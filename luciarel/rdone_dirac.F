*  Import one-electron integrals from DIRAC environment
*
      subroutine rdone_dirac(ONEINT,LUONE,DIRACONE,MORB,NDIM,
     &                       ISPINFREE,IBOSYM,IDUMAR,DUMAR)
*
      implicit real*8 (A-H,O-Z)
*
#include "input.inc"
*
      character*14 CHA(8)
      common/CECORE/ECORE
*
      character*7 DIRACONE
      logical BREIT
*
      dimension IDUMAR(*),IBOSYM(*),NSTR(2)
      dimension ONEINT((2*NDIM)**2*2),DUMAR(*)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
      open(LUONE,FILE=DIRACONE,FORM='UNFORMATTED')
      read(LUONE) MORB,BREIT,ECORE_DIR,NFSYM,NZ
      read(LUONE) NSYMRP,(CHA(I),I=1,NSYMRP),(IDUMAR(I),I=1,NSYMRP),
     &                                         (NSTR(I),I=1,NFSYM)
      read(LUONE)
      read(LUONE)
      read(LUONE) (IDUMAR(ID),IDUMAR(ID),DUMAR(ID),ID=1,MORB,1),
     &            (IBOSYM(I),I=1,MORB,1)
      if (NTEST.ge.3) then
        write(6,*) 
        write(6,*) 'Testing read MRCONEE'
        write(6,*) 'MORB = ',MORB
        write(6,*) 'BREIT = ',BREIT
        write(6,*) 'ECORE_DIR = ',ECORE_DIR
        write(6,*) 'Number of spinors per fermion symmetry:'
        write(6,'(2I5)') (NSTR(I),I=1,NFSYM)
        write(6,*) 'Read boson symmetries of spinors:'
        do IB=1,MORB,1
          write(6,'(A11,I5,A3,I12)') 'Spinor no. ',IB,' : ',IBOSYM(IB)
        end do
      end if
*
      if (2*NDIM.ne.MORB) then
        write(6,*) 'Orbital dimensions of '
        write(6,*) 'DIRAC and input do not match.'
        write(6,*) 'DIRAC : ',MORB,' spinors'
        write(6,*) 'Input : ',NDIM,' orbitals'
        call abend2('Quitting.')
      end if
*
* Read the integrals
      NDIMSQ = MORB*MORB
      if (ISPINFREE.eq.1) then
        read(LUONE) (ONEINT(I),TDUM,I=1,MORB*MORB)
      else if (ISPINFREE.eq.0) then
        read(LUONE) (ONEINT(I),ONEINT(I+NDIMSQ),I=1,MORB*MORB)
*
* Set IBOSYM array so that U->alpha and B->beta spinors
*              To get correct Fermion symmetry labeling in
*              LUCIAREL we need to specify same(!) gerade boson
*              irrep for all gerade orbitals and same ungerade
*              boson irrep for all ungerade orbitals (note: which
*              gerade or ungerade boson irrep we choose doesn't
*              matter).
        IOFFFS = 0
        IRROFF = 0
        do IFSYM = 1,NFSYM,1
          do ISP = 1+IOFFFS,NSTR(IFSYM)+IOFFFS,1
            IBOSYM(ISP) = 0 + IRROFF
          end do
*
          if (NTEST.ge.3) then
            write(6,*) ' Reset IBOSYM array, ferm symm. :',IFSYM
            write(6,*) (IBOSYM(I),I=1+IOFFFS,NSTR(IFSYM)+IOFFFS,1)
          end if
*
          IOFFFS = NSTR(IFSYM)
          IRROFF = 2
          if (NZ.eq.1) IRROFF = 1
        end do
*
      end if
*
      if (NTEST.ge.5) then
        write(6,*)
        MORB_L = MORB
        write(6,*) 'Real One-electron ints from dirac:'
        call wrtmat(ONEINT(1),MORB_L,MORB_L,MORB_L,MORB_L)
        if (ISPINFREE.eq.0) then
          write(6,*) 'Imaginary One-electron ints from dirac:'
          call wrtmat(ONEINT(1+NDIMSQ),MORB_L,MORB_L,MORB_L,MORB_L)
        end if
      end if  
*
      ECORE = ECORE_DIR
* 
      return
      end
