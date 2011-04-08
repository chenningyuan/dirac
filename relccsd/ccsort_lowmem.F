!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

C  $Id$
C
#ifdef MOD_UNRELEASED
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CCSETI_LM (ICALC,EPS,NSP,IREPSP,IREPSPI)
C
      use memory_allocator
      use allocator_parameters, only : klongint, kreal
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Driver to generate the sorted integral files, using a tradeoff
C     between memory usage and number of passes for sorting
C
C     Menu driven routine. Depending on value of ICALC
C
C       1) Calculate memory requirements for the different modules
C       2) Produce desired integral files or rebuild Fock matrix
C          - OOOO, VOOO, VVOO and VOVO files
C          - VOVV file
C          - VVVV file
C          - Fock matrix
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : me
C
C---------------Common Blocks--------------------------------------
C
#include "files.inc"
#include "diis.inc"
#include "inpt.inc"
#include "symm.inc"
#include "eqns.inc"
#include "waio.h"
#include "complex.inc"
#if defined (VAR_MPI)
#include "mpif.h"
#endif
C
C---------------Calling variables--------------------------------------
C
      REAL*8 EPS(*)
      INTEGER IREPSP(NSP,4),IREPSPI(NSP,MXREP,2)
C
C---------------Local variables--------------------------------------
C
      LOGICAL ALLOC_ME,SRTALL
      INTEGER*8 NVT
Caspg, replacing memget
      integer,      allocatable, target :: a_kr(:),     a_index(:), 
     &                                     a_indk(:),   a_indl(:), 
     &                                     a_ibuf(:)
      real(kind=kreal), allocatable, target :: a_cbuf(:),   a_buf1(:), 
     &                                     a_buf2(:),   a_foo(:), 
     &                                     a_fvv(:),    a_fvo(:),
     &                                     a_ints(:),   a_intsa(:)
C
      real(kind=kreal)    :: total_mem_new = 0
      integer(kind=klongint) :: current_mem_use = 0
      integer(kind=klongint) :: start_mem_use   = 0
      integer(kind=klongint) :: mxcorr
      integer        :: freemem = 0
C
C---------------Executable code--------------------------------------
C
      call allocator_get_words_inuse(start_mem_use)
      ALLOC_ME = ICALC.GE.2
C
C  Read integrals & make fock matrices
C
      SRTALL = .TRUE.
C
C  MP2 calculations need less integrals
C  MP: DHF-ADC also needs all integral classes !
C  MP: The eqns variable hereby turns on full sorting in this case.
C
      IF (EQNS.NE.'FOCKSP') SRTALL = DOCCSD.OR.DOMP2GOLD
      IF (EQNS.EQ.'RELADC') SRTALL = .TRUE.
C
C  Allocate core for integral sorting step
C
      call alloc(a_kr,2*NSP+2, id="kr")                   
      IF (EQNS.EQ.'FOCKSP') THEN                          
         call alloc(a_INDEX,NSP*NSP*16, id="index")       
      ELSE                                                
         call alloc(a_INDEX,NSP*NSP*4, id="index")        
      ENDIF                                               
      call alloc(a_INDK,NSP*NSP    , id="indk")           
      call alloc(a_INDL,NSP*NSP    , id="indl")           
      call alloc(a_CBUF,NSP*NSP*RCW, id="cbuf")           
      IF (EQNS.EQ.'FOCKSP') THEN                          
         call alloc(a_IBUF,NSP*NSP*7*256, id="ibuf")      
      ELSE                                                
         call alloc(a_IBUF,NSP*NSP*7*16, id="ibuf" )      
      ENDIF                                               
      call alloc(a_BUF1,NSP*NSP, id="buf1" )              
      call alloc(a_BUF2,NSP*NSP, id="buf2" )              
      call alloc(a_FOO,NFOO*RCW, id="foo" )               
      call alloc(a_FVO,NFVO*RCW, id="fvo" )               
      call alloc(a_FVV,NFVV*RCW, id="fvv" )               

Caspg, replacing memget
C
C     After sorting and antisymmetrizing the OOOO,VOOO and VVOO integrals
C     parts of the above arrays are not needed anymore and could be used
C     for buffering the VOVV and VVVV integrals as well. but this is not
C     a critical shortage of memory so far.
C
C     The maximum necessary size is
C
      if (srtall) then
         NVT = MAX(NV5,NV6)
         NVT = MAX(NVT,INT(NV4,8))
         NVT = MAX(NVT,INT(JOOOO(2),8))
         NVT = MAX(NVT,INT(JVOOO(2),8))
C
C     now we get what is the maximum allocatable dimension
         call allocator_get_maxsize(freemem,kind(a_ints))
C     and rescale the result so that two buffers (ints, intsa can be allocated
          NVF = freemem / (2 * RCW) 
      else
         NVT = NV4
         NVF = NVT
      endif

C
C  Determine number of available buffers with respect to free memory
C  Records are defined for integers, calculate how many real or
C  complex words we can store on these in MINBUF
C  -----------------------------------------------------------------
C
      MINBUF = NWORDS_ON_RECORD / (RCW * IRW)
      IF (NVT.LT.NVF) THEN
         NBUF = NVT / MINBUF + 1
      ELSE
         NBUF = NVF / MINBUF 
      ENDIF
      IF (NBUF.LT.1) NBUF = 1
      NVT = NBUF * MINBUF
C
C     We must be able to store at least a NV4 size array, check that.
C     Since we use this array as work in MAKEFM we also need to store the
C     first representation of the sorted (non-triangular) OOOO and VOOO ints
C
      NVTMIN = MAX(NVT,INT(NV4,8))
      NVTMIN = MAX(NVTMIN,JOOOO(2))
      NVTMIN = MAX(NVTMIN,JVOOO(2))
      NVTMIN = NVTMIN + MINBUF 

C  we now adjust the buffer size to be a multiple of the buffer size
      if (MOD(NVTMIN,MINBUF).NE.0) then
         IEXTRA = MOD(NVTMIN,MINBUF)
         NVTMIN = NVTMIN+IEXTRA
         NVTMIN = ((NVTMIN/MINBUF)+1)*MINBUF
      endif

C
C     in this new scheme we only work with one array that will hold 
C     the different classes, one at a time

      call alloc(a_ints, NVTMIN*RCW, id="ints")
      call alloc(a_intsa,NVTMIN*RCW, id="intsa")

      call allocator_get_words_inuse(current_mem_use)
      mxcorr = current_mem_use - start_mem_use
C
      IF (ICALC.EQ.1) THEN
         WRITE (IW,*) ' '
         WRITE (IW,*) 'Note: Low(er)-memory strategy active for RELCCSD'
         WRITE (IW,*) ' '
         WRITE (IW,*) '      Now the integral sorting step will use'
         WRITE (IW,*) '      information about the amount of memory'
         WRITE (IW,*) '      available on this system and try seting'
         WRITE (IW,*) '      the amount of buffering space accordingly.'
         WRITE (IW,*) '      Due to that the number of passes for the'
         WRITE (IW,*) '      different integral classes may increase,'
         WRITE (IW,*) '      and therefore higher I/O loads are likely.'
         WRITE (IW,*) ' '
         WRITE (IW,1001) "reading+sorting integrals",mxcorr
         WRITE (IW,1001) "integral sorting buffers only",(2*NVTMIN*RCW)
         WRITE (IW,*) ' '
         RETURN
      ENDIF
C
      CALL RDINTS_LM(IPRNT,NSP,a_INDEX,a_KR,a_INDK,a_INDL,a_CBUF,a_IBUF,
     &               NVTMIN,a_BUF1,a_BUF2,IREPSP,IREPSPI,EPS,
     &               a_FOO,a_FVO,a_FVV,a_ints,a_intsa,BREIT,SRTALL,
     &               INTERFACE)

      call dealloc(a_kr)
      call dealloc(a_INDEX)
      call dealloc(a_INDK)
      call dealloc(a_INDL)
      call dealloc(a_CBUF)
      call dealloc(a_IBUF) 
      call dealloc(a_BUF1)
      call dealloc(a_BUF2)
      call dealloc(a_FOO)
      call dealloc(a_FVO)
      call dealloc(a_FVV)
      call dealloc(a_ints)
      call dealloc(a_intsa)
C
      RETURN
1001  FORMAT (" Memory for ",A," :",T50,I10," 8-byte words")
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RDINTS_LM(IPRNT_LOC,NSP,INDEX,KR,INDK,INDL,CBUF,IBUF,
     &                     NVT,BUF1,BUF2,IREPSP,IREPSPI,
     &                     EPS,FOO,FVO,FVV,INTS,INTSA,BREIT_LOC,SRTALL,
     &                     INTERFACE_LOC)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Reads 1-electron and 2-electron integrals
C     Constructs Fock matrix elements
C     This version reads from files created by TMOONE & ROTRAN or DIRAC
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : me
C
C---------------Common Blocks--------------------------------------
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "complex.inc"
#include "results.inc"
#include "files.inc"
#include "ccpar.inc"
#include "waio.h"
#include "inpt.inc"
#include "dcbgen.h"
#include "dcbham.h"
#if defined (VAR_MPI)
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
C
C---------------Calling variables--------------------------------------
C
      INTEGER NVT
      INTEGER NSP
      INTEGER KR(-NSP:NSP),INDK(*),INDL(*)
      INTEGER INDEX(NSP,NSP,4,*),IBUF(7,*)
      INTEGER IREPSP(NSP,4),IREPSPI(NSP,MXREP,2)
      REAL*8 FOO(*),FVO(*),FVV(*)
      REAL*8 BUF1(*),BUF2(*),CBUF(RCW,*)
      REAL*8 EPS(*)
      REAL*8 INTS(RCW,*),INTSA(RCW,*)
      LOGICAL BREIT_LOC,SRTALL
      CHARACTER*10 INTERFACE_LOC
C
C---------------Local variables--------------------------------------
C
      REAL*8 E0,E1,ORIG_TOTSCF
      CHARACTER*10 DATEX,TIMEX*8
      LOGICAL FOCKSP,TOBE
      INTEGER*8 NF,NL,IND,ICOUNT,ITRSGN
      CHARACTER*10 FNODE,LFNAM
C
C
C---------------Executable code--------------------------------------
C
C
      IF (SRTALL) THEN
         WRITE (IW,1004) 'OOOO',NV1,'VOOO',NV2,'VVOO',NV3,'VOVO',NV4,
     &                   'VOVV',NV5,'VVVV',NV6
      ELSE
         WRITE (IW,1004) 'OOOO',NV1,'VOOO',NV2,'VVOO',NV3,'VOVO',NV4
      ENDIF
C
      CALL FLSHFO (IW)
      FOCKSP = EQNS.EQ.'FOCKSP'
C
      IF (FOCKSP) THEN
         CALL FKINDEX (NSP,IREPSP,IREPSPI,INDEX)
      ELSE
         CALL MKINDEX (NSP,IREPSP,IREPSPI,INDEX)
      ENDIF
C
C  Check buffer size for passing of VOVV and VVVV integrals
C
      IF (MOD(NVT*IRW*RCW,NWORDS_ON_RECORD).NE.0) then
          write (*,*) 'invalid buffer size in rdints'
          write (*,*) 'NVT',NVT
          write (*,*) 'IRW',IRW
          write (*,*) 'RCW',RCW
          write (*,*) 'NWORDS_ON_RECORD',NWORDS_ON_RECORD
          write (*,*) 'MOD(NVT*IRW*RCW,NWORDS_ON_RECORD)',
     &                 MOD(NVT*IRW*RCW,NWORDS_ON_RECORD)
          CALL QUIT("INVALID BUFFER SIZE IN RDINTS")
      endif
C
C  INITIALIZE READING OF COULOMB INTEGRALS
C
      IF (MYPROC.EQ.0) THEN
         FNODE = 'MDCINT    '
      ELSE
         FNODE = LFNAM('MDCINT')
      ENDIF
      OPEN(MDINT,FILE=FNODE,FORM='UNFORMATTED')
      REWIND(MDINT)
      IF (INTERFACE_LOC(1:7).EQ.'MOLFDIR') THEN
         IF (MYPROC.EQ.MASTER) THEN
C           The oldest MOLFDIR file format
            READ (MDINT,ERR=10000,END=10) DATEX,TIMEX,NKR,
     &       (IDUM,I=1,4*NKR),(KR(I),KR(-I),I=1,NKR)
            GO TO 11
C           The newer MOLFDIR file format = identical to older DIRAC
   10       READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,
     &       (KR(I),KR(-I),I=1,NKR)
   11       CONTINUE
         ENDIF
      ELSEIF (INTERFACE_LOC(1:5).EQ.'DIRAC') THEN
         IF (MYPROC.EQ.MASTER.OR.INTERFACE_LOC(6:6).EQ.'6') THEN
         READ (MDINT,ERR=10000,END=10000) DATEX,TIMEX,NKR,
     &                                    (KR(I),KR(-I),I=1,NKR)
         ENDIF
         MREADER = MASTER
      ELSE
         WRITE (IW,'(//1X,3A)') 'Interface ',INTERFACE_LOC,
     & ' not defined'
         CALL QUIT('Adapt RDINTS')
      ENDIF
c     write (*,*) 'debug: fnode ',fnode,'  mreader',mreader
C
C     In the older schemes, the slaves have no knowledge of nkr and kr, broadcast
C     and check for consistency afterwards
C
#if defined (VAR_MPI)
      IF (NMPROC .GT. 1) THEN
         CALL MPI_BCAST(NKR,1,my_MPI_INTEGER,
     &                  MASTER,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(KR,2*NSP+1,my_MPI_INTEGER,
     &                  MASTER,MPI_COMM_WORLD,IERR)
      END IF
#endif
      IF (2*NKR.NE.NSP) 
     &   CALL QUIT("INCONSISTENT MRCONEE AND MDCINT FILES")
      ITRSGN = 1
C
C  Initialize writing to sorted integral files
C
      CALL WAIO_OPEN(ITAPT+0)
C
C  We may have already sorted these integrals when restarting, in that case jump ahead
C
      IF (IMSTAT(1).GT.1.AND..NOT.FOCKSP) THEN
         WRITE (IW,1010)
         GOTO 150
      ENDIF
C
C     RE-ENTRY POINT IN CASE THE GAUNT INTERACTION HAS TO BE ADDED
C
    1 ICOUNT = 0
C
C     GET THE OOOO, VOOO, VVOO AND VOVO INTEGRALS, one at a time...
C
C     first the oooo, symmetric and antisymmetric parts
      call integral_pass(interface_loc,'OOOO',ints,intsa,nvt,nv1,1,1,
     &                   icount,kr,nsp,index,itrsgn,indk,indl,
     &                             cbuf,ibuf,mreader,fnode)
C
C     then the vooo, symmetric and antisymmetric parts
    2 ICOUNT = 0
      call integral_pass(interface_loc,'VOOO',ints,intsa,nvt,nv2,2,2,
     &                   icount,kr,nsp,index,itrsgn,indk,indl,
     &                             cbuf,ibuf,mreader,fnode)
C
C     then the vvoo, symmetric and antisymmetric parts
    3 ICOUNT = 0
      call integral_pass(interface_loc,'VVOO',ints,intsa,nvt,nv3,4,3,
     &                   icount,kr,nsp,index,itrsgn,indk,indl,
     &                             cbuf,ibuf,mreader,fnode)
C
C     then the vovo, symmetric and antisymmetric parts
    4 ICOUNT = 0
      call integral_pass(interface_loc,'VOVO',ints,intsa,nvt,nv4,6,4,
     &                   icount,kr,nsp,index,itrsgn,indk,indl,
     &                             cbuf,ibuf,mreader,fnode)
C
C     This is a restart point : save the files and update status.
C
  150 CONTINUE
      IMSTAT(1) = 3
      IF(EQNS.NE.'RELADC') CALL CCDUMP
C
C  We may have already sorted these integrals when restarting, in that case jump ahead
C
      IF (IMSTAT(2).GT.1.AND..NOT.FOCKSP) THEN
         WRITE (IW,1011)
         GOTO 550
      ENDIF
C
C     Determine # of passes necessary for treating VOVV integrals
C     on a specific node
C
      IF (.NOT.SRTALL) GOTO 600
C
C     then the vovv, symmetric and antisymmetric parts
    5 ICOUNT = 0
      call integral_pass(interface_loc,'VOVV',ints,intsa,nvt,nv5par,14,
     &                   5,icount,kr,nsp,index,itrsgn,indk,indl,
     &                             cbuf,ibuf,mreader,fnode)
C
C     then the vovv, symmetric and antisymmetric parts
    6 ICOUNT = 0
      call integral_pass(interface_loc,'VVVV',ints,intsa,nvt,nv6par,16,
     &                   6,icount,kr,nsp,index,itrsgn,indk,indl,
     &                             cbuf,ibuf,mreader,fnode)

C
  500 CLOSE(MDINT)
C
C  NOW PASS THE BREIT INTEGRAL FILE IF DESIRED
C  Breit is only used with MOLFDIR in this version
C
      IF (BREIT_LOC.AND.ITRSGN.EQ.1) THEN
         WRITE (IW,1000) 'Coulomb',DATEX,TIMEX,ICOUNT
         OPEN(MDINT,FILE='MDBINT',FORM='UNFORMATTED')
         REWIND(MDINT)
         IF (MYPROC.EQ.MASTER) THEN
            READ (MDINT,ERR=10001,END=10001) DATEX,TIMEX,NKR2
            IF (NKR2.NE.NKR) GOTO 10003
         ENDIF
         ITRSGN = -1
         GOTO 1
      ELSEIF (BREIT_LOC) THEN
         WRITE (IW,1000) 'Breit',DATEX,TIMEX,ICOUNT
      ELSE
         WRITE (IW,1000) 'Coulomb',DATEX,TIMEX,ICOUNT
      ENDIF
C
C     This is a restart point : save the files and update status.
C
  550 CONTINUE
      IMSTAT(2) = 3
      IF(EQNS.NE.'RELADC') CALL CCDUMP
C
  600 WRITE (IW,1006)
C     since we don't keep the integrals, we set up the fock matrix here
      CALL GETFM (INDEX,NSP,EPS,BUF1,BUF2,ints,intsa,FOO,FVO,FVV,E0,E1)

C
C     This is a restart point : save the files and update status.
C
      CALL WAIO_SAVE(ITAPT+0)
      IMSTAT(3) = 3
      IF(EQNS.NE.'RELADC') CALL CCDUMP
C
      ESCF = E0 + E1 
C
C     SK - a first step towards a consistent treatment of 
C     NORECMP == .TRUE. + USEOE == .TRUE. 
C     when we run with the MOLMF(3) hamiltonian
      IF(NORECMP .and. USEOE) THEN
        WRITE(IW,'(/A,/A,/A,/A)') 
     & ' ***********************************************************',
     & ' * Info from RDINTS: No recomputation of SCF energy.       *',
     & ' * Value is taken from the file DFCOEF instead.            *',
     & ' ***********************************************************'
        ORIG_TOTSCF = 0.0D0
        IF( MYPROC .eq. MASTER )THEN
          INQUIRE(FILE='DFCOEF',EXIST=TOBE)
          IF(.NOT.TOBE) THEN
              WRITE(IW,'(A)')
     &          'RDINTS: Coefficient file DFCOEF not found !'
              CALL QUIT('RDINTS: Coefficient file DFCOEF not found !')
          ENDIF
          CALL OPNFIL(LUCOEF,'DFCOEF','OLD','RDINTS')
          CALL REACMO_no_work(LUCOEF,'DFCOEF',DUMMY,DUMMY,IDUMMY,
     &                        ORIG_TOTSCF,1)
          CLOSE(LUCOEF,STATUS='KEEP')
        ENDIF
#if defined VAR_MPI
        CALL MPI_BCAST(ORIG_TOTSCF,1,MPI_REAL8,MASTER,
     &                 MPI_COMM_WORLD,IERR)
#endif
        ESCF = ORIG_TOTSCF
        WRITE(IW,1013) ESCF
      ELSE 
        WRITE(IW,1003) ECORE,E0,E1,ESCF,ESCF+ECORE
        ESCF = ESCF + ECORE
      END IF
      ETOT = ESCF
C
 1000 FORMAT (/' Reading ',A7,' integrals :'/' File date :',T18,A10
     &/' File time :',T20,A8,/' # of integrals',I16)
 1003 FORMAT (//' Nuclear repulsion + core energy :',T40,F25.15
     &/' Zero order electronic energy :',T40,F25.15 
     &/' First order electronic energy :',T40,F25.15 
     &/' Electronic energy :',T40,F25.15 
     &/' SCF energy :',T40,F25.15) 
 1004 FORMAT (/' Expanding and sorting integrals to unique types :'/
     & (' Type ',A4,' :',I10,' integrals'))
 1005 FORMAT (/' Need ',I6,' passes to sort ',A4,' integrals')
 1006 FORMAT (/' Finished sorting of integrals')
 1010 FORMAT (/' Skipping integral sorting step 1')
 1011 FORMAT (/' Skipping integral sorting step 2')
 1012 FORMAT (/' Skipping Fock-matrix build')
 1013 FORMAT (//' SCF energy :',T40,F25.15)
      RETURN
10000 CALL QUIT('ERROR READING HEADER MDCINT')
10001 CALL QUIT('ERROR READING HEADER MDBINT')
10003 CALL QUIT('MDCINT AND MDBINT FILES ARE INCOMPATIBLE')
10010 CALL QUIT('ERROR READING INTEGRALS FROM MDCINT')
10011 CALL QUIT('ERROR READING INTEGRALS FROM MDBINT')
      END

C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
      subroutine integral_pass(interface,integral_class,ints,intsa,nvt,
     &                         nvpar,curr_typ,io_unit,
     &                         icount,kr,nsp,index,itrsgn,indk,indl,
     &                                   cbuf,ibuf,mreader,fnode)
      IMPLICIT INTEGER (A-Z)
C
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "complex.inc"
#include "results.inc"
#include "files.inc"
#include "ccpar.inc"
#include "waio.h"
#if defined (VAR_MPI)
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif

      INTEGER nonzr,nze
      INTEGER NSP, ICOUNT
      INTEGER KR(-NSP:NSP),INDK(*),INDL(*)
      INTEGER INDEX(NSP,NSP,4,*),IBUF(7,*)
      CHARACTER*10 FNODE

      LOGICAL BREIT,SRTALL
      CHARACTER*10 INTERFACE
      REAL*8 CBUF(RCW,*)

      character(len=4) :: integral_class ! label for this integral class
      integer:: curr_typ ! integral type we are currently processing
      integer:: io_unit, nvt, nvpar, itrsgn
      real(kind=8) :: ints(RCW,*), intsa(RCW,*)

C     first open file for this integral class
      CALL WAIO_OPEN(ITAPT+io_unit)
C
C     Determine # of passes necessary for treating a given integral class 
C     on a specific node
C
      NPASS = (NVPAR-1) / NVT + 1
      NSECT = 1
      NL    = 0
c     write (*,*) 'npass,nvt,nvpar',npass,nvt,nvpar
      WRITE (IW,1005) NPASS,integral_class
      CALL FLSHFO (IW)
C
C   start of IPASS loop for this integral class
C
C ----------------------------
      DO IPASS = 1, NPASS
C ----------------------------

      NF = NL + 1
      NL = MIN(NL+NVT,NVPAR)
      N  = NL - NF + 1
      NSEC = N*IRW*RCW / NWORDS_ON_RECORD
      IF (MOD(N*IRW*RCW,NWORDS_ON_RECORD).NE.0) NSEC = NSEC + 1
C
      IF (MYPROC.EQ.MASTER.OR.INTERFACE(6:6).EQ.'6') THEN
         CLOSE (MDINT,STATUS='KEEP')
         IF (BREIT.AND.ITRSGN.EQ.-1) THEN
             OPEN(MDINT,FILE='MDBINT',FORM='UNFORMATTED')
         ELSE
             OPEN(MDINT,FILE=FNODE,FORM='UNFORMATTED')
         ENDIF
         READ (MDINT)
         MREADER = MASTER
      ENDIF
C
      CALL XCOPY (nvt,A0,0,ints ,1)
      CALL XCOPY (nvt,A0,0,intsa,1)
C
c     write (*,*) 'curr typ',curr_typ
c     write (*,*) 'nze before calling get2ints',nze
  200 continue 
      if (curr_typ.lt.10) then
         if (npass.gt.1) 
     &   call quit('multipass sort not implemented for this class')
       
         CALL GET2INTS (MDINT,KR,NSP,INDEX,ITRSGN,1,INDK,INDL,
     &                  NONZR,NZE,CBUF,IBUF,MREADER)
      else
         CALL GET2INTS (MDINT,KR,NSP,INDEX,ITRSGN,curr_typ,INDK,INDL,
     &                  NONZR,NZE,CBUF,IBUF,MREADER)
      endif
c     write (*,*) 'nze after calling get2ints',nze
      
      IF (NONZR.LT.0) THEN
         IF (MREADER.LT.NMPROC-1.AND.INTERFACE(6:6).EQ.'6') THEN
             MREADER = MREADER + 1
             GOTO 200
         ELSE
             GOTO 210
         ENDIF
      ELSE
         ICOUNT = ICOUNT + NONZR
      ENDIF
      DO 202 INZE = 1, NZE
            ITYP = IBUF(1,INZE)
            IREP = IBUF(2,INZE) 
            IND1 = IBUF(3,INZE)
            IND2 = IBUF(4,INZE)
            ITRS = IBUF(5,INZE)
            JTRS = IBUF(6,INZE)
            I    = IBUF(7,INZE)

            if (ITYP.EQ.1.and.curr_typ.eq.1) THEN
               IND = IOOOOTT(IREP) + (IND1-1)*NOOT(IREP) + IND2
               ints(1,IND)  = CBUF(1,I) * ITRS
               IF (CARITH) ints(2,IND) = CBUF(2,I) * JTRS
            elseif (ITYP.EQ.-1.and.curr_typ.eq.1) THEN
               IND = IOOOOTT(IREP) + (IND1-1)*NOOT(IREP) + IND2
               intsa(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) intsa(2,IND) = CBUF(2,I) * JTRS
            elseif (ITYP.EQ.2.and.curr_typ.eq.2) THEN
               IND = IVOOOT(IREP)  + (IND1-1)*NVO(IREP)  + IND2
               ints(1,IND)  = CBUF(1,I) * ITRS
               IF (CARITH) ints(2,IND) = CBUF(2,I) * JTRS
            elseif (ITYP.EQ.-2.and.curr_typ.eq.2) THEN
               IND = IVOOOT(IREP)  + (IND1-1)*NVO(IREP)  + IND2
               intsa(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) intsa(2,IND) = CBUF(2,I) * JTRS
            elseif (ITYP.EQ.4.and.curr_typ.eq.4) THEN
               IND = IVVOOTT(IREP) + (IND1-1)*NVVT(IREP) + IND2
               ints(1,IND)  = CBUF(1,I) * ITRS
               IF (CARITH) ints(2,IND) = CBUF(2,I) * JTRS
            elseif (ITYP.EQ.-4.and.curr_typ.eq.4) THEN
               IND = IVVOOTT(IREP) + (IND1-1)*NVVT(IREP) + IND2
               intsa(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) intsa(2,IND) = CBUF(2,I) * JTRS
            elseif (ITYP.EQ.6.and.curr_typ.eq.6) THEN
               IND = IVOVO(IREP)   + (IND1-1)*NVO(IREP)  + IND2
               ints(1,IND)  = CBUF(1,I) * ITRS
               IF (CARITH) ints(2,IND) = CBUF(2,I) * JTRS
            elseif (ITYP.EQ.-6.and.curr_typ.eq.6) THEN
               IND = IVOVO(IREP)   + (IND1-1)*NVO(IREP)  + IND2
               intsa(1,IND) = CBUF(1,I) * ITRS
               IF (CARITH) intsa(2,IND) = CBUF(2,I) * JTRS

            elseif (iabs(ityp).eq.14) then
               IND = IVOVVT(IREP)  + (IND1-1)*NVO(IREP)  + IND2
            elseif (iabs(ityp).eq.16) then
               IND = IVVVVTT(IREP) + (IND1-1)*NVVT(IREP) + IND2
            endif

C for classes 5 and 6 we need some extra work in determining the index..

            if (iabs(ityp).eq.14) then

                IF (ITYP.EQ.curr_typ) THEN
                   IF (IND.GT.IDIST(3,2,IREP).AND.
     &                 IND.LE.IDIST(4,2,IREP)) THEN
C                now check relative offset if we need it in
C                the current pass
                      IND = IND -IDIST(3,2,IREP) + IDIST(5,2,IREP)
                      IF (IND.GE.NF.AND.IND.LE.NL) THEN
                         IND = IND - NF + 1
c                 write (*,*) 'ind + ',ind
                         ints(1,IND)  = CBUF(1,I) * ITRS
                         IF (CARITH) ints(2,IND) = CBUF(2,I) * JTRS
                      ENDIF
                   ENDIF
                elseif (ITYP.EQ.-curr_typ) THEN
                   IF (IND.GT.IDIST(3,2,IREP).AND.
     &                 IND.LE.IDIST(4,2,IREP)) THEN
C             Check absolute address and decide if this integral goes
C             to this node
                      IND = IND -IDIST(3,2,IREP) + IDIST(5,2,IREP)
                      IF (IND.GE.NF.AND.IND.LE.NL) THEN
                         IND = IND - NF + 1
c                 write (*,*) 'ind - ',ind
                         intsa(1,IND) = CBUF(1,I) * ITRS
                         IF (CARITH) intsa(2,IND) = CBUF(2,I) * JTRS
                      ENDIF
                   ENDIF 
                endif

             elseif (iabs(ityp).eq.16) then

                IF (ITYP.EQ.curr_typ) THEN
                   IF (IND.GT.IDIST(3,1,IREP).AND.
     &                 IND.LE.IDIST(4,1,IREP)) THEN
C                now check relative offset if we need it in
C                the current pass
                      IND = IND -IDIST(3,1,IREP) + IDIST(5,1,IREP)
                      IF (IND.GE.NF.AND.IND.LE.NL) THEN
                         IND = IND - NF + 1
c                 write (*,*) 'ind + ',ind
                         ints(1,IND)  = CBUF(1,I) * ITRS
                         IF (CARITH) ints(2,IND) = CBUF(2,I) * JTRS
                      ENDIF
                   ENDIF
                elseif (ITYP.EQ.-curr_typ) THEN
                   IF (IND.GT.IDIST(3,1,IREP).AND.
     &                 IND.LE.IDIST(4,1,IREP)) THEN
C             Check absolute address and decide if this integral goes
C             to this node
                      IND = IND -IDIST(3,1,IREP) + IDIST(5,1,IREP)
                      IF (IND.GE.NF.AND.IND.LE.NL) THEN
                         IND = IND - NF + 1
c                 write (*,*) 'ind - ',ind
                         intsa(1,IND) = CBUF(1,I) * ITRS
                         IF (CARITH) intsa(2,IND) = CBUF(2,I) * JTRS
                      ENDIF
                   ENDIF 
                endif

             endif

  202 CONTINUE
      GOTO 200
C
C     ADD THE DIRECT AND EXCHANGE INTEGRALS
C
  210 continue
      if (curr_typ.lt.10) then
         CALL XAXPY (nvpar,-A1,intsa,1,ints,1)
      else
         CALL XAXPY (N,-A1,intsa,1,ints,1)
      endif
C
C     ADD BREIT AND COULOMB INTEGRALS
C
      IF (BREIT.AND.ITRSGN.LT.0) THEN
c        write (*,*) 'before rread in pass',ipass
         if (curr_typ.lt.10) then
            CALL RREAD(ITAPT+io_unit,intsa,nvpar*RCW*IRW,1)
            CALL XAXPY(nvpar,A1,intsa,1,ints,1)
         else
            CALL RREAD(ITAPT+io_unit,intsa,N*RCW*IRW,NSECT)
            CALL XAXPY(N,A1,intsa,1,ints,1)
         endif
      ENDIF
c     write (*,*) 'before rwrit in pass',ipass
      if (curr_typ.lt.10) then
         CALL RWRIT(ITAPT+io_unit,ints,nvpar*RCW*IRW,1)
      else
         CALL RWRIT(ITAPT+io_unit,ints,N*RCW*IRW,NSECT)
         NSECT = NSECT + NSEC
      endif
C
C  end of IPASS LOOP
C ----------------------
      ENDDO
C ----------------------
c     write (*,*) 'before waio_save',io_unit
      CALL WAIO_SAVE(ITAPT+io_unit)
      WRITE(IW,*) integral_class,' sorting done.'
c     write (*,*) 'debug, icount',icount
      CALL FLSHFO (IW)
C
 1005 FORMAT (/' Need ',I6,' passes to sort ',A4,' integrals')
      RETURN
10000 CALL QUIT('ERROR READING HEADER MDCINT')
10001 CALL QUIT('ERROR READING HEADER MDBINT')
10003 CALL QUIT('MDCINT AND MDBINT FILES ARE INCOMPATIBLE')
10010 CALL QUIT('ERROR READING INTEGRALS FROM MDCINT')
10011 CALL QUIT('ERROR READING INTEGRALS FROM MDBINT')
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#else
      subroutine ccseti_lm
      end
#endif
