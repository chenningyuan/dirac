      subroutine dist_2eints_rel(NLIST,IUB,ISM,IRI,IPRTRA)
*
* Distribute relativistic integrals into
* double group format
*
* Modified good old non-relativistic routine
* dist_2eints
*
*---------------------------------------------
* Timo Fleig, December 2000
*---------------------------------------------
*
* Input :
*
* NLIST : Number of two electron integral lists
*
* IUB :
* IUB(INDEX,ILIST) = 1 => Index INDEX  in list ILIST is unbarred
*                  = 2 => Index INDEX  in list ILIST is   barred
*
* ISM :
* ISM(1,ILIST) = 1 => Symmetry between index 1 and 2 in list ILIST
*              = 0 => No Symmetry between index 1 and 2 in list ILIST
*
* ISM(2,ILIST) = 1 => Symmetry between index 3 and 4 nd  in list ILIST
*              = 0 => No Symmetry between index 3 and 4  in list ILIST
*
* ISM(3,ILIST) = 1 => Symmetry between index (12) and (34) in list ILIST
*              = 0 => No Symmetry between  index (12) and (34) in  ILIST
*
* Note : at the moment all symmetry flags are supposed to be zero !!
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "cintfo.inc"
#include "glbbas.inc"
*
      character*1 YINDEX
      dimension YINDEX(4)
      INTEGER IUB(4,*),ISM(3,*)
*
      NTESTL = 2
      NTEST = max(IPRTRA,NTESTL)
*
      if (NTEST.ge.1) then
        write(6,*)
        write(6,*) ' ================================ '
        write(6,*) '  Welcome to DIST_2EINTS_DIRAC    '
        write(6,*) ' ================================ '
        write(6,'(A,I3,A)') '    processing ',NLIST,' LISTS'
      end if
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIST2E')
*
* Some additional information:
*. Block holding largest block of integrals
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
        DO IOBSM = 1, NSMOB
          MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
        END DO
      END DO
      INTBLK = MXTSOB_L**4
      CALL MEMMAN(KLINT2,INTBLK,'ADDL  ',2,'KLINT2')
*
*
*. Distribute to double group integral lists
*  Loop over non-redundant lists
      do JLIST=1,NLIST,1
        if (NTEST.ge.2) then
          write(6,*) '+++++++++++++++++++++'
          write(6,'(A,I3)') '  This is list  ',JLIST
          write(6,*) '+++++++++++++++++++++'
        end if
*
* Barred/unbarred type of index (1 is unbarred)
        I_UB = IUB(1,JLIST)
        J_UB = IUB(2,JLIST)
        K_UB = IUB(3,JLIST)
        L_UB = IUB(4,JLIST)
*
*. Loop over double group symmetries 
*  determine corresponding point group symmetry
*  and fetch relevant integral block.
        do LFSYM=NIRR_PN+1,NSMOB,1
         do KFSYM=NIRR_PN+1,NSMOB,1
          do JFSYM=NIRR_PN+1,NSMOB,1
*
*. Double group symmetry
            LKFSYM = IDBGMULT(KFSYM,LFSYM)
            JIFSYM = IDBGMULT(1,IADJSYM(LKFSYM))
            IFSYM = IDBGMULT(JIFSYM,IADJSYM(JFSYM))
*
*. IFSYM, JFSYM, KFSYM and LFSYM are now known
            if (NTEST.ge.3) then
              write(6,*)
              write(6,*)
              write(6,'(A)') 'Double group index symmetries:'
              write(6,'(A,4I3)') ' IFSYM, JFSYM, KFSYM, LFSYM',
     &        IFSYM,JFSYM,KFSYM,LFSYM
            end if
*
* Determine boson symmetries of indices
            IBSYM = ISPSF_MO(IFSYM,I_UB)
            JBSYM = ISPSF_MO(JFSYM,J_UB)
            KBSYM = ISPSF_MO(KFSYM,K_UB)
            LBSYM = ISPSF_MO(LFSYM,L_UB)
*
            if (NTEST.ge.3) then
              write(6,'(A)') 'Corresponding point group symmetries:'
              write(6,'(A,4I3)') ' IBSYM, JBSYM, KBSYM, LBSYM',
     &        IBSYM,JBSYM,KBSYM,LBSYM
            end if
*
* Number of orbitals in this symmetry
            NITOT = NGSSHT(IBSYM)
            NJTOT = NGSSHT(JBSYM)
            NKTOT = NGSSHT(KBSYM)
            NLTOT = NGSSHT(LBSYM)
            if (NTEST.ge.3) then
              write(6,'(A)') 'Index dimensions:'
              write(6,'(A,4I3)') ' NITOT, NJTOT, NKTOT, NLTOT',
     &        NITOT,NJTOT,NKTOT,NLTOT
            end if

            NGSOBCT = 0
*
            DO LGAS = 1, NGAS
             NL = NGSSH(LBSYM,LGAS)
             DO LLGAS = 1, LGAS
               IF(LLGAS.EQ.1) THEN
                 LOFF = 1
               ELSE
                 LOFF = LOFF + NGSSH(LBSYM,LLGAS)
               END IF
             END DO
*
*
             DO KGAS = 1, NGAS
              NK = NGSSH(KBSYM,KGAS)
              DO KKGAS = 1, KGAS
                IF(KKGAS.EQ.1) THEN
                  KOFF = 1
                ELSE
                  KOFF = KOFF + NGSSH(KBSYM,KKGAS)
                END IF
              END DO
*
*
              DO JGAS = 1, NGAS
               NJ = NGSSH(JBSYM,JGAS)
               DO JJGAS = 1, JGAS
                 IF(JJGAS.EQ.1) THEN
                  JOFF = 1
                ELSE
                  JOFF = JOFF + NGSSH(JBSYM,JJGAS)
                END IF
               END DO
*
*
               DO IGAS = 1, NGAS
                NI = NGSSH(IBSYM,IGAS)
                DO IIGAS = 1, IGAS
                  IF(IIGAS.EQ.1) THEN
                    IOFF = 1
                  ELSE
                    IOFF = IOFF + NGSSH(IBSYM,IIGAS)
                  END IF
                END DO
*
*
                IGSOBCT = 0
                JGSOBCT = 0
                KGSOBCT = 0
                LGSOBCT = 0
*
                do IIGAS = 1,IGAS-1,1
                   IGSOBCT = IGSOBCT + NGSSH(IBSYM,IIGAS)
                end do
*
                do JJGAS = 1,JGAS-1,1
                   JGSOBCT = JGSOBCT + NGSSH(JBSYM,JJGAS)
                end do
                JGSOBCT = JGSOBCT * NITOT
*
                do KKGAS = 1,KGAS-1,1
                   KGSOBCT = KGSOBCT + NGSSH(KBSYM,KKGAS)
                end do
                KGSOBCT = KGSOBCT * NITOT * NJTOT
*
                do LLGAS = 1,LGAS-1,1
                   LGSOBCT = LGSOBCT + NGSSH(LBSYM,LLGAS)
                end do
                LGSOBCT = LGSOBCT * NITOT * NJTOT * NKTOT
*
                NGSOBCT = IGSOBCT + JGSOBCT + KGSOBCT + LGSOBCT
*
*
                if (NI.eq.0.or.NJ.eq.0.or.NK.eq.0.or.NL.eq.0) goto 666
*
                if (NTEST.ge.3) then
                  write(6,*)
                  write(6,'(A,4I4)') ' IGAS, JGAS, KGAS, LGAS ',
     &            IGAS, JGAS, KGAS, LGAS
                  write(6,'(A,4I5)') ' NI,NJ,NK,NL', NI,NJ,NK,NL
                  if (NTEST.ge.5) then
                    write(6,'(A,4I4)') 
     &                    'IGSOBCT,JGSOBCT,KGSOBCT,LGSOBCT ='
     &                    ,IGSOBCT,JGSOBCT,KGSOBCT,LGSOBCT
                    write(6,'(A,I5)') 'NGSOBCT =  ',NGSOBCT
                  end if
                end if
*
*. Fetch relativistic integrals from full real/imaginary lists
*  all supported environments (MOLCAS, DIRAC, ...) implemented.
                if (NTEST.ge.2) then
                   write(6,*) 'Input for GETINTN:'
                   write(6,*) 'IBSYM =  ',IBSYM
                   write(6,*) 'JBSYM =  ',JBSYM
                   write(6,*) 'KBSYM =  ',KBSYM
                   write(6,*) 'LBSYM =  ',LBSYM
                end if
*
                call getintn(WORK(KLINT2),
     &             IGAS,IBSYM,JGAS,JBSYM,KGAS,KBSYM,LGAS,LBSYM,
     &             0,0,0,0,IRI)
C                    getintn(XINT,
C                           ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &                      IXCHNG,IKSM,JLSM,ICOUL,IRI)
*
*(IJ!KL) is in WORK(KLINT2) as full (IK!JL) matrix
*
                if (NTEST.ge.2) then
                  write(6,*) 'What has GETINTN fetched?'
                  call wrtmat(WORK(KLINT2),NI*NK*NJ*NL,1,NI*NK*NJ*NL,1)
                end if
*
* Offset for double group output block
                IOFF_DBG = IOFF_DBG_LIST(JLIST,IFSYM,JFSYM,KFSYM,LFSYM)
                if (NTEST.ge.3) 
     &              write(6,'(A,I6)') 'IOFF_DBG is  ',IOFF_DBG
*
*. And now : Distribute to double group list
                do L=LOFF,LOFF+NL-1
                  do K=KOFF,KOFF+NK-1
                    do J=JOFF,JOFF+NJ-1
                      do I=IOFF,IOFF+NI-1
*
                        IREL = I - IOFF + 1
                        JREL = J - JOFF + 1
                        KREL = K - KOFF + 1
                        LREL = L - LOFF + 1
*
*. Address in input list (has form (IK!JL))
                        IJKL_IN = ((LREL-1)*NJ+JREL-1)*NI*NK
     &                           + (KREL-1)*NI
     &                           + IREL
*
*. Address in output list (has usual structure (IJ!KL))
                        IJKL_OUT = IOFF_DBG - 1
     &                    + ((LREL-1)*NKTOT+KREL-1)*NITOT*NJTOT
     &                    + (JREL-1)*NITOT+IREL
     &                    + NGSOBCT
*
*. and now process integral
                        WORK(KINT2D(IRI,JLIST)-1+IJKL_OUT) =
     &                  WORK(KLINT2-1+IJKL_IN)
                        if (NTEST.ge.5) then
                          write(6,'(A,I6)') 'IJKL_IN = ',IJKL_IN
                          write(6,'(A,1F18.12)') 
     &                          'WORK(KLINT2-1+IJKL_IN) =',
     &                           WORK(KLINT2-1+IJKL_IN)
                          write(6,'(A,I6)') 'IJKL_OUT = ',IJKL_OUT
                          write(6,'(A,1F18.12)') 
     &                          'WORK(KINT2D(1,JLIST)-1+IJKL_OUT) =',
     &                           WORK(KINT2D(1,JLIST)-1+IJKL_OUT)
                        end if
                      end do
                    end do
                  end do
                end do
*               ^ End of loop over orbital of given TS
 666            continue
               end do
              end do
             end do
            end do
*           ^ End of loop over gasspaces
          end do
         end do
        end do
*       ^ End of loop over symmetries
      end do
*     ^End of loop over output lists
*
      if (NTEST.ge.1) then
        WRITE(6,*)
        WRITE(6,*) ' Integral lists over spinors            '
        if (IRI.eq.1) write(6,*) '  (real part) '
        if (IRI.eq.2) write(6,*) '  (imaginary part) '
        WRITE(6,*) ' ---------------------------------------'
        WRITE(6,*) ' Distributed 2-electron integrals       '
        WRITE(6,*) ' ---------------------------------------'
        WRITE(6,*)
        do ILIST =1,NLIST,1
          WRITE(6,*)
          LEN_T = L2D(ILIST)
          IOFF = KINT2D(IRI,ILIST)
          do IL=1,4,1
            if (IUB(IL,ILIST).eq.1) YINDEX(IL) = 'U'
            if (IUB(IL,ILIST).eq.2) YINDEX(IL) = 'B'
          end do
          write(6,*) 'List and length :'
          write(6,'(I2,3X,A,4A1,A,3X,I6)') 
     &          ILIST,
     &          '(',(YINDEX(M),M=1,4,1),')',
     &          LEN_T
          write(6,*) ' ----------- '
          call wrtmat(WORK(IOFF),LEN_T,1,LEN_T,1)
        end do
      end if
*
      call memman(IDUM,IDUM,'FLUSM ',IDUM,'DIST2E')
*
      return
      end
