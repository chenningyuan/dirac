! $Id$
!
! subroutine  : GetFDEPotential_FromFile
!
! authors     : Andre Gomes, Christoph Jacob (VUA)
!
! description : this subroutine will read an embedding potential from 
!               file and construct a matrix representation of that 
!               potential in the ao basis currently in use, via numerical
!               integration. 
!
!               the routine expects that the grid file has the following
!               format:
!
!               line        one : number of grid point
!               lines to to n-1 : quadruplet of numbers, where the
!                                 first three are the (x,y,z) coords
!                                 and the last is the value of the
!                                 embedding potential at that point
!                                 times the weight of the quadrature
!                                 [= v(i)*w(i)]
!               line          n : a non-negative integer, like  -42 =)
!
!               the default name for the grid file is EMBPOT, but that
!               can be overriden via a menu option. see routine fdeinp_
!               for further details
!
! arguments   :
!
!   1. fdemat : matrix to hold the embedding potential. can be the
!               one-electron fock matrix, or any other matrix one
!               wishes to use... 
!   2.  lwork : 
!   3.   work : workspace array 
!
! in the end  : the embedding potential is added to fdemat. 
!
      SUBROUTINE GetFDEPotential_FromFile(fdemat,work,lwork)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "dftcom.h"
#include "ccom.h"
#include "inforb.h"
#include "fde.h"
      integer p,q
      dimension work(lwork)
#if defined (PRG_DIRAC)
      dimension fdemat(nbast,nbast)
#else
!     maybe better to make implicit size for dalton, or in general?
!     dimension fdemat(*)
      dimension fdemat(nnbasx)
#endif
      logical debug, embpot_found
!
      call qenter('GetFDEPotential_FromFile')
#include "memint.h"

      ndimension   = nbast 
      debug        = .false.

      if (debug) write (LUPRI,*) 'debug: before opening FDE pot file'

!     before proceeding, we check if the FDE  potential is available, 
      write (LUPRI,*) ' '
      write (LUPRI,*) ' FDE Potential will be read from file: ',
     &                 FDEPOT(1:20)
      write (LUPRI,*) ' '
      inquire(FILE=FDEPOT,EXIST=embpot_found)
      if (.not.embpot_found) then
         call quit ('Fatal error: FDE potential file not present')
      endif

      if (debug) then
!        for testing purposes, we will use this file to dump the grid we read.
         open(lutesting,
     &        FORM   = 'FORMATTED',
     &        STATUS = 'UNKNOWN',
     &        ACCESS = 'SEQUENTIAL',
     &        FILE   =  'TSTPOT')
      endif

      open(lufdepot,
     &       FORM   = 'FORMATTED',
     &       STATUS = 'UNKNOWN',
     &       ACCESS = 'SEQUENTIAL',
     &       FILE   =  FDEPOT)

!     we start processing the grid file by first getting the number of point in 
!     the grid. the size of the grid is then used to allocate workspace for 
!     reading the actual grid, along with the embedding potential at each point
      rewind lufdepot
      read (lufdepot,*) npoints
      if (debug) write (LUPRI,*) 'debug: #points to be read:',npoints

      if (npoints.gt.0) then
         CALL MEMGET('REAL',KX,NPOINTS,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KY,NPOINTS,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KZ,NPOINTS,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KW,NPOINTS,WORK,KFREE,LFREE)
      else
         call quit ('Fatal error: empty FDE potential file')
      endif

      CALL RYO_GRID(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),
     &              lufdepot,NPOINTs)
      close (lufdepot)

      if (debug) then
         call header ('1-electron matrix, before fde update',-1)
         call printdmat(1,ndimension,fdemat)
      endif

!     we allocate the variables to hold the value of the basis functions
!     at the grid ponts, and also some work arrays
      call memget ('REAL',kgao ,  ndimension,work,kfree,lfree)
      call memget ('REAL',kbuf ,4*ndimension,work,kfree,lfree)
      call memget ('REAL',kncnt,4*ndimension,work,kfree,lfree)

      nder   = 0
      iprint = 0
      irepop = 0

      if (debug) write (lutesting,'(I10)') npoints

      do i = 1, npoints              ! begin numerical integration
         ioff = i-1

!       we get the value of all gaussians at the current grid point with 
!       routines from dft
#if defined (PRG_DIRAC)
         call getsos(work(kgao),dummy,dummy,dummy,
     &               work(kcnt),
     &               work(kx+ioff),
     &               work(ky+ioff),
     &               work(kz+ioff),
     &               work(kbuf),
     &               work,lwork,ndimension,0,.false.,iprfde)
#else
         call dftaos(work(kgao),dummy,dummy,dummy,dummy,
     &               work(kcnt),
     &               work(kx+ioff),
     &               work(ky+ioff),
     &               work(kz+ioff),
     &               work,lwork,ndimension,
     &               .false.,.false.,DFTHRI,iprdfe)
#endif

      if (debug)                                 ! dump of grid points 
     &   write (lutesting,'(4(E25.18,2X))') work(kx+ioff),work(ky+ioff),
     &                                      work(kz+ioff),work(kw+ioff)

!         only the upper triangular part of the one electron embedding 
!         potential matrix is calculated.
          do p = 1, ndimension
            do q = p, ndimension
               
!              first we get the value of the orbital pair (gamma_pq) at point i
!              and multiply it by the embedding potential 
               gamma_pq      = work(kgao+(p-1))*work(kgao+(q-1))
               ep_pq_update  = work(kw+ioff)*gamma_pq

#if !defined (PRG_DIRAC)
!              dalton will use packed storage so we do the array indexing 
!              accordingly. something to be careful about. pq index is given 
!              by the formula for packed storage of the upper triangle.
               pq = p + q*(q-1)/2
               fdemat(pq)  = fdemat(pq) + ep_pq_update 
#else 
!              in dirac, we can work with the full matrix
               fdemat(p,q) = fdemat(p,q) + ep_pq_update 
               fdemat(q,p) = fdemat(p,q)
#endif

            enddo ! loop over q
         enddo ! loop over p
      enddo ! loop over points

      if (debug) then
         call header ('1-electron matrix, after FDE update',-1)
         call printdmat(1,ndimension,fdemat)
      endif

      call memrel ('fde-grid',work,kx,kx,kfree,lfree)
      if (debug) close (lutesting)
      call qexit('GetFDEPotential_FromFile')

      end


      SUBROUTINE FDEINI()
#include "implicit.h"
#include "priunit.h"
#include "fde.h"
!     set defalts for common /FDE/ 
!     filename for embedding potential file
      FDEPOT = 'EMBPOT'
!     print level
      IPRFDE = 0
      END


      SUBROUTINE FDEINP(WORD,FDEHAM)
#include "implicit.h"
#include "priunit.h" 
      PARAMETER (NTABLE = 4)
C
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
#include "fde.h"
C
      logical fdeham,newdef
      DATA TABLE /'.EMBPOT','.PRINT ','.XXXXXX','.XXXXXX'/
C
#include "memint.h"

#if defined (PRG_DIRAC)
      iunit = LUPRI
#else
      iunit = LUPRI
#endif

      NEWDEF = (WORD .EQ. '*FDE   ')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            write (*,*) 'processing line: ',word
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GOTO ( 1, 2, 3, 4), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in FDEINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in FDEINP.')
 1             CONTINUE
! read filename for the embedding potential
               write (*,*) 'reading potential filename' 
               READ(LUCMD,'(A30)') FDEPOT 
               write (*,*) 'file is ',fdepot
               GO TO 100
 2             CONTINUE
! set printlevel 
               write (*,*) 'reading potential filename' 
               READ(LUCMD,*) IPRFDE 
               write (*,*) 'print level is ',iprfde
               GO TO 100
 3             CONTINUE
! add action here
               GO TO 100
 4             CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in FDEINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in FDEINP.')
            END IF
      END IF
  300 CONTINUE

!     options are echoed
      IF(.NOT.FDEHAM) GOTO 999
      write (lupri,*) ' '
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)')
     &   ' Frozen Density Embedding Calculation Setup'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A,I5)')  ' * Print level: ',IPRFDE
      write(LUPRI,'(A,A60)') ' * FDE Potential read from file: ',FDEPOT
      write (lupri,*) ' '
  999 CONTINUE
      end

!
! auxiliary routines
!
#if !defined (PRG_DIRAC) /* e.g. for dalton */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* RYO_GRID */
      SUBROUTINE RYO_GRID(X,Y,Z,W,LUNIT,NR_POINTS)
C=======================================================================
C     read your own grid from file
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                           last revision: february 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D3 = 3.0D0,D4 = 4.0D0)

      DIMENSION X(NR_POINTS),
     &          Y(NR_POINTS),
     &          Z(NR_POINTS),
     &          W(NR_POINTS)

      DO I = 1,NR_POINTS
        READ(LUNIT,*) X(I), Y(I), Z(I), W(I)
      ENDDO

      RETURN
      END
#endif
!
! in case we use dalton 2.0 or dirac without enabling ERI
! we include printdmat.
!
#if !defined (MOD_ERI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* PRINTDMAT */
C aspg
C taken from eri/eri2ftuv.F, added here for convenience
      SUBROUTINE PRINTDMAT( NSIM, NORBT, D )
#include "implicit.h"
#include "priunit.h"

      DIMENSION D(*)

      DO 200 ISIM=1,NSIM
        WRITE( LUPRI, * ) 'MATRIX NUMBER : ', ISIM
          DO 100 I=1,NORBT
            WRITE( LUPRI, 400 ) (D((I-1)*NORBT+J),J=1,NORBT)
 100      CONTINUE
 200  CONTINUE

 400  FORMAT(6F16.12)

      RETURN
      END
C
#endif
