! arguments:
!
!    lufdepot:  unit number of the file containing the embedding potential
!    fdep1e:  matrix to hold the embedding potential. to be added to the one-electron fock matrix
!    work    :  workspace array 
!
      SUBROUTINE GetEmbeddingPotentialFromFile(fdep1e,initmp,work,lwork)
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

      integer lufdepot
      dimension work(lwork), fdep1e(nbast,nbast)
      logical debug, fdepot_is_here

      debug = .true.
      ndimension = nbast 
      fdepot_is_here = .false.

      lutesting=46
      open(lutesting,
     &       FORM   = 'FORMATTED',
     &       STATUS = 'UNKNOWN',
     &       ACCESS = 'SEQUENTIAL',
     &       FILE   =  'TSTPOT')

#include "memint.h"

! zero out matrix to hold the embedding potential, since it willl be accumulated
! pointwise from the grid

      if (initmp) then
         call memget ('REAL',kfdptmp,ndimension*ndimension,
     &                work,kfree,lfree)
         call dzero (work(kdfptmp),ndimension*ndimension)
      endif
c     call dzero (fdep1e,ndimension*ndimension)

! open file containing the embedding potential

      FDEPOT = 'EMBPOT'
      lufdepot = 45
      write (*,*) 'before opening embedding file'
      call flush()

      inquire(FILE=FDEPOT,EXIST=fdepot_is_here)

      if (fdepot_is_here.eq..false.) then
         call quit ('fatal error: embedding potential file not present')
      else
         write (*,*) 'embedding potential file found.'
      endif
         
      open(lufdepot,
     &       FORM   = 'FORMATTED',                                     
     &       STATUS = 'UNKNOWN',                                       
     &       ACCESS = 'SEQUENTIAL',
     &       FILE   =  FDEPOT)

! we start processing the grid file by first getting the number of point in the
! grid

      rewind lufdepot 
      read (lufdepot,*) npoints
      write (*,*) 'number of points to find...',npoints
      call flush()

      if (npoints.gt.0) then
         CALL MEMGET('REAL',KX,NPOINTS,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KY,NPOINTS,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KZ,NPOINTS,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KW,NPOINTS,WORK,KFREE,LFREE)
      else
         call quit ('fatal error: empty embedding potential file')
      endif 

      CALL RYO_GRID(WORK(KX),WORK(KY),WORK(KZ),WORK(KW),
     &              lufdepot,NPOINTs)
      close (lufdepot)

! the size of the grid is then used to allocate workspace for reading the
! actual grid, along with the embedding potential at each point

      call memget ('REAL',kgao ,  ndimension,work,kfree,lfree)
      call memget ('REAL',kbuf ,4*ndimension,work,kfree,lfree)
      call memget ('REAL',kncnt,4*ndimension,work,kfree,lfree)

      nder  = 0
      dolnd = .false.
      dogga = .false.
      iprint = 0
      irepop = 0

      write (lutesting,'(I10)') npoints
! now we are ready to build the embedding one electron matrix.

      do i = 1, npoints

! we get the value of all gaussians at the current grid point 

         ioff = i-1
#if defined (PRG_DIRAC)
         call getsos(work(kgao),dummy,dummy,dummy,
     &               work(kcnt),
     &               work(kx+ioff),
     &               work(ky+ioff),
     &               work(kz+ioff),
     &               work(kbuf),
     &               work,lwork,ndimension,0,.false.,iprint)
#else
         call dftaos(work(kgao),dummy,dummy,dummy,dummy,
     &               work(kcnt),
     &               work(kx+ioff),
     &               work(ky+ioff),
     &               work(kz+ioff),
     &               work,lwork,ndimension,
     &               false.,.false.,DFTHRI,IPRINT)
#endif
c        write (*,'(I10,4F20.10)') i,work(kx+ioff),work(ky+ioff),
      write (lutesting,'(4(E25.18,2X))') work(kx+ioff),work(ky+ioff),
     &                 work(kz+ioff),work(kw+ioff)

! only upper triangular part of the one electron embedding potential
! matrix is calculated.

          do p = 1, ndimension
            do q = p, ndimension
               
               gamma_pq      = work(kgao+(p-1))*work(kgao+(q-1))
               ep_pq_update  = work(kw+ioff)*gamma_pq

               fdep1e(p,q) = fdep1e(p,q) + ep_pq_update 
               fdep1e(q,p) = fdep1e(p,q)
            enddo
         enddo
      enddo ! loop over points


      if (debug) then
         call header ('1-electron embedding matrix',-1)
! does not know....
!        if (initmp) then
!           call printdmat(1,ndimension,work(kfdptmp))
!        else
!           call printdmat(1,ndimension,fdep1e)
!        endif
      endif
 
      if (initmp) then
!copies the embedding potential to the 1e fock matrix
          do p = 1, ndimension
            poff = (p-1)*ndimension
            do q = p, ndimension
               qoff = poff + q
               fdep1e(p,q) = fdep1e(p,q) + work(kfdptmp+qoff) 
               fdep1e(q,p) = fdep1e(p,q)
            enddo
         enddo
         call memrel ('fde-1',work,kfdptmp,kfdptmp,kfree,lfree)
      else
         call memrel ('fde-1',work,kx,kx,kfree,lfree)
      endif

      close (lutesting)
      end


      subroutine FDEINP(WORD,WORK,LWORK)
#include "fde.h"
      write (*,*) ' '
      write (*,*) 'Activating Frozen Density Embedding menu'
      write (*,*) ' '
      FDEPOT = 'EMBPOT'
      end

#if !defined (PRG_DIRAC)
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
