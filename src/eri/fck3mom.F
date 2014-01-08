!dalton_copyright_start
!      Copyright (c) 2013 by the authors of Dalton (see below).
!      All Rights Reserved.
!   
!      The source code in this file is part of
!      "Dalton, a molecular electronic structure program,
!       Release DALTON2013 (2013), see http://daltonprogram.org"
!   
!      This source code is provided under a written licence and may be
!      used, copied, transmitted, or stored only in accord with that
!      written licence.
!   
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the licence holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the licence holder may not take any part of Dalton,
!      including modified files, with him/her, unless that person has
!      obtained his/her own licence.
!   
!      For further information, including how to get a licence, see:
!         http://daltonprogram.org
!dalton_copyright_end

#ifdef MOD_FMM
      SUBROUTINE SETMM_OR_JFF(DMAT,NDMAT,DAUX,IODBCLS,IODBTCH,RODBTCH,
     &                 WORK,LWORK,IPRINT1,taskstr,DJFARFIELD)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h" 
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#ifdef MOD_FMM
#include "cbifmm.h" 
#endif
#include "nuclei.h"
#include "inforb.h"
C
#ifdef MOD_PBC
#include "pbc.h"
#endif
#ifdef MOD_DNF
#include "densfit.h"
#endif
      DIMENSION DMAT(*),DAUX(*),IODBCLS(5,*),IODBTCH(2,*),RODBTCH(7,*),
     &     WORK(*), DJFARFIELD(*)
      character*(*) taskstr

C     This routine can be called to do three different tasks:
C     1. taskstr == 'FMM' means ordinary FMM
C     2. taskstr == 'PBCINDEXOVL' means that moments for the individual
C        overlap distributions \Omega^{0l}_{ab} associated with
C        J^{0l}_{ab} should be generated.
C        
C        The moments are then immediately used to compute the far-field
C        contribution to the J^{0l} matrix after which they are thrown
C        away. These contributions are added to J^{0l} so this matrix
C        needs to be passed as input (it will only be referenced if
C        taskstr == 'PBCINDEXOVL').
C     3. taskstr == 'PBCDUMMYOVL' means that moments for the overlaps
C        \Omega^{0j}_{cd} should be generated and summed up. I.e.,
C        only electronic moments corresponding to a fixed j are accounted
C        for.
C     4. taskstr == 'PBCNUCLMMOM' means that the moments from the
C        nuclei (and fictitious charges) are added to the total
C        moments of the cell. The total electronic moments are needed
C        to determine the fictitious charges, so this option must
C        be used no sooner than after all electronic contributions
C        have been added up.
C
C     These moments will be used to calculate the Coulomb interaction
C
C     J^{0l;FAR}_{ab} = \sum_m \sum_{cdj} D^{0j}_{cd} *
C                          * (\Omega^{-C,l-C}_{ab}|\Omega^{0j}_{cd})
C                     = q^{abl} T^{lattice} Q^{cell}
C
C     where C is the center of the cell that contains \Omega^{0l}_{ab}.
C     Each overlap between primitive GTOs gets its individual
C     translation C (just as in the integral evaluation for the near-field).
C     Assumption:  The cell that is represented explicitly is centered
C                  at 0. /ErikT
      if (taskstr .eq. 'FMM') then
C        get order of expansion from FMM common block
         LMAX_GEN = LMAX
#ifdef MOD_PBC
      else if ((taskstr .eq. 'PBCINDEXOVL') .or.
     &         (taskstr .eq. 'PBCDUMMYOVL') .or.
     &         (taskstr .eq. 'PBCNUCLMMOM')) then
         if (.not. (dopbc .and. pbcactive)) then
            call quit('Called to do PBC task with PBC disabled.')
         end if
C        get order of expansion from PBC common block
         LMAX_GEN = qlmax_pbc
#endif
      else
         call quit('Invalid task in SETMM.')
      end if

      LTOP = MAX( LMAX_GEN, MAXAONHK-1 )

      MAXMMNHK = LTOP + 1
      MAXMMKHK = MAXMMNHK * 2 - 1
      MAXMMKCK = MAXMMNHK * ( MAXMMNHK + 1 ) / 2  

      LUINTM = -1
      CALL GPOPEN(LUINTM,'MM_DATA','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
C
C  This is just to force SPHMAT to generate spherical 
C  transformation matrices to the sufficient order:
C

      LSIZECSQE = 0 
  
      DO L=0,LTOP
       LSIZECSQE = LSIZECSQE + (2*L+1)*(L+1)*(L+2)/2
      END DO

      KLMNPWR  = 1
      KLMNSYM  = KLMNPWR  + ( 3*MAXMMNHK*MAXMMKCK - 1)/IRAT + 1
      KCSQ     = KLMNSYM  + (16*MAXMMNHK*MAXMMKCK - 1)/IRAT + 1
      KCSQE    = KCSQ     + NCSQ1*NCSQ2
      KICSQADR = KCSQE    + LSIZECSQE
      KWRK     = KICSQADR + LMAX_GEN + 1

      LNEED = KWRK - 1

      IF ( LNEED .GT. LWORK ) THEN
         WRITE (LUPRI,'(/,2(2X,A,I15,/))')
     &      'Memory requirements total     : ',LNEED,
     &      'Memory available              : ',LWORK
         CALL QUIT('Insufficient memory in SETMM (FCK3)')
      END IF

      LWRK  = LWORK - 1 

C
C    Prepare carteesian-spherical transformation matrices 
C    and index array LMN for packing carteesian components
C    =====================================================
C

      CALL CARPOW

      CALL SPHMAT(WORK(KCSQ),NCSQ1,NCSQ2,0,.FALSE.,
     &            WORK(KWRK),LWRK,IPRINT1)
 
      CALL SETCSQE( WORK(KICSQADR), WORK(KCSQE), LMAX_GEN,
     &              WORK(KWRK),LWRK,IPRINT1)

C
C     Calculate the multipole moments and write to file
C     =================================================
C
      CALL SETMM1(NDISTR,NMOMNT,DMAT,DAUX,NDMAT,IODBCLS,IODBTCH,
     &            RODBTCH,WORK(KCSQ),WORK(KCSQE),WORK(KICSQADR),
     &            WORK(KLMNPWR),WORK(KLMNSYM),MAXMMNHK,MAXMMKCK,
     &            WORK(KWRK),LWRK,LUINTM,IPRINT1,taskstr,DJFARFIELD)

      WRITE (LUINTM) -1,0, 0,0, 0,0,0, 0, 0D0, 0D0,0D0,0D0, 0D0,0D0
C
C     Nuclear information
C     =================== 
C
      DO I = 1, NUCIND
         WRITE (LUINTM) CHARGE(I), CORD(1,I), CORD(2,I), CORD(3,I)
      END DO
C
      CALL GPCLOSE(LUINTM,'KEEP')

C     Write counts to a separate file for fast access.
      LUINTM = -1
      CALL GPOPEN(LUINTM,'MM_CNTS','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      WRITE (LUINTM) LMAX_GEN, NBAST, NDISTR, NUCIND, NMOMNT,
cSimen Added for density fitting
     & NBASEAUX, DOFIT .AND. LHSFIT, DOFIT .AND. RHSFIT
cSimen
#if 0
      WRITE (LUPRI,'(/,2X,A,2I8,2X,F12.6)')
     &     ' Number of nonzero distributions in SETMM (FCK3):',NDISTR
      WRITE (LUPRI,'(2X,A,2I12,2X,F12.6)')
     &     ' Number of moments written in SETMM: (FCK3):',NMOMNT
#endif
C
      CALL GPCLOSE(LUINTM,'KEEP')

      RETURN
      END 

      SUBROUTINE SETMM1(NDISTR,NMOMNT,
     &                  DMAT,DAUX,NDMAT,IODBCLS,IODBTCH,RODBTCH,CSQ,
     &                  CSQE,ICSQADR,LMNPWR,LNMSYM,MAXMMNHK,MAXMMKCK,
     &                  WORK,LWORK,LUINTM,IPRINT1,taskstr,DJFARFIELD)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
#ifdef MOD_FMM
#include "cbifmm.h"
#endif
C
#include "nuclei.h"
#ifdef MOD_PBC
#include "pbc.h"
#endif
#include "maxmom.h"

C     Nfict should always be 4 - it has been introduced as a
C     parameter only for readability / ErikT & Filip
      parameter ( Nfict = 4 )
      DIMENSION DMAT(*),DAUX(*),IODBCLS(5,*),IODBTCH(2,*),RODBTCH(7,*),
     &          WORK(*),
     &          LMNPWR(MAXMMKCK,MAXMMNHK,3),
     &          LMNSYM(MAXMMKCK,MAXMMNHK,2,0:7),CSQ(*),CSQE(*),
     &          ICSQADR(*),KCCF(2,2),
     &          DJFARFIELD(*)
      dimension tvec1(3), tvec2(3)
      dimension dipmom(3), fictq(Nfict), fictpos(3,Nfict)
      character*(*) taskstr
      integer sph_off, car_off
      logical reset_flag
C
#ifdef MOD_PBC
C     declare a logical external function:
      logical pbc_3dims_active
#endif
C
      IELECTR = 1
C
C  here is maxmmkcksize ... sum( (j+1)*(j+2)/2 )
C          maxmmkhksize ... sum( 2*j+1 )
C

      MAXMMKCKSIZE = MAXMMNHK * ( MAXMMNHK + 1 ) * ( MAXMMNHK + 2 ) / 6
      MAXMMKHKSIZE = MAXMMNHK * MAXMMNHK

C
C  SMM...  - multipole moment integrals
C      HCP   - Herminte Cartesion primitive
C      CP    - Cartesian primitive
C      CC    - Cartesion contracted 
C      SC    - spherical contracted
C
      KETIJ    = 1
      KPREF12  = KETIJ    + MAXODPRF * MAXODNHK*(MAXAONHK*MAXAONHK)*3*2
      KIPR12   = KPREF12  + MAXODPRF
      KCOOR12  = KIPR12   + MAXODPRF * 2 / IRAT + 1 
      KEXP12   = KCOOR12  + MAXODPRF * 3
      KHPI     = KEXP12   + MAXODPRF * 3
      KP1      = KHPI     + MAXODPRF
      KP2      = KP1      + MAXODPRF * 3
      KCOORP   = KP2      + MAXODPRF * 3
      KSMMHCP  = KCOORP   + MAXODPRF * 3
      KSMMCP   = KSMMHCP  + MAXODPRF * MAXMMNHK      * MAXMMNHK * 3
      KSMMCC   = KSMMCP   + MAXODPRF * (MAXAOKCK**2) * MAXMMKCKSIZE 
      KSMMSC   = KSMMCC   + MAXODCTF * (MAXAOKCK**2) * MAXMMKCKSIZE
      KINDDSTR = KSMMSC   + MAXODCTF * (MAXAOKHK**2) * MAXMMKHKSIZE
      KLAST    = KINDDSTR + MAXODCTF * (MAXAOKHK**2)

      LNEED = KLAST - 1

      IF ( LNEED .GT. LWORK ) THEN
         WRITE (LUPRI,'(/,2(2X,A,I15,/))')
     &      'Memory requirements total     : ',LNEED,
     &      'Memory available              : ',LWORK
         CALL QUIT('Insufficient memory in SETMM (FCK3)')
      END IF

      NMOMNT = 0 
      NDISTR = 0

#ifdef MOD_PBC
      if (taskstr .eq. 'PBCNUCLMMOM') then
C        Here we add the nuclear moments to the total moments.
C        It is assumed that the electronic moments have already been
C        computed, so that we now add the last contributions
C        and do error checking on the total charge and dipole moment.
         JMAXE = qlmax_pbc

         if (JMAXE .gt. MXQNM) then
            write(LUPRI,'(A,I5,A,I5,A)') 'LMAX for PBC is ',LMAX,
     &                  ' but ',MXQNM,' is the maximum.'
            write(LUPRI,*) 'Get around this by changing MXQNM in',
     &                  ' the include file maxmom.h'
            call quit('LMAX is too large.')
         end if

        write(LUPRI,*) 'debug: Charge from electronic moments is ',
     &        cell_mulmom(1)
        write(LUPRI,*) '  el. dip.mom. = ',cell_mulmom(2),
     &       cell_mulmom(3),cell_mulmom(4)
C        Loop over nuclei in cell 0 and add up their contributions
C        to the total multipole moments that represent cell 0.
         reset_flag = .true.
         do inuc = 1, NUCIND
            write(LUPRI,*) 'debug: n',inuc,JMAXE,CHARGE(inuc)
            call flshfo(LUPRI)
C           Nuclei data in CHARGE(*) and CORD(1:3,*).
            call calc_pcharge_mom(WORK(KSMMCC),CHARGE(inuc),
     &                            CORD(1,inuc),JMAXE,reset_flag)
            reset_flag = .false.
         end do

         write(LUPRI,*) 'debug: Charge from Cartesian nuc.mom. is ',
     &        WORK(KSMMCC)

C        Determine the total dipole moment and solve for the
C        sizes of the fictitious charges necessary to cancel
C        the dipole moment.
         dipmom(1) = WORK(KSMMCC+1) + cell_mulmom(4)    !3
         dipmom(2) = WORK(KSMMCC+2) + cell_mulmom(2)    !4
         dipmom(3) = WORK(KSMMCC+3) + cell_mulmom(3)    !2
         if (pbc_3dims_active()) then
            write(LUPRI,*) 'Cartesian dipole moment of unit cell '//
     &          '(excluding fictitious charges):'
            write(LUPRI,'(A,3G12.5,A)') '    dip.mom. = (',
     &          dipmom(1),dipmom(2),dipmom(3),')'
C           We only need fictitious charges for 3D periodic systems.
C           N.B. They can be used for other systems too, but this
C           code works only for 3D periodicity, and they are not
C           needed for lower periodicity, anyway.
            call pbc_calc_fictcharges(dipmom)
            call pbc_get_fictcharges(fictq,fictpos)
C           Loop over fictitious charges and collect their moments.
            do I = 1, Nfict
               write(LUPRI,*) 'debug: f',I,JMAXE
               call calc_pcharge_mom(WORK(KSMMCC),fictq(I),
     &                            fictpos(1,I),JMAXE,.false.)
            end do
         else
            write(LUPRI,*) 'Cartesian dipole moment of unit cell:'
            write(LUPRI,'(A,3G12.5,A)') '    dip.mom. = (',
     &          dipmom(1),dipmom(2),dipmom(3),')'
         end if
C        Transform to spherical basis
         call dzero(WORK(KSMMSC),(1+JMAXE)**2)
         sph_off = 0
         car_off = 0
C        N.B. Here the index of ICSADR range from 1 to LMAX+1,
C             elsewhere it ranges from 0 to LMAX. / ErikT
         KCSQ1  = ICSQADR(0 + 1)
         KCSQ2  = ICSQADR(0 + 1)
         do je = 0, JMAXE
            KCSQE = ICSQADR(je + 1)
            NHKTE = je + 1
            KHKTE = NHKTE * 2 - 1
            KCKTE = NHKTE * ( NHKTE + 1 ) / 2
            call MMCAR2MMSPH(WORK(KSMMSC+sph_off),WORK(KSMMCC+car_off),
     &                  CSQ(KCSQ1), CSQ(KCSQ2), CSQE(KCSQE),
     &                  1,1,KHKTE,1,1,KCKTE,1,IPRINT1 ) 
            sph_off = sph_off + KHKTE
            car_off = car_off + KCKTE
         end do

         write(LUPRI,*) 'debug: Charge from spherical nuc.mom. is ',
     &        WORK(KSMMSC)

C        Add the contributions to the multipole moments of cell 0.
         call pbc_addupcellmom(WORK(KSMMSC),DMAT,
     &              -1,-1,0,0,JMAXE,1,1,
     &              -1,-1,.true.,IPRINT1)
C        We now have all the multipole moments for a unit cell in the
C        common block variable cell_mulmom. It only remains to put in
C        some prefactors to make the moments compatible with the W and
C        T tensors in the PBC-FMM code.
         call pbc_redefine_q(cell_mulmom,JMAXE)
C
         if (pbc_3dims_active()) then
           write(LUPRI,*) 'Spherical dipole moment of unit cell '//
     &          '(including fictitious charges):'
           write(LUPRI,'(A,3G12.5,A)') '    dip.mom. = (',
     &          cell_mulmom(1),cell_mulmom(2),cell_mulmom(3),')'
C          Error check: The total dipole moment (and charge!)
C          should now be zero
           if (abs(cell_mulmom(1)) .gt. mm_q_loss_thres) then
             write(LUPRI,'(A,F10.6)')
     &             'Error: Cell multipole moments represent'//
     &             ' a system with total charge ',cell_mulmom(1)
             write(LUPRI,'(A,D)') '       Maximum allowed charge'//
     &            ' loss is ',mm_q_loss_thres
             call quit('PBC-FMM: non-zero charge')
             call quit('PBC-FMM: non-zero charge')
           else if (abs(cell_mulmom(2)) .gt. 1.0D-8) then
             call quit('PBC-FMM: non-zero dip.mom. (z)')
           else if (abs(cell_mulmom(3)) .gt. 1.0D-8) then
             call quit('PBC-FMM: non-zero dip.mom. (x/2)')
           else if (abs(cell_mulmom(4)) .gt. 1.0D-8) then
             call quit('PBC-FMM: non-zero dip.mom. (y/2)')
           end if
C          Set total charge and dipole moment of the unit cell to zero
           cell_mulmom(1) = 0.0D0 ! charge
           cell_mulmom(2) = 0.0D0 ! dip.mom. (z)
           cell_mulmom(3) = 0.0D0 ! dip.mom. (x/2)
           cell_mulmom(4) = 0.0D0 ! dip.mom. (y/2)
         else
           if (abs(cell_mulmom(1)) .gt. mm_q_loss_thres) then
             write(LUPRI,'(A,F10.6)')
     &             'Error: Cell multipole moments represent'//
     &             ' a system with total charge ',cell_mulmom(1)
             write(LUPRI,'(A,D)') '       Maximum allowed charge'//
     &            ' loss is ',mm_q_loss_thres
             call quit('PBC-FMM: non-zero charge')
           end if
C          Set total charge of the unit cell to zero
           cell_mulmom(1) = 0.0D0    ! charge
         end if
C        now exit - no electronic contributions should be done with this
C        option
         return
      end if
#endif /* MOD_PBC */
      
      DO ICLS = 1, NODBCLS

       ICENT1  = IODBCLS( 1, ICLS ) 
       ICENT2  = IODBCLS( 2, ICLS ) 
       JMAX12  = IODBCLS( 3, ICLS ) - 1 
       NHKT12  = IODBCLS( 3, ICLS ) 
       KSTBTCH = IODBCLS( 4, ICLS ) 
       NBTCH   = IODBCLS( 5, ICLS ) 

       if (taskstr .eq. 'FMM') then
          tvec1(1) = 0.0
          tvec1(2) = 0.0
          tvec1(3) = 0.0
          tvec2(1) = 0.0
          tvec2(2) = 0.0
          tvec2(3) = 0.0
#ifdef MOD_PBC
       else if (taskstr .eq. 'PBCINDEXOVL') then
C         The overlap is between an orbital in cell 0 and
C         an orbital in cell l, but we will later move the overlap
C         center so that it is in cell 0 (the system is translation
C         invariant so this only affects the division into near-field
C         and far-field parts).
          tvec1(1) = 0.0    ! the 1st index orbital is in cell 0
          tvec1(2) = 0.0
          tvec1(3) = 0.0
          tvec2(1) = rspc_vec_std(1)  ! the 2nd index orb. is in cell l
          tvec2(2) = rspc_vec_std(2)
          tvec2(3) = rspc_vec_std(3)
       else if (taskstr .eq. 'PBCDUMMYOVL') then
          tvec1(1) = 0.0    ! the 1st dummy orbital is in cell 0
          tvec1(2) = 0.0
          tvec1(3) = 0.0
          tvec2(1) = nf2_vec_std(1)   ! the 2nd dummy orb. is in cell j
          tvec2(2) = nf2_vec_std(2)
          tvec2(3) = nf2_vec_std(3)
#endif /* MOD_PBC */
       else
          call quit('Invalid taskstr in SETMM1.')
       end if
       CALL GETR12( NSYMOPR, NSYMOPS, ICENT1, ICENT2, 2,
     &      tvec1, tvec2 )

       DO IBTCH = KSTBTCH, ( KSTBTCH+NBTCH-1)

        NMOM = 0 

        IAOB1 = IODBTCH( 1, IBTCH ) 
        IAOB2 = IODBTCH( 2, IBTCH ) 

        CALL GETAOB(NPRF1,NCTF1,NTYP1,NHKT1,JMAX1,KHKT1,KCKT1,
     &              KCSQ1,KEXP1,KCCF(1,1),IAOB1)

        CALL GETAOB(NPRF2,NCTF2,NTYP2,NHKT2,JMAX2,KHKT2,KCKT2,
     &              KCSQ2,KEXP2,KCCF(1,2),IAOB2)

        NCTF12 = NCTF1 * NCTF2

        CALL SCREENPRIM12(NPRF12,RODBTCH(1,IBTCH),WORK(KIPR12),.TRUE.,
     &                    NPRF1,NCTF1,NTYP1,NPRF2,NCTF2,NTYP2,
     &                    EXPBTCH(KEXP1),EXPBTCH(KEXP2),
     &                    CCFBTCH(KCCF(1,1)),CCFBTCH(KCCF(2,1)),
     &                    CCFBTCH(KCCF(1,2)),CCFBTCH(KCCF(2,2)),
     &                    R122(0,0),COOR1(1,0),COOR2(1,0),1,IPRINT1)

C        WRITE(6,*) 'data:',NPRF12, NCTF12
C        WRITE(LUPRI,*) 'data:',NPRF12, NCTF12

        CALL SETETIJ(WORK(KETIJ),WORK(KPREF12),
     &               R12(1,0,0),R122(0,0),
     &               WORK(KIPR12),EXPBTCH(KEXP1),EXPBTCH(KEXP2),
     &               WORK(KCOOR12),WORK(KEXP12),WORK(KHPI),WORK(KP1),
     &               WORK(KP2),NPRF1,NHKT1,JMAX1,
     &               NPRF2,NHKT2,JMAX2,NPRF12,NHKT12,IELECTR,IPRINT1)

        J1 = JMAX1 
        J2 = JMAX2 
        if (taskstr .eq. 'FMM') then
C          get order of expansion from FMM common block
           JMAXE = LMAX
#ifdef MOD_PBC
        else if ((taskstr .eq. 'PBCINDEXOVL') .or.
     &          (taskstr .eq. 'PBCDUMMYOVL')) then
C          get order of expansion from PBC common block
           JMAXE = qlmax_pbc
        else
#endif /* MOD_PBC */
           call quit('Invalid task in SETMM.')
        end if

C       Some error-checking here / ErikT and Filip
        if (JMAXE .gt. MXQNM) then
           write(LUPRI,'(A,I5,A,I5,A)') 'LMAX for PBC is ',LMAX,
     &          ' but ',MXQNM,' is the maximum.'
           write(LUPRI,*) 'Get around this by changing MXQNM in',
     &          ' the include file maxmom.h'
           call quit('LMAX is too large.')
        end if

        DO ICOMP1 = 1, NTYP1
         IF ( NTYP1 .EQ. 2 ) CALL SETJKKK(J1,KCSQ1,KCKT1,KHKT1,ICOMP1)

         DO ICOMP2 = 1, NTYP2 
          IF ( NTYP2 .EQ. 2 ) CALL SETJKKK(J2,KCSQ2,KCKT2,KHKT2,ICOMP2)

C
C     Calculate the multipole moment integrals
C     ========================================
C           

          CALL SETMM2(WORK(KSMMSC), WORK(KSMMCC), WORK(KSMMCP), 
     &         WORK(KSMMHCP), WORK(KETIJ),WORK(KIPR12),
     &         WORK(KPREF12),EXPBTCH(KEXP1),EXPBTCH(KEXP2),
     &         LMNPWR,WORK(KCOORP),COOR1(1,0),COOR2(1,0),
     &         WORK(KEXP12),CSQ,CSQE,ICSQADR,
     &         CCFBTCH(KCCF(ICOMP1,1)),CCFBTCH(KCCF(ICOMP2,2)),
     &         RODBTCH(1,IBTCH),MAXMMNHK,MAXMMKCK,
     &         J1,JMAX1,KHKT1,KCKT1,NPRF1,NCTF1,
     &         J2,JMAX2,KHKT2,KCKT2,NPRF2,NCTF2,
     &         JMAXE,JMAX12,NPRF12,NCTF12,IPRINT1,taskstr)

          if (taskstr .eq. 'FMM') then
C            Save in file that will be read by the FMM code
             CALL PRINTMM(WORK(KSMMSC),DMAT,NDMAT,DAUX,RODBTCH(1,IBTCH),
     &                 WORK(KINDDSTR),IAOB1,IAOB2,J1,J2, JMAX1,
     &                 JMAX2,JMAXE,NCTF1,NCTF2,ICOMP1,ICOMP2,
     &                 IBTCH, NDISTR,NMOM,LUINTM,IPRINT1)
#ifdef MOD_PBC
          else if (taskstr .eq. 'PBCINDEXOVL') then
C            Instead of storing the moments q^{abl} of the index overlaps,
C            we add their contribution
C                  J^{0l;far}_{ab} = q^{abl} Tlattice Qcell
C            to the matrix DJFARFIELD.
             call pbc_add_to_JFF(WORK(KSMMSC),
     &              IAOB1,IAOB2,J1,J2,JMAXE,NCTF1,NCTF2,
     &              ICOMP1,ICOMP2,IPRINT1,DJFARFIELD)
C
C     No need to store moments anymore / ErikT & Filip
C
C            Store the moments of the index overlaps individually.
C             call pbc_printmm(WORK(KSMMSC),
C     &              IAOB1,IAOB2,J1,J2,JMAXE,NCTF1,NCTF2,
C     &              ICOMP1,ICOMP2,LUFILE,IPRINT1)
          else if (taskstr .eq. 'PBCDUMMYOVL') then
C            Add the contributions to the multipole moments of cell 0.
             call pbc_addupcellmom(WORK(KSMMSC),DMAT,
     &              IAOB1,IAOB2,J1,J2,JMAXE,NCTF1,NCTF2,
     &              ICOMP1,ICOMP2,.false.,IPRINT1)
#endif /* MOD_PBC */
          else
             call quit('Invalid task in SETMM1')
          end if
         END DO
        END DO

       IF ( NMOM .GT. 0 ) THEN
        NMOMNT = NMOMNT + NMOM 
       END IF

       END DO
      END DO

      RETURN
      END 

C     This subroutine adds to cmm the Cartesian multipole moments
C     (expansion center: 0) from a point charge with charge q and
C     position pos. / ErikT
      subroutine calc_pcharge_mom(cmm,q,pos,lmax,reset_flag)
#include "implicit.h"   
#include "priunit.h"   
#include "pbc.h"
      dimension cmm(*), pos(*)
      integer lmax
      logical reset_flag
      dimension pos_pow(3,max_qlmax_pbc+1)
      integer i,j,k,ltot,index_ijk,icomp

      write(LUPRI,'(A,F10.4,A,3F10.4,A)') 'calc_pcharge_mom: q = ',q,
     &     ' at pos = (',pos(1),pos(2),pos(3),')'
      call flshfo(LUPRI)


      if (lmax .lt. 0) then
         call quit('calc_pcharge_mom: negative lmax.')
      else if (lmax .gt. max_qlmax_pbc) then
         call quit('calc_pcharge_mom: lmax is too large.')
      end if


C     Compute x^n, y^n, z^n for n = 0,...,lmax
      pos_pow(1,1) = 1.0D0
      pos_pow(2,1) = 1.0D0
      pos_pow(3,1) = 1.0D0
      do i = 1,lmax
         do icomp = 1,3
            pos_pow(icomp,i+1) = pos_pow(icomp,i) * pos(icomp)
         end do
      end do

C     The components are generated in the order (0 0 0), (1 0 0),
C     (0 1 0), (0 0 1), (2 0 0), (1 1 0), (1 0 1), (0 2 0), ...

C     Loop over Cartesian total "ang. mom."
      index_ijk = 1
      do ltot = 0, lmax, 1
         do i = ltot, 0, -1
            do j = ltot - i, 0, -1
               k = ltot - i - j
               if (reset_flag) then
                  cmm(index_ijk) = q * pos_pow(1,i+1)
     &                 * pos_pow(2,j+1) * pos_pow(3,k+1)
               else
                  cmm(index_ijk) = cmm(index_ijk) + q * pos_pow(1,i+1)
     &                 * pos_pow(2,j+1) * pos_pow(3,k+1)
               end if
               index_ijk = index_ijk + 1
            end do
         end do
      end do

      return
      end

      SUBROUTINE SETMM2(SMMSC,SMMCC,SMMCP,SMMHCP,ETIJ,IPR12,
     &                  PREF12,EXPBT1,EXPBT2,LMNPWR,
     &                   COORP,COORC,COORD,EXP12,CSQ,CSQE,
     &                   ICSQADR,
     &                   CCF1,CCF2,RODP,MAXMMNHK,MAXMMKCK,
     &                   J1,JMAX1,KHKT1,KCKT1,NPRF1,NCTF1,
     &                   J2,JMAX2,KHKT2,KCKT2,NPRF2,NCTF2,
     &                   JMAXE,JMAX12,NPRF12,NCTF12,IPRINT1,taskstr)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxmom.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#ifdef MOD_FMM
#include "cbifmm.h"
#endif
#include "pi.h"
#include "xyzpow.h"

      character*(*) taskstr
      DIMENSION SMMSC(NCTF1,NCTF2,KHKT1,KHKT2,*),
     &          SMMCC(NCTF1,NCTF2,KCKT1,KCKT2,*),
     &          SMMCP(NPRF12,KCKT1,KCKT2,*),
     &          SMMHCP(NPRF12,0:JMAXE,0:JMAXE,3),
     &          ETIJ(NPRF12,0:JMAX12,0:JMAX1,0:JMAX2,3,2),
     &          PREF12(NPRF12),EXPBT1(NPRF1),EXPBT2(NPRF2),
     &          LMNPWR(MAXMMKCK,MAXMMNHK,3),
     &          COORP(NPRF12,3),COORC(3),COORD(3),EXP12(NPRF12,3),
     &          CSQ(*),CSQE(*),ICSQADR(0:JMAXE),
     &          CCF1(NPRF1,NCTF1),CCF2(NPRF2,NCTF2), 
     &          IPR12(2,NPRF12),RODP(5)

      INTEGER E, T,T1,T2,T3
      integer nf_pos_index(3)
      dimension pos_ovl(3), pos_nf_lat(3), pos_nf_std(3)

      NHKT1 = J1    + 1
      NHKT2 = J2    + 1
      NHKTE = JMAXE + 1

      MAXNHKTE = JMAXE + 1

      MAXKCKTESIZE = MAXNHKTE * ( MAXNHKTE + 1 ) * ( MAXNHKTE + 2 ) / 6
      MAXKHKTESIZE = MAXNHKTE * MAXNHKTE      

      CALL DZERO(SMMHCP,NPRF12*(JMAXE+1)*(JMAXE+1)*3)
      CALL DZERO(SMMCP,NPRF12*KCKT1*KCKT2*MAXKCKTESIZE)
      CALL DZERO(SMMCC,NCTF1 *NCTF2*KCKT1*KCKT2*MAXKCKTESIZE)
      CALL DZERO(SMMSC,NCTF1 *NCTF2*KHKT1*KHKT2*MAXKHKTESIZE)

      DO I = 1, NPRF12
       IPRFA = IPR12(1,I)
       IPRFB = IPR12(2,I)

       EXPA  = EXPBT1(IPRFA)
       EXPB  = EXPBT2(IPRFB)

       EXPPI = EXP12(I,3)

       SMMHCP00 =  sqrt( PI * EXPPI )

C
C    We will calculate the moments relative to the 
C    center of the contracted OD ( stored in RODP ).
C
       if (taskstr .eq. 'FMM') then
          DO IX = 1, 3
             COORP(I,IX) = (COORC(IX)*EXPA+COORD(IX)*EXPB)*EXPPI
     &                     - RODP(IX)
             SMMHCP(I,0,0,IX) = SMMHCP00
             IF ( IPRINT1 .GE. 5 ) THEN
                WRITE(LUPRI,'(A,2I5,6F12.6)') 'COORP - ',I,IX,EXPA,EXPB, 
     &               EXPPI,RODP(IX), COORC(IX)*EXPA , COORD(IX)*EXPB
             END IF
          END DO
#ifdef MOD_PBC
       else if (taskstr .eq. 'PBCINDEXOVL') then
C         move the index overlap center into cell 0
          DO IX = 1, 3
             pos_ovl(IX) = (COORC(IX)*EXPA+COORD(IX)*EXPB)*EXPPI
             SMMHCP(I,0,0,IX) = SMMHCP00
          END DO
          call pbc_nf_center(pos_ovl, nf_pos_index)
          pos_nf_lat(1) = dfloat(nf_pos_index(1))
          pos_nf_lat(2) = dfloat(nf_pos_index(2))
          pos_nf_lat(3) = dfloat(nf_pos_index(3))
          call pbc_lat2std_coord(pos_nf_lat, pos_nf_std)
          COORP(I,1) = pos_ovl(1) - pos_nf_std(1)
          COORP(I,2) = pos_ovl(2) - pos_nf_std(2)
          COORP(I,3) = pos_ovl(3) - pos_nf_std(3)
       else if (taskstr .eq. 'PBCDUMMYOVL') then
C          In the PBC case we want moments w.r.t. to 0 (center of
C          the central cell). / ErikT
           DO IX = 1, 3
              COORP(I,IX) = (COORC(IX)*EXPA+COORD(IX)*EXPB)*EXPPI
              SMMHCP(I,0,0,IX) = SMMHCP00
           END DO
        else
#endif /* MOD_PBC */
           call quit('Invalid task in SETMM2.')
        end if
      END DO

C
C    Generate the primitive hermite multipole-moment integrals by recursion
C    ======================================================================

      DO IX = 1, 3       
       DO E = 1, JMAXE
     
        DO I = 1, NPRF12 
         SMMHCP(I,E,0,IX) = COORP(I,IX)        * SMMHCP(I,E-1,0  ,IX) 
     &                    + 0.5D0 * EXP12(I,3) * SMMHCP(I,E-1,1  ,IX)
         SMMHCP(I,E,E,IX) =    E               * SMMHCP(I,E-1,E-1,IX) 
        END DO

        DO T = 1, E - 1
         DO I = 1, NPRF12 
          SMMHCP(I,E,T,IX) =     T             * SMMHCP(I,E-1,T-1,IX) 
     &                    +   COORP(I,IX)      * SMMHCP(I,E-1,T  ,IX)
     &                    + 0.5D0 * EXP12(I,3) * SMMHCP(I,E-1,T+1,IX)
         END DO
        END DO
       END DO
      END DO

      IF ( IPRINT1 .GE. 5 ) THEN
       DO I=1,NPRF12
        WRITE(LUPRI,'(A,I5,3F12.6)') 'I PCX PCY PCZ :', 
     &              I, COORP(I,1), COORP(I,2), COORP(I,3)
       END DO

       DO IX=1,3
        DO E=0,JMAXE
         DO T=0,E
          DO I=1,NPRF12
           WRITE(LUPRI,'(A,4I4,F12.6)') 'SMMHCP',
     &                  I,E,T,IX,SMMHCP(I,E,T,IX)
          END DO
         END DO
        END DO       
       END DO
      END IF

C
C    Contract the hermite mm integrals with the expansion coefficients
C    =================================================================

      DO KC1 = 1, KCKT1

       L1 = NHKT1 - ISTEP(KC1)
       M1 = MVAL(KC1)
       N1 = NVAL(KC1)

       DO KC2 = 1, KCKT2

        L2 = NHKT2 - ISTEP(KC2)
        M2 = MVAL(KC2)
        N2 = NVAL(KC2)

        KCE0 = 1
        DO JE = 0, JMAXE 
         NHKTE   = JE + 1 
         KCKTE   = NHKTE * ( NHKTE + 1 ) / 2
         KCELAST = KCE0 + KCKTE - 1

         DO KCE = KCE0, KCELAST

          LE = NHKTE - ISTEP(KCE-KCE0+1)
          ME = MVAL(KCE-KCE0+1)
          NE = NVAL(KCE-KCE0+1)

          DO T1 = 0, MIN( JMAX12, NHKTE-1 )
          DO T2 = 0, MIN( JMAX12, NHKTE-1 )
          DO T3 = 0, MIN( JMAX12, NHKTE-1 )

           DO IPRF12 = 1, NPRF12 
            SMMCP(IPRF12,KC1,KC2,KCE) = SMMCP(IPRF12,KC1,KC2,KCE)
     &                                + PREF12(IPRF12)
     &             * SMMHCP(IPRF12,LE,T1,1) * ETIJ(IPRF12,T1,L1,L2,1,1)
     &             * SMMHCP(IPRF12,ME,T2,2) * ETIJ(IPRF12,T2,M1,M2,2,1)
     &             * SMMHCP(IPRF12,NE,T3,3) * ETIJ(IPRF12,T3,N1,N2,3,1)
           END DO
          END DO
          END DO
          END DO

C
C     Transform to the contracted basis
C     =================================


          DO IPRF12 = 1, NPRF12
           IPRF1 = IPR12(1,IPRF12)
           IPRF2 = IPR12(2,IPRF12)

           DO ICTF1 = 1, NCTF1
            DO ICTF2 = 1, NCTF2
             SMMCC(ICTF1,ICTF2,KC1,KC2,KCE) = 
     &               SMMCC(ICTF1,ICTF2,KC1,KC2,KCE)
     &            +  SMMCP(IPRF12,KC1,KC2,KCE)
     &            * CCF1(IPRF1,ICTF1)*CCF2(IPRF2,ICTF2)
            END DO
           END DO
          END DO

          IF ( IPRINT1 .GE. 5 ) THEN
           DO ICTF1 = 1, NCTF1
            DO ICTF2 = 1, NCTF2
             WRITE(LUPRI,'(A,5I5,F12.6)') 
     &          'ICTF1,ICTF2, KCE,KC1,KC2, SMMCC:',ICTF1,ICTF2,
     &          KCE,KC1,KC2,SMMCC(ICTF1,ICTF2,KC1,KC2,KCE)
            END DO
           END DO
          END IF

         END DO
         KCE0 = KCE0 + KCKTE

        END DO
       END DO
      END DO 
  
      IF ( IPRINT1 .GE. 5 ) THEN        
       KCE0 = 1  
       DO JE = 0, JMAXE 
        NHKTE   = JE + 1 
        KCKTE   = NHKTE * ( NHKTE + 1 ) / 2
        KCELAST = KCE0 + KCKTE - 1

        DO KCE = KCE0, KCELAST
  
         DO KC1 = 1, KCKT1
          DO KC2 = 1, KCKT2
           WRITE(LUPRI,'(A,4I5,2F12.6)') 'JE KCE-rel KC1 KC2 ', 
     &           JE, KCE-KCE0+1, KC1, KC2, SMMCC(1,1,KC1,KC2,KCE),
     &            SMMCP(1,KC1,KC2,KCE)

          END DO
         END DO
        END DO

        KCE0 = KCE0 + KCKTE
       END DO
      END IF

C
C     Transform to the spherical basis
C     ================================
C

      NCTF12 = NCTF1 * NCTF2

      KCSQ1  = ICSQADR(J1)
      KCSQ2  = ICSQADR(J2)

      KHE = 1
      KCE = 1
 
      DO JE = 0, JMAXE

       KCSQE = ICSQADR(JE)

       NHKTE = JE + 1
       KHKTE = NHKTE * 2 - 1
       KCKTE = NHKTE * ( NHKTE + 1 ) / 2

       CALL MMCAR2MMSPH(SMMSC(1,1,1,1,KHE),SMMCC(1,1,1,1,KCE),
     &                  CSQ(KCSQ1), CSQ(KCSQ2), CSQE(KCSQE),
     &                  KHKT1,KHKT2,KHKTE,
     &                  KCKT1,KCKT2,KCKTE,NCTF12,IPRINT1 ) 

       KHE = KHE + KHKTE
       KCE = KCE + KCKTE
      END DO

      RETURN
      END 


      SUBROUTINE MMCAR2MMSPH(SMMSC,SMMCC,CSQ1,CSQ2,CSQE,
     &                       KHKT1,KHKT2,KHKTE,
     &                       KCKT1,KCKT2,KCKTE,NCTF12,IPRINT1 )
#include "implicit.h"
#include "priunit.h"

      DIMENSION SMMSC(NCTF12,KHKT1,KHKT2,KHKTE),
     &          SMMCC(NCTF12,KCKT1,KCKT2,KCKTE),
     &          CSQ1(KHKT1,KCKT1), CSQ2(KHKT2,KCKT2),
     &          CSQE(KCKTE,KHKTE)


      IF ( IPRINT1 .GE. 5 ) THEN 
       WRITE(LUPRI,'(A,6I5)') 'C2S: KHKT(1,2,E), KCKT(1,2,E)',
     &        KHKT1,KHKT2,KHKTE,KCKT1,KCKT2,KCKTE

       WRITE(LUPRI,*) 'CSQ1:'
       CALL OUTPUT(CSQ1,1,KCKT1,1,KHKT1,KCKT1,KHKT1,1,LUPRI)

       WRITE(LUPRI,*) 'CSQ2:'
       CALL OUTPUT(CSQ2,1,KCKT2,1,KHKT2,KCKT2,KHKT2,1,LUPRI)

       WRITE(LUPRI,*) 'CSQE:'
       CALL HEADER('CSQE',-1)
       CALL OUTPUT(CSQE,1,KCKTE,1,KHKTE,KCKTE,KHKTE,1,LUPRI)
      END IF     

      DO IS1=1,KHKT1
      DO IS2=1,KHKT2
      DO ISE=1,KHKTE

       DO IC1=1,KCKT1
       DO IC2=1,KCKT2
       DO ICE=1,KCKTE

        DO I = 1, NCTF12
         SMMSC(I,IS1,IS2,ISE) = SMMSC(I,IS1,IS2,ISE)
     &           + SMMCC(I,IC1,IC2,ICE) * CSQ1(IS1,IC1)
     &             * CSQ2(IS2,IC2) * CSQE(ICE,ISE)

        END DO

       END DO
       END DO
       END DO

      END DO
      END DO
      END DO

      RETURN
      END

C     This subroutine resets the array of multipole moments
C     that represents the charge density attributed to the central cell.
      subroutine pbc_reset_cellmom()
      implicit none
#include "pbc.h"
      integer i
      do i = 1,max_qsiz_pbc
         cell_mulmom(i) = 0.0D0
      end do
      return
      end

C     This subroutine adds up electronic multipole moments.
C     It is supposed to be called several times so that the
C     final sum constitute the multipole moments of overlaps
C     between dummy orbitals in the central cell with dummy orbitals
C     in arbitrary cells j.
C
C     Note that DMAT = D^{0j} so be sure to send in the density
C     matrix corresponding to the right j.
      SUBROUTINE pbc_addupcellmom(SMMSC,DMAT,
     &       IAOB1,IAOB2,J1,J2,JMAXE,NCTF1,NCTF2,
     &       ISPCOMP1,ISPCOMP2,nuc_flag,IPRINT1)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
C
#include "pbc.h"

      dimension SMMSC(NCTF1,NCTF2,2*J1+1,2*J2+1,*)
      dimension DMAT(NBASE,NBASE)
      logical nuc_flag

      integer jmindex

      if (nuc_flag) then
         if (J1 .ne. 0 .or. J2 .ne. 0 .or.
     &       NCTF1 .ne. 1 .or. NCTF2 .ne. 1) then
            call quit('pbc_addupcellmom: Programmer error - '//
     &                'invalid arguments passed.')
         end if
         KHE0 = 0
         jmindex = 1
         DO JE=0,JMAXE
            DO M=-JE,JE
C              columns are ordered as 0, 1, -1, 2, -2, etc.
               IF(M.GT.0) THEN
                  KHE = KHE0 + 2*M
               ELSE
                  KHE = KHE0 -2*M+1
               END IF
C              The nuclear charges must already be included in the
C              moments, so here we just add up.
               SPHINT = SMMSC(1,1,1,1,KHE)
C              Add to the cell's moments
               cell_mulmom(jmindex) = cell_mulmom(jmindex) + SPHINT
               jmindex = jmindex + 1
            END DO
            KHE0 = KHE0 + ( 2*JE+1 ) 
         END DO
         return
      end if

      KHKT1 = 2*J1+1
      KHKT2 = 2*J2+1

      DO ICTF1 = 1, NCTF1
       ISTAO1 = INDXBTCH(KNDXBTCH(ISPCOMP1,IAOB1,1)+ICTF1-1,1)

       DO ICTF2 = 1, NCTF2
        ISTAO2 = INDXBTCH(KNDXBTCH(ISPCOMP2,IAOB2,1)+ICTF2-1,1)
        
        DO ISP1=1, KHKT1
         IAO1 = ISTAO1+ISP1-1
         DO ISP2=1, KHKT2
          IAO2 = ISTAO2 + ISP2 - 1

          KHE0 = 0
          jmindex = 1
          DO JE=0,JMAXE
           DO M=-JE,JE
C                     columns are ordered as 0, 1, -1, 2, -2, etc.
             IF(M.GT.0) THEN
                KHE = KHE0 + 2*M
             ELSE
                KHE = KHE0 -2*M+1
             END IF
C            Multiply the purely "geometric" moment of the contr. overlap
C            distributions with (minus) the density matrix to account for
C            the (negative) electric charge.
             SPHINT = -DMAT(IAO1,IAO2)*SMMSC(ICTF1,ICTF2,ISP1,ISP2,KHE)
C            Add to the cell's moments
             cell_mulmom(jmindex) = cell_mulmom(jmindex) + SPHINT
             jmindex = jmindex + 1
           END DO
           KHE0 = KHE0 + ( 2*JE+1 ) 
          END DO

         END DO    ! ISP2
        END DO    ! ISP1
       END DO    ! ICTF2
      END DO    ! ICTF1

      RETURN
      END 

#ifdef MOD_PBC
C     This subroutine is based on PRINTMM, with modifications to
C     adapt it to the needs of the PBC code. / ErikT
      SUBROUTINE pbc_add_to_JFF(SMMSC,
     &       IAOB1,IAOB2,J1,J2,JMAXE,NCTF1,NCTF2,
     &       ISPCOMP1,ISPCOMP2,IPRINT1,DJFARFIELD)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "inforb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
C
#include "pbc.h"
      DIMENSION SMMSC(NCTF1,NCTF2,2*J1+1,2*J2+1,*),
     &          DJFARFIELD(NBAST,NBAST)

      KHKT1 = 2*J1+1
      KHKT2 = 2*J2+1

      DO ICTF1 = 1, NCTF1
       ISTAO1 = INDXBTCH(KNDXBTCH(ISPCOMP1,IAOB1,1)+ICTF1-1,1)

       DO ICTF2 = 1, NCTF2
        ISTAO2 = INDXBTCH(KNDXBTCH(ISPCOMP2,IAOB2,1)+ICTF2-1,1)
        
        DO ISP1=1, KHKT1
         IAO1 = ISTAO1+ISP1-1
         DO ISP2=1, KHKT2
          IAO2 = ISTAO2 + ISP2 - 1
          KHE0 = 0
          jmindex = 1
          DO JE=0,JMAXE
           DO M=-JE,JE
C                     columns are ordered as 0, 1, -1, 2, -2, etc.
             IF(M.GT.0) THEN
                KHE = KHE0 + 2*M
             ELSE
                KHE = KHE0 -2*M+1
             END IF
C            The minus sign accounts for the electronic charge.
             ovl_mulmom(jmindex) = -SMMSC(ICTF1,ICTF2,ISP1,ISP2,KHE)
             jmindex = jmindex + 1
           END DO
           KHE0 = KHE0 + ( 2*JE+1 ) 
          END DO

C         We are almost ready to use the multipole moments, but they are
C         not compatible with the local moments due their prefactors.
C         We put in the right prefactors here.
          call pbc_redefine_q(ovl_mulmom, JMAXE)
C         Contract with the far-field local moments to get the
C         far-field contribution to the Fock matrix.          
          call pbc_contract_qLlat(tmp_jff,ovl_mulmom,qlmax_pbc)
          DJFARFIELD(IAO1,IAO2) = DJFARFIELD(IAO1,IAO2) + tmp_jff

Careful / ErikT
C          if (min(IAO1,IAO2) .eq. 2 .and. max(IAO1,IAO2) .eq. 3) then
C             write(LUPRI,'(2(A,I3),A,D,A)') ' J^{FF}(',IAO1,',',IAO2,
C     &            ') = ',tmp_jff,'  (jmindex vs. KHE???)'
C          end if
Careful

         END DO    ! ISP2
        END DO    ! ISP1
       END DO    ! ICTF2
      END DO    ! ICTF1

      RETURN
      END 

C     This subroutine is based on PRINTMM, with modifications to
C     adapt it to the needs of the PBC code. / ErikT
      SUBROUTINE pbc_printmm(SMMSC,
     &       IAOB1,IAOB2,J1,J2,JMAXE,NCTF1,NCTF2,
     &       ISPCOMP1,ISPCOMP2,LUFILE,IPRINT1)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"

      DIMENSION SMMSC(NCTF1,NCTF2,2*J1+1,2*J2+1,*)
C Remember to put in an error check for the completeness of the 
C file (e.g., number of moments written).

      KHKT1 = 2*J1+1
      KHKT2 = 2*J2+1

      DO ICTF1 = 1, NCTF1
       ISTAO1 = INDXBTCH(KNDXBTCH(ISPCOMP1,IAOB1,1)+ICTF1-1,1)

       DO ICTF2 = 1, NCTF2
        ISTAO2 = INDXBTCH(KNDXBTCH(ISPCOMP2,IAOB2,1)+ICTF2-1,1)
        
        DO ISP1=1, KHKT1
         IAO1 = ISTAO1+ISP1-1
         DO ISP2=1, KHKT2
          IAO2 = ISTAO2 + ISP2 - 1
          WRITE (LUFILE) IAO1,IAO2,JMAXE

          KHE0 = 0
          DO JE=0,JMAXE
           DO M=-JE,JE
C                     columns are ordered as 0, 1, -1, 2, -2, etc.
             IF(M.GT.0) THEN
                KHE = KHE0 + 2*M
             ELSE
                KHE = KHE0 -2*M+1
             END IF

             SPHINT = SMMSC(ICTF1,ICTF2,ISP1,ISP2,KHE)
C            The minus sign accounts for the electronic charge.
             WRITE (LUFILE) -SPHINT
           END DO
           KHE0 = KHE0 + ( 2*JE+1 ) 
          END DO
         END DO    ! ISP2
        END DO    ! ISP1
       END DO    ! ICTF2
      END DO    ! ICTF1

      RETURN
      END 
#endif /* MOD_PBC */

      SUBROUTINE PRINTMM(SMMSC,DMAT,NDMAT,DAUX,RODP,INDDSTR,
     &       IAOB1,IAOB2,J1,J2,JMAX1,JMAX2,JMAXE,NCTF1,NCTF2,
     &       ISPCOMP1,ISPCOMP2,IBTCH,NDISTR,NMOM,LUINTM,IPRINT1)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "fck3aob.h"
#include "fck3odb.h"
#include "fck3thr.h"
#ifdef MOD_DNF
#include "densfit.h"
#endif

      DIMENSION SMMSC(NCTF1,NCTF2,2*J1+1,2*J2+1,*),
     &          INDDSTR(NCTF1,NCTF2,2*J1+1,2*J2+1),
     &          DMAT(NBASE,NBASE),RODP(5),DAUX(NBASEAUX)
      LOGICAL AUXDENS1, AUXDENS2

      ERFAC = ERFCIV(THRNCL)

      KHKT1 = 2*J1+1
      KHKT2 = 2*J2+1

      CALL IZERO( INDDSTR, NCTF1*NCTF2*KHKT1*KHKT2 )

      AUXDENS1 = ((KCLSBTCH(IAOB1).EQ.3) .AND. KCLSBTCH(IAOB2).EQ.0)
      AUXDENS2 = ((KCLSBTCH(IAOB1).EQ.0) .AND. KCLSBTCH(IAOB2).EQ.3)
      IOFF1    = 0
      IOFF2    = 0

      KHE0 = 0
      DO JE=0,JMAXE
       DO M=-JE,JE
C                 columns are ordered as 0, 1, -1, 2, -2, etc.
        IF(M.GT.0) THEN
         KHE = KHE0 + 2*M
        ELSE
         KHE = KHE0 -2*M+1
        END IF

        DO ICTF1 = 1, NCTF1
         ISTAO1 = INDXBTCH(KNDXBTCH(ISPCOMP1,IAOB1,1)+ICTF1-1,1)

         DO ICTF2 = 1, NCTF2
          ISTAO2 = INDXBTCH(KNDXBTCH(ISPCOMP2,IAOB2,1)+ICTF2-1,1)
          DO ISP1=1, KHKT1
           IAO1 = ISTAO1+ISP1-1
           DO ISP2=1, KHKT2
            IAO2 = ISTAO2 + ISP2 - 1

            IF ( ( IAO1 .LE. IAO2 ) .OR. ( ( IAOB1 .NE. IAOB2 ) .AND. 
     &           (IAO2 .LE. IAO1) .AND. ( ISTAO2 .LE. ISTAO1 ) ) ) THEN

             SPHINT = SMMSC(ICTF1,ICTF2,ISP1,ISP2,KHE)

             IF ( NDMAT .GT. 0 ) THEN 
              IF ( IAO1 .NE. IAO2 ) THEN 
                IF (AUXDENS1) THEN
                  IF (LHSFIT) THEN 
                    DENS12 = 0.d0
                    IOFF1  = NBASE
                  ELSE
                    DEN12 = DAUX(IAO1-NBASE)
                  ENDIF
                ELSEIF (AUXDENS2) THEN
                  IF (LHSFIT) THEN 
                    DENS12 = 0.d0
                    IOFF2  = NBASE
                  ELSE
                    DEN12 = DAUX(IAO2-NBASE)
                  ENDIF
                ELSE 
                  DEN12 = DMAT(IAO1,IAO2) + DMAT( IAO2, IAO1 )
                ENDIF
              ELSE 
               DEN12 = DMAT(IAO1,IAO2) 
               IF ( IAOB1 .NE. IAOB2 ) THEN
                 SPHINT = SPHINT * 2 
               END IF
              END IF
             ELSE 
              DEN12 = 0.0D0 
             END IF

             IF ( ABS( SPHINT ) .GT. THRNCL ) THEN

               NMOM = NMOM + 1

              IF ( INDDSTR( ICTF1, ICTF2, ISP1,ISP2) .EQ. 0 )
     &         THEN 
               NDISTR = NDISTR + 1 
               INDDSTR( ICTF1, ICTF2, ISP1,ISP2) = NDISTR
              END IF

              WRITE (LUINTM) JE,M, IAO1-IOFF1,IAO2-IOFF2, J1, J2, IBTCH,
     &                  INDDSTR( ICTF1, ICTF2, ISP1,ISP2),
     &                  RODP(4),RODP(1),RODP(2),RODP(3),
     &                  -SPHINT, DEN12
c              WRITE (LUPRI,'(A,8I4,A,4F10.6,E12.6,F10.6)') 
c     &        'L M IA IB LA LB ICNT I', JE,M, IAO1,IAO2, J1,J2,
c     &            IBTCH, INDDSTR( ICTF1, ICTF2, ISP1,ISP2),
c     &           ' P(M,X,Y,Z) S D', RODP(4),
c     &            RODP(1),RODP(2),RODP(3),-SPHINT, DEN12
             END IF
            END IF
           END DO
          END DO
         END DO
        END DO

       END DO

       KHE0 = KHE0 + ( 2*JE+1 ) 
      END DO

      RETURN
      END 


      SUBROUTINE SETCSQE( ICSQADR, CSQE, LMAX,
     &                    WORK,LWORK,IPRINT1)
#include "implicit.h"
#include "priunit.h"
      DIMENSION ICSQADR(0:LMAX),CSQE(*),WORK(*)

      KOFF = 1 
      DO L=0,LMAX
       ICSQADR(L) = KOFF
       KHKT = 2*L + 1
       KCKT = (L+1)*(L+2)/2

       CALL SPHCOM(L,CSQE(KOFF),KHKT,KCKT,0,0,
     &             WORK,LWORK,IPRINT1)
       KOFF = KOFF + KHKT*KCKT
      END DO

      RETURN
      END 
#else
C to keep intel compiler happy if FMM is not included
      SUBROUTINE DUMFCK3MOM
      RETURN
      END
#endif /* MOD_FMM */
