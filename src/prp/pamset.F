!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

! define task symbols for CALL DIRAC_PARCTL( task )

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck getfck */
      SUBROUTINE GETFCK(FMO,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Get Fock-matrix in MO-basis; use coefficients for MO-transformation
C
C     Written by T.Saue, Sep 16 1996
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dcbham_mod, only: SPINFR
!      use dcbbas_mod, only: N2BBASXQ, NTBAS
!      use dcborb_mod, only: N2ORBXQ, NASHT, NCMOTQ, NORBT
!      use dgroup_mod, only: IPQTOQ, NZ
#include "implicit.h"
C
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      LOGICAL TOBE
      DIMENSION FMO(*),WORK(*)
C
      CALL QENTER('GETFCK')
      KFRSAV = KFREE
C
C     First check if Fock-matrix is on file PAMFCK
C
      LUBUF = 22
      INQUIRE(FILE='PAMFCK',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'PAMFCK','OLD','GETFC1')
C
C       Read FCmo and FVmo
C
        CALL READT(LUBUF,N2ORBXQ,FMO)
        IF (NASHT .GT. 0)
     &     CALL READT(LUBUF,N2ORBXQ,FMO(1+N2ORBXQ))
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
C
         IF (NASHT .EQ. 0) THEN
            CALL MEMGET2('REAL','FAO',KFAO,N2BBASXQ,WORK,KFREE,LFREE)
         ELSE
            CALL MEMGET2('REAL','FAO',KFAO,2*N2BBASXQ,WORK,KFREE,LFREE)
         END IF
         CALL MEMGET2('REAL','CMO',KCMO,NCMOTQ,WORK,KFREE,LFREE)
         IF(SPINFR) THEN
           CALL MEMGET2('INTE','IBEIG',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
         ELSE
           KBEIG = KFREE
         ENDIF
C
         CALL GETFC1(FMO,WORK(KFAO),WORK(KCMO),WORK(KBEIG),
     &             IPRINT,WORK,KFREE,LFREE)
C
         CALL OPNFIL(LUBUF,'PAMFCK','NEW','GETFC1')
         CALL WRITT(LUBUF,N2ORBXQ,FMO)
         IF (NASHT .GT. 0) 
     &      CALL WRITT(LUBUF,N2ORBXQ,FMO(1+N2ORBXQ))
         CLOSE(LUBUF,STATUS='KEEP')
         CALL MEMREL('GETFCK',WORK,1,KFRSAV,KFREE,LFREE)
      ENDIF
C
C     Print section
C
      IF(IPRINT.GE.5) THEN
          CALL HEADER('GETFCK:Fock (FCmo) matrix in MO-basis',-1)
          CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,
     &                NZ,IPQTOQ(1,0),LUPRI)
          IF (NASHT .GT. 0) THEN
             CALL HEADER('GETFCK:Fock (FVmo) matrix in MO-basis',-1)
             CALL PRQMAT(FMO(1+N2ORBXQ),NORBT,NORBT,NORBT,NORBT,
     &                   NZ,IPQTOQ(1,0),LUPRI)
          END IF
      ENDIF
C
      CALL QEXIT('GETFCK')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GETFC1 */
      SUBROUTINE GETFC1(FMO,FAO,CMO,IBEIG,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C     
C     Get Fock-matrix in MO-basis
C
C     Written by T.Saue, Sep 16 1996
C
C***********************************************************************
!      use dcbgen_mod, only: PARCAL
!      use dcbdhf_mod, only: ESOLVE, ESOLVN, IFCKOP, IHRMOP, INTDEF, INTFLG, ISYMOP, NFMAT
!      use dcbham_mod, only: IPRSOL, SOLVEN, SPINFR
!      use dgroup_mod, only: IPQTOQ, NFSYM, NZ
!      use dcbbas_mod, only: I2BASX, N2BBASX, N2BBASXQ, NFBAS, NTBAS
!      use dcborb_mod, only: I2ORBX, ICMOQ, IORB, N2ORBXQ, NASHT, NORB, NORBT
       use dirac_cfg
       use dft_cfg
       use xcint_main
       use fde_mod
       use num_grid_gen
       use num_grid_cfg
#ifdef HAS_PCMSOLVER
      use pcm_linear_response, only: compute_mep_asc, pcm_oper_ao_driver
#endif
#ifdef HAS_PELIB
      use pelib, only: pelib_master
      use pelib_dirac_interface, only: use_pelib, pelib_ifc_fock
#endif

      use priunit_mod, only: LUPRI
      use dirac_partask_mod, only: FDE_PAR, XCINT_PAR
#include "implicit.h"
      INTEGER :: I
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0, DP5 = 0.50D00,D2=2.0D0)
C
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      type(fde_import) :: itmp
CTOBS: Note the use of dcbdhf ! Should be generalized !!!
      LOGICAL TOBE
      CHARACTER*6 F2NAME
      DIMENSION FMO(*),FAO(*),CMO(*),IBEIG(*),WORK(*)
      real(8), allocatable :: xc_matrix(:, :, :)
#ifdef HAS_PCMSOLVER
      real(8), allocatable :: fock_pcm(:, :, :)
#endif

C
      CALL QENTER('GETFC1')
      LUBUF = 22
C
C     1. Get coefficients
C     ===================
C
      IF(SPINFR) THEN
        IOPT = 10
      ELSE
        IOPT = 2
      ENDIF
      IF (IPRINT.GE.10) IOPT = IOPT + 1
      CALL REACMO(LUBUF,'DFCOEF',CMO,DUM,IBEIG,TOTERG,IOPT)
      NFMAT = 1
      IF (NASHT .GT. 0) NFMAT = 2
C
C     1. Get two-electron Fock matrix in AO-basis
C     ===========================================
C
      INQUIRE(FILE='DFFCK2',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'DFFCK2','OLD','GETFC1')
        WRITE(LUPRI,'(/A)')
     &    '* GETFCK: reading old two-electron Fock matrix'
        IF (NASHT .EQ. 0) THEN
           CALL REAFCK(LUBUF,FMO,.TRUE.,1)
         ELSE
           CALL REAFCK(LUBUF,FMO,.TRUE.,2)
         END IF
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
        INTFLG = INTDEF
Chjaaj MAERKE 01-Feb-2001: this is the best we can do now; 
C      this assumes the user uses the same input for restart
C      with *SCF specified as in the .DHF calculation such that
C      any user-specified change of INTFLG (e.g. no SSSS) is used.
C >>>  THIS MUST BE CHANGED FOR MCSCF !!!!
C
        ILLINT = IAND(INTFLG,1)
        ISLINT = IAND(INTFLG/2,1)
        ISSINT = IAND(INTFLG/4,1)
        IGTINT = IAND(INTFLG/8,1)
        WRITE(LUPRI,'(/A,4I2)')
     &    '*** WARNING *** GETFCK: No 2-Fock found.'//
     &    ' Regenerating with .INTFLG =',ILLINT,ISLINT,ISSINT,IGTINT

C
C       Closed shell Hartree - Fock or open shell Hartree-Fock
C       ======================================================
CTROND: This must be generalized !!!
C       Consider it done...
        DO IMAT = 1,NFMAT
C          Totally symmetric operator
           ISYMOP(IMAT) = 1
C          Fock matrix type
           IFCKOP(IMAT) = 1
C          Hermitian operator
           IHRMOP(IMAT) = 1
        END DO
        INQUIRE(FILE='DFDENS',EXIST=TOBE)
        IF(TOBE) THEN
          CALL OPNFIL(LUBUF,'DFDENS','OLD','GETFC1')
          CALL READNS(LUBUF,FAO)
          CLOSE(LUBUF,STATUS='KEEP')
        ELSE
          CALL GENDEN(FAO,CMO,1,IPRINT)
          IF (NASHT .GT. 0) 
     &       CALL GENDEN(FAO(1+N2BBASXQ),CMO,2,IPRINT)
        ENDIF
        CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,FMO,FAO,NFMAT,
     &       INTFLG,IPRINT,WORK(KFREE),LFREE)

      ENDIF


!       construct xc contribution
!       =========================

        if (dirac_cfg_dft_calculation) then
           call genden(fao, cmo, 1, iprint)
           if (dft_cfg_overlap_diagnostic) then    
              call reset_num_grid()
              call stash_zipgrid(num_grid_cfg_zipgrid)       
           end if
           call generate_num_grid(fao)
           if (dft_cfg_overlap_diagnostic) then    
              call unstash_zipgrid(num_grid_cfg_zipgrid)       
           end if
#ifdef VAR_MPI
            if (parcal) call dirac_parctl( XCINT_PAR )
#endif
           allocate(xc_matrix(ntbas(0), ntbas(0), nz))
           xc_matrix = 0.0d0
           call integrate_xc(xc_mat_dim               = ntbas(0),
     &                        xc_nz                    = nz,
     &                        xc_dmat_0                = fao,
     &                        xc_nr_dmat               = 0,
     &                        xc_nr_fmat               = 1,
     &                        xc_fmat                  = xc_matrix,
     &                        xc_do_potential          = .true.,
     &                        xc_do_overlap_diagnostic = 
     &                           dft_cfg_overlap_diagnostic)
           call daxpy(n2bbasxq, 1.0d0, xc_matrix, 1, fmo, 1)
           deallocate(xc_matrix)
        end if

        if (dirac_cfg_fde) then
           call fde_get_import_info(itmp)
!
! here we should not skip the following call, even if we don't
! want the fde contributions to response, since it will be
! necessary to start by adding the contributions from the potential 
! to the fmat
!
           if (itmp%im_frozen.and.itmp%im_update_vemb) then
              if (.not.dirac_cfg_dft_calculation) then
                 call genden(fao, cmo, 1, iprint)
              endif

!             call fde_set_nz(nz)
#ifdef VAR_MPI
            if (parcal) call dirac_parctl( FDE_PAR )
#endif
              call fde_calculate_emb_pot_mat(ntbas(0),fao,fmo)

              write (*,*) 'FDE Fock contributions (re)calculated'
           else if (itmp%im_vemb) then
              write (*,*) 'FDE Fock contributions (CHECKPOINT) added'
           endif
        endif
C
C     3a. Solvent: get solvent contribution
C     =====================================
C
      IF(SOLVEN) THEN
C        CALL DENSTY(FAO,CMO,IPRINT)        
CJTH check this
        NF = 1
        CALL SOLFCK(FMO,FAO,NF,ESOLVE,ESOLVN,WORK(KFREE),LFREE,IPRSOL)
      ENDIF          

#ifdef HAS_PCMSOLVER
        if (dirac_cfg_pcm) then
           call genden(fao, cmo, 1, iprint)
! Generate the total MEP and ASC, we are doing response
! so pcm_energy_driver is not called...
! fao is the AO density matrix
           call compute_mep_asc(fao, work(kfree), lfree)
           write(lupri,'(a)') '* Calling pcm_oper_ao_driver to get '//
     &                    'J + X(D) contribution.'
           allocate(fock_pcm(ntbas(0), ntbas(0), nz))
           fock_pcm = 0.0d0
! fmo is the Fock matrix in AO basis, fao is the density matrix
           call pcm_oper_ao_driver(fock_pcm, 'TotASC'//char(0), 
     &          work(kfree), lfree)
           call daxpy(n2bbasx, -1.0d0, fock_pcm, 1, fmo, 1)
!     print *, "FOCK_PCM"
!     call prqmat(fock_pcm, ntbas(0), ntbas(0), ntbas(0), ntbas(0), &
!                      nz, ipqtoq(1,0), 6)
           deallocate(fock_pcm)
        endif 
#endif

#ifdef HAS_PELIB
! edh: As the solvent methods before us (see above) we need
! to construct the (static) PE contributions from environment
! before starting a response calculation
      IF (use_pelib()) THEN
        call genden(fao, cmo, 1, iprint)
        call pelib_ifc_fock(fao, fmo, tot_nrg, el_nrg)
      ENDIF
#endif 

C
C     3. Get one-electron Fock matrix
C     ===============================
C
      INQUIRE(FILE='DFFCK1',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'DFFCK1','OLD','GETFC1')
        WRITE(LUPRI,'(/A)')
     &    '* GETFCK: reading old one-electron Fock matrix'
        CALL REAFCK(LUBUF,FAO,.TRUE.,1)
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
        WRITE(LUPRI,'(/A)')
     &    '*** WARNING *** GETFCK: No 1-Fock found. Regenerating.'
        CALL ONEFCK(FAO,IPRINT,WORK(KFREE),LFREE)
      ENDIF
C
C     4. Add 1- and 2-electron matrix and transform to MO-basis
C     =========================================================
C
      CALL DAXPY(N2BBASXQ,D1,FMO,1,FAO,1) 
      IF (NASHT .GT. 0) THEN
         CALL DCOPY(N2BBASXQ,FMO(1+N2BBASXQ),1,FAO(1+N2BBASXQ),1)
C
C        Scale active Fock matrix with 0.5 (all Fock matrices are
C        (multiplied with 2 in TWOFCK)
C
         CALL DSCAL(N2BBASXQ,DP5,FAO(1+N2BBASXQ),1)
      END IF
      IOFF_MO = 0
      IOFF_AO = 0
      DO IMAT = 1,NFMAT
        CALL DZERO(FMO(1+IOFF_MO),N2ORBXQ)
        DO 10 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 10
        IF(SPINFR) THEN
          CALL QBTRANS(0,'AOMO','S',D0,
     &      NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &      FAO(I2BASX(I,I)+1+IOFF_AO),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &      FMO(I2ORBX(I,I)+1+IOFF_MO),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &      CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &      IBEIG(IORB(I)+1),
     &      CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &      IBEIG(IORB(I)+1),
     &      WORK(KFREE),LFREE,IPRINT)
        ELSE
          CALL QTRANS('AOMO','S',D0,
     &      NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &      FAO(I2BASX(I,I)+1+IOFF_AO),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &      FMO(I2ORBX(I,I)+1+IOFF_MO),NORBT,NORBT,NZ,IPQTOQ(1,0),
     &      CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &      CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &      WORK(KFREE),LFREE,IPRINT)
        ENDIF
 10     CONTINUE
        IOFF_MO = IOFF_MO + N2ORBXQ
        IOFF_AO = IOFF_AO + N2BBASXQ
      END DO
      IF(IPRINT.GE.5) THEN
        CALL HEADER('GETFC1:Fock (FCmo) matrix in MO-basis',-1)
        IF (SPINFR) THEN 
         WRITE(LUPRI,'(4X,A/)') 
     &   'Applied boson symmetry based AO_to_MO transformation !'
        ENDIF                       
        CALL PRQMAT(FMO,NORBT,NORBT,NORBT,NORBT,NZ,
     &    IPQTOQ(1,0),LUPRI)
        IF (NASHT .GT. 0) THEN
           CALL HEADER('GETFC1:Fock (FVmo) matrix in MO-basis',-1)
           CALL PRQMAT(FMO(1+N2ORBXQ),NORBT,NORBT,NORBT,NORBT,NZ,
     &       IPQTOQ(1,0),LUPRI)
        END IF
      ENDIF
C
      CALL QEXIT('GETFC1')
      RETURN
C
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'GETFC1: Coefficient file not found !'
      CALL QUIT('GETFC1: Coefficients not found !')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck REDPAR */
      SUBROUTINE REDPAR(IREDJOB,GPCI,GPOE,GPOP,JXOPE,JXOPP,
     &                  NZXOPE,NZXOPEQ,NZXOPP,NZXOPPQ,NZXOPT,NZXOPTQ,
     &                  NZVAR,NZVARQ,NZCONF,NZCONFQ,THRCOM,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Reduce number of parameters by deleting parameters
C     for which the corresponding gradient element is
C     below threshold THRCOM
C
C     Written by T.Saue - July 22 1997
C
C     Last revision: S. Knecht - Jan 2010
C     added control variable IREDJOB:
C     bit 1 of IREDJOB: reduction in configurational vector
C     bit 2 of IREDJOB: reduction in e-e orbital rotations
C     bit 3 of IREDJOB: reduction in e-p orbital rotations
C                    
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dgroup_mod, only: NZ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I, II, J
      REAL(dp), PARAMETER :: D0 = 0.0D0, DC = 100.0D0
C
#include "dgroup.h"
      DIMENSION GPCI(*),GPOE(*),GPOP(*),WORK(*),
     &          JXOPE(2,NZXOPE),JXOPP(2,NZXOPP)
      LOGICAL LBIT
C
      KFRSAV = KFREE
      IF(IPRINT.GE.4) CALL HEADER('Output from REDPAR',-1)        

C
C     Configurational part
C     ====================
C
      IF(LBIT(IREDJOB,1)) THEN
        WRITE(LUPRI,'(3A)')' * REDPAR: reduction in conf vector ',
     &                     ' not implemented; thus no reduction ',
     &                     ' in conf vector.'
        WRITE(LUPRI,'(1A)')' * Continuing with orbital rotations.'
      ENDIF
C
C     e-e orbital rotations
C     =====================
C
      IF(LBIT(IREDJOB,2)) THEN
        IF(NZXOPE.GT.0) THEN
          CALL MEMGET('REAL',KBUF,NZXOPEQ,WORK,KFREE,LFREE)
          II = 0
          DO I = 1,NZXOPE
            DTST = D0
            IZOFF = 0
            DO IZ = 1,NZ
              DTST = DTST + GPOE(I+IZOFF)*GPOE(I+IZOFF)  
              IZOFF = IZOFF + NZXOPE
            ENDDO
            IF(SQRT(DTST).GE.THRCOM) THEN
              II = II + 1
              JXOPE(1,II) = JXOPE(1,I)
              JXOPE(2,II) = JXOPE(2,I)
              IOFF = KBUF-1
              JOFF = 0
              DO IZ = 1,NZ
                WORK(IOFF+II) = GPOE(JOFF+I)
                IOFF = IOFF + NZXOPE
                JOFF = JOFF + NZXOPE
              ENDDO
            ENDIF
          ENDDO
          IF(II.LT.NZXOPE) THEN
            IOFF = KBUF
            JOFF = 1
            DO IZ = 1,NZ
              CALL DCOPY(II,WORK(IOFF),1,GPOE(JOFF),1)
              IOFF = IOFF + NZXOPE
              JOFF = JOFF + II
            ENDDO
            IF(IPRINT.GE.4) THEN
              DTST = DC*(NZXOPE-II)/NZXOPE
              WRITE(LUPRI,'(A,I10,A,F5.1,A)') 
     &              'e-e rotations: ',II,' (reduction: ',DTST,'%)'
              WRITE(LUPRI,'(2I6)') (JXOPE(1,J),JXOPE(2,J),J=1,II)
              CALL PRQMAT(GPOE,II,1,II,1,NZ,1,LUPRI)
            ENDIF
            NZXOPE  = II
            NZXOPEQ = NZXOPE*NZ
          ENDIF
          CALL MEMREL('REDPAR.ee',WORK,1,KFRSAV,KFREE,LFREE)
        ENDIF
      ENDIF
C
C     e-p orbital rotations
C     =====================
C
      IF(LBIT(IREDJOB,3)) THEN
        IF(NZXOPP.GT.0) THEN
          CALL MEMGET('REAL',KBUF,NZXOPPQ,WORK,KFREE,LFREE)
          II = 0
          DO I = 1,NZXOPP
            DTST = D0
            IZOFF = 0
            DO IZ = 1,NZ
              DTST = DTST + GPOP(I+IZOFF)*GPOP(I+IZOFF)  
              IZOFF = IZOFF + NZXOPP
            ENDDO
            IF(SQRT(DTST).GE.THRCOM) THEN
              II = II + 1
              JXOPP(1,II) = JXOPP(1,I)
              JXOPP(2,II) = JXOPP(2,I)
              IOFF = KBUF-1
              JOFF = 0
              DO IZ = 1,NZ
                WORK(IOFF+II) = GPOP(JOFF+I)
                IOFF = IOFF + NZXOPP
                JOFF = JOFF + NZXOPP
              ENDDO
            ENDIF
          ENDDO
          IF(II.LT.NZXOPP) THEN
            IOFF = KBUF
            JOFF = 1
            DO IZ = 1,NZ
              CALL DCOPY(II,WORK(IOFF),1,GPOP(JOFF),1)
              IOFF = IOFF + NZXOPP
              JOFF = JOFF + II
            ENDDO
            IF(IPRINT.GE.4) THEN
              DTST = DC*(NZXOPP-II)/NZXOPP
              WRITE(LUPRI,'(A,I10,A,F5.1,A)') 
     &              'e-p rotations: ',II,'(reduction: ',DTST,'%)'
              WRITE(LUPRI,'(2I6)') (JXOPP(1,J),JXOPP(2,J),J=1,II)
            ENDIF
            NZXOPP  = II
            NZXOPPQ = NZXOPP*NZ
          ENDIF
          CALL MEMREL('REDPAR.ep',WORK,1,KFRSAV,KFREE,LFREE)
        ENDIF      
      ENDIF      
C
      NZXOPT   = NZXOPE + NZXOPP
      NZXOPTQ  = NZXOPT * NZ
C
C     Total number of parameters
C     ==========================
C
      NZVAR    = NZCONF  + NZXOPT
      NZVARQ   = NZCONFQ + NZXOPTQ
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DEFORB */
      SUBROUTINE DEFORB(IORBCL,NIND,INDSTR,WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     This routine reads orbital strings and classifies
C     and indexes orbitals
C
C     Written by T.Saue Oct 3 1997
C
C***********************************************************************
!      use dcborb_mod, only: NORBT
!      use dcbbas_mod, only: MXFBAS
#include "implicit.h"
C
#include "dcborb.h"
#include "dcbbas.h"
      CHARACTER INDSTR(3,2)*72
      DIMENSION IORBCL(4,NORBT),NIND(2,3,2),WORK(*)
C
      CALL QENTER('DEFORB')
      KFRSAV = KFREE
C
C     Allocate buffer areas
C     =====================
C
      CALL MEMGET2('INTE','NIND',KIND ,MXFBAS,WORK,KFREE,LFREE)
      CALL DEFOR1(IORBCL,NIND,INDSTR,WORK(KIND),
     &            WORK,KFREE,LFREE,IPRINT)
      CALL MEMREL('DEFORB',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DEFORB')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck DEFOR1 */
      SUBROUTINE DEFOR1(IORBCL,NIND,INDSTR,IND,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     This routine reads orbital strings and classifies
C     and indexes orbitals
C
C     Three types are defined:
C     ITYP = 1   inactive
C     ITYP = 2   active
C     ITYP = 3   secondary
C
C     NIND(1,ITYP,IFRP) - number of positronic orbitals 
C     NIND(2,ITYP,IFRP) - number of electronic orbitals
C       ...of type ITYP and fermion ircop IFRP  
C
C     IORB is the absolute index of the orbital
C     IORBCL(1,IORB) - ITYP
C     IORBCL(2,IORB) - index within fermion irrep
C     IORBCL(3,IORB) - fermion irrep IFRP
C     IORBCL(4,IORB) - absolute index for given type ITYP
C                      and electron/positron
C
C     Written by T.Saue Oct 3 1997
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dcborb_mod, only: IORB, NESH, NORBT, NPSH
!      use dcbbas_mod, only: MXFBAS, NFBAS
!      use dgroup_mod, only: NFSYM
#include "implicit.h"
      INTEGER :: I, II, J
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      CHARACTER INDSTR(3,2)*72
      DIMENSION IORBCL(4,NORBT),IND(MXFBAS),
     &          NIND(2,3,2),WORK(KFREE)
C
      NDIM = NORBT*4
      CALL IZERO(IORBCL,NDIM)
C
C     Unpack strings
C     ==============
C
      DO ITYP = 1,3
        KOFE = 0
        KOFP = 0
        DO IFRP = 1,NFSYM
          NVEC = 1
          CALL NUMLST(INDSTR(ITYP,IFRP),IND,
     &                NFBAS(IFRP,0),-NFBAS(IFRP,2),NFBAS(IFRP,1),
     &                IFRP,NVEC)
          CALL ORBCNT(IND,NVEC,NPSH(IFRP),NESH(IFRP),
     &                NIND(1,ITYP,IFRP),NIND(2,ITYP,IFRP))
C
C         Absolute indices of positronic orbitals
C
          IOFF = NPSH(IFRP)+IORB(IFRP)+1
          DO I = 1,NIND(1,ITYP,IFRP)
            J           = IOFF+IND(I)
            IORBCL(1,J) = ITYP
            IORBCL(2,J) = NPSH(IFRP)+IND(I)+1
            IORBCL(3,J) = IFRP
            KOFP = KOFP + 1
            IORBCL(4,J) = KOFP
            
          ENDDO
C
C         Absolute indices of electronic orbitals
C
          IOFF = NPSH(IFRP)+IORB(IFRP)
          JOFF = NIND(1,ITYP,IFRP)
          DO I = 1,NIND(2,ITYP,IFRP)
            II = I + JOFF
            J = IOFF + IND(II) 
            IORBCL(1,J) = ITYP
            IORBCL(2,J) = IND(II)
            IORBCL(3,J) = IFRP
            KOFE = KOFE + 1
            IORBCL(4,J) = KOFE
          ENDDO
        ENDDO
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.3) THEN
        CALL HEADER('Output from DEFORB',-1)
        WRITE(6,'(A)')
     &     '* Defined types:',
     &     '  ITYP = 1   inactive',
     &     '  ITYP = 2   active',
     &     '  ITYP = 3   secondary'
        WRITE(6,'(/3X,A)')
     &      'IFRP IND  ITYP ITYP_IND'
        DO I = 1,NORBT
          IF(IORBCL(1,I).NE.0) THEN
             WRITE(LUPRI,'(4I5)')
     &     IORBCL(3,I),IORBCL(4,I),IORBCL(1,I),IORBCL(2,I)
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Ostring */
      SUBROUTINE OSTRING(INDSTR,SKIPEE,SKIPEP,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Set orbital strings         
C     Stolen by T.Saue Aug 24 2000
C
C      Revisions: Aug 2002 hjaaj
C                 May 2003 MI
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dgroup_mod, only: FREP, NFSYM
!      use dcborb_mod, only: NASH, NESH, NISH, NOCC, NORB, NPSH, NSSH
#include "implicit.h"
      INTEGER :: I, J
C
#include "dgroup.h"
#include "dcborb.h"
C
      CHARACTER*72 INDSTR(3,2)
      LOGICAL SKIPEE,SKIPEP,ORBDEF(3)
      DIMENSION WORK(*),IR(2,2),NR(2)

      CALL QENTER('OSTRING')
C
      DO I = 1,3
         ORBDEF(I) = .FALSE.
         DO J = 1, NFSYM
            ORBDEF(I) = ORBDEF(I) .OR. (INDSTR(I,J) .NE. ' ')
         END DO
      END DO
C
C     Inactive orbitals...
C
      IF (.NOT. ORBDEF(1)) THEN
         DO I = 1,NFSYM
            INDSTR(1,I) = ' '
            IF (NISH(I).EQ.1) THEN
               WRITE(INDSTR(1,I),'(I1)') 1
            ELSEIF (NISH(I).GT.1) THEN
               WRITE(INDSTR(1,I),'(I6,A2,I6)') 1,'..',NISH(I)
            ENDIF
         ENDDO
      END IF
C
C     Active orbitals...
C
      IF (.NOT. ORBDEF(2)) THEN
         DO I = 1,NFSYM
            INDSTR(2,I) = ' '
            IF (NASH(I).EQ.1) THEN
               WRITE(INDSTR(2,I),'(I6)') NISH(I)+1
            ELSEIF (NASH(I).GT.1) THEN
               WRITE(INDSTR(2,I),'(I6,A2,I6)')
     &               NISH(I)+1,'..',NISH(I)+NASH(I)
            ENDIF
         ENDDO
      END IF
C
C     Secondary orbitals...
C
      IF (.NOT. ORBDEF(3)) THEN
        DO I = 1,NFSYM
          INDSTR(3,I) = ' '
          IT=1
          IF(.NOT.SKIPEP .AND. NPSH(I).GT.0) THEN
            IF (NPSH(I) .EQ. 1) THEN
              WRITE(INDSTR(3,I),'(I2)') -1
              IT=3
            ELSE
              WRITE(INDSTR(3,I),'(I6,A4)') -NPSH(I),'..-1'
              IT=11
            ENDIF
          ENDIF
          IF(.NOT.SKIPEE .AND. NSSH(I).GT.0) THEN
            IF (IT .GT. 1) THEN ! e-p rotations were added above
              WRITE(INDSTR(3,I)(IT:IT),'(A1)') ','
              IT=IT+1
            ENDIF
            IF(NSSH(I).EQ.1) THEN
              WRITE(INDSTR(3,I)(IT:IT+5),'(I6)') NOCC(I)+1
            ELSE
              WRITE(INDSTR(3,I)(IT:IT+13),'(I6,A2,I6)')
     &             NOCC(I)+1,'..',NESH(I)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C
      nvec_inactive = 0
      nvec_active   = 0
C
C      Definition of orbital classes
C      -----------------------------
C
      IF(IPRINT.GE.0) THEN
        WRITE(LUPRI,'(/4X,A)') '* Orbital classes:'
      ENDIF
      DO I = 1,NFSYM
        IF(SKIPEP) THEN
          NR(I)   = NESH(I)
          IR(1,I) = 1
          IR(2,I) = NESH(I)
        ELSE
          NR(I)   =  NORB(I)
          IR(1,I) = -NPSH(I)
          IR(2,I) =  NESH(I)
        ENDIF
        IF(IPRINT.GE.0) THEN
           WRITE(LUPRI,'(4X,A,A3)') 'Fermion ircop:',FREP(I)
        ENDIF
C
C       Inactive orbitals
C       -----------------
C
        NVEC = 0
        CALL  NUMLST(INDSTR(1,I),IDUMMY,NR(I),
     &             IR(1,I),IR(2,I),I,NVEC)
        nvec_inactive = nvec_inactive + nvec
        IF(IPRINT.GE.0) THEN
          IF(NVEC.EQ.0) THEN
            WRITE(LUPRI,'(4X,A)') 
     &         '- Inactive orbitals  : No electrons'
          ELSE
            WRITE(LUPRI,'(4X,A,A56)') 
     &         '- Inactive orbitals  : ', INDSTR(1,I)
          ENDIF
        ENDIF
C
C       Active orbitals
C       ---------------
C
        NVEC = 0
        CALL  NUMLST(INDSTR(2,I),IDUMMY,NR(I),
     &             IR(1,I),IR(2,I),I,NVEC)
        nvec_active = nvec_active + nvec
        IF(IPRINT.GE.0) THEN
          IF(NVEC.EQ.0) THEN
            WRITE(LUPRI,'(4X,A)') 
     &         '- Active orbitals    : No electrons'
          ELSE
            WRITE(LUPRI,'(4X,A,A56)') 
     &         '- Active orbitals    : ', INDSTR(2,I)
          ENDIF
        ENDIF
C
C       Secondary orbitals
C       ------------------
C
        NVEC = 0
        CALL  NUMLST(INDSTR(3,I),IDUMMY,NR(I),
     &       IR(1,I),IR(2,I),I,NVEC)
        IF(IPRINT.GE.0) THEN
          IF(NVEC.EQ.0) THEN
            WRITE(LUPRI,'(4X,A)') 
     &         '- Secondary orbitals : No electrons'
          ELSE
            WRITE(LUPRI,'(4X,A,A56)') 
     &         '- Secondary orbitals : ', INDSTR(3,I)
          ENDIF
        ENDIF 
      ENDDO

      if ((nvec_inactive + nvec_active) == 0) then 
         print *, 'error: nvec_inactive + nvec_active is zero'
         print *, 'please specify the orbital occupation explicitly'
         call quit('please specify occupation explicitly')
      end if

      CALL QEXIT('OSTRING')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck xrspar */
      SUBROUTINE XRSPAR(XRS_INDSTR,XRS_SKIPEE,XRS_SKIPEP,IPRINT)
C***********************************************************************
C
C     Initialize configurational and orbital parameters for
C     linear response calculation
C
C     Written by T. Saue July 18 1997
C
C***********************************************************************

      use priunit_mod, only: LUPRI
!      use dgroup_mod, only: NZ
!      use dcborb_mod, only: NORBT
!      use dcbxrs_mod, only: INDSTR, JOPSY, NZCONF, NZCONFQ, NZHOPE, NZVAR, NZVARQ, NZXOPE, NZXOPEQ, NZXOPP, NZXOPPQ, NZXOPT, NZXOPTQ, SKIPEE, SKIPEP
      use orbital_rotation_indices, only: setxop
      use conf_parameters, only: allocate_conf_array

      implicit none

      integer, intent(in) :: iprint
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbxrs.h"
C
      LOGICAL GASWF, XRS_SKIPEE, XRS_SKIPEP
      CHARACTER*72 XRS_INDSTR(3,2)
      integer, allocatable :: norot(:)
C
      CALL QENTER('XRSPAR')
C
      SKIPEE = XRS_SKIPEE
      SKIPEP = XRS_SKIPEP
      INDSTR = XRS_INDSTR
C     
C     Configurational parameters
C     ==========================
C
      NZCONF   = 0
      NZCONFQ  = NZCONF *NZ
C     ... I guess it should be *MZ where MZ = MIN(2,NZ) ??? /hjaaj mar 2003

      call allocate_conf_array(nzconfq)


C     Orbital parameters
C     ==================
C
C===TODO start
C     hjaaj June 2001: NOROT should be 
C     taken from common block dcbnrt.h
C     but is not defined for DHF/DFT yet (only for KRMC now)
C     hjaaj June 2001: GASWF should be based on either
C     MCSCF common block or that this WF is not KRMC
C
      allocate(norot(norbt))
      norot = 0
      GASWF = .FALSE.
      CALL SETXOP(INDSTR,XRS_SKIPEE,XRS_SKIPEP,GASWF,JOPSY,NOROT,
     &            NZXOPE,NZHOPE,NZXOPP,
     &            IPRINT)
      deallocate(norot)
C===TODO end


      NZXOPEQ  = NZXOPE *NZ
      NZXOPPQ  = NZXOPP *NZ
      NZXOPT   = NZXOPE + NZXOPP
      NZXOPTQ  = NZXOPT *NZ
C
C     Total number of parameters
C     ==========================
C
      NZVAR    = NZCONF + NZXOPT
      NZVARQ   = NZCONFQ + NZXOPTQ
C
      IF (IPRINT .GE. 3) THEN
         WRITE (LUPRI,'(/A,3I10)')
     &   'No. of configuration, e-e, and e-p variables:',
     &   NZCONF,NZXOPE,NZXOPP
      END IF
C
      CALL QEXIT('XRSPAR')
      END

      SUBROUTINE set_prop_flags_ccsd()
C
C---------------Description--------------------------------------------
C
C   set flags to perform property computations if lambda is activated
C
C   Johann Pototschnig, Fall 2019

C---------------Common Blocks--------------------------------------
C
      use mxcent_mod, only: MXCENT, MXCOOR
!      use dcbprp_mod, only: DOEXP, HDF5IO, MXPRP_WF, NPRP_WF, PRP_WF
#include "implicit.h"
      INTEGER :: I
!#include "dgroup.h"
#include "dcbprp.h"

      character(4)  :: WF
C
      WF='CCSD'
C     Switch on evaluation of property output
      DOEXP = .TRUE.
C     disable data transfer via the hdf5 file(s), until
C     support is stable for all functionality related to ccsd, such as
C     mp2fnos
      HDF5IO = .FALSE.
C     Add this density matrix to the list to be evaluated
      IPRP_WF = 0
      DO I = 1, NPRP_WF
         IF (PRP_WF(I).EQ.WF) IPRP_WF = I
      END DO
      IF (IPRP_WF.EQ.0) THEN
         IF (NPRP_WF .GE. MXPRP_WF) THEN
             IPRP_WF = 1 ! Unlikely that this happens, just skip the HF print out in this case
         ELSE
             NPRP_WF = NPRP_WF + 1
             IPRP_WF = NPRP_WF
         END IF
         PRP_WF(IPRP_WF) = WF
      END IF

      END
