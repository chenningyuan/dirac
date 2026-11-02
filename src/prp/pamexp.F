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

      SUBROUTINE KUPRPEXP(WF,EXPVAL,WORK,LWORK,state)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIEXP/ at the KU level
C
C     Written by Chima Chibueze Februari 2022
C     Last revision April 2022 - cchibueze
C
C***********************************************************************
#include "implicit.h"
C
      DIMENSION EXPVAL(*),WORK(LWORK)
      CHARACTER WF*4
      INTEGER, OPTIONAL :: state
      INTEGER :: state_local
C
C
      CALL QENTER('PRPEXP')
#include "memint.h"
      if (present(state)) then
          state_local = state
      else
          state_local = 1
      end if
      CALL KUPRPEX1(WF,EXPVAL,state)
      CALL QEXIT('PRPEXP')
C
      RETURN
      END
      
      SUBROUTINE KUPRPEX1(WF,EXPVAL,state)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Chima Chibueze Februari 2022
C     Last revision April 2022 - cchibueze
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use mxcent_mod, only: MXCENT
!      use dcborb_mod, only: NESH
!      use dcbexp_mod, only: LEXPP, NEXPP
!      use dcbxpr_mod, only: IPRPSYM, IPRPTIM, PRPNAM
!      use dcbprl_mod, only: IPRLTYP
!      use dcbham_mod, only: NONREL, X2C
!      use dcbgen_mod, only: KCHARG, LUX2C
      use labeled_storage
      use checkpoint
      use dircmo
!     use x2c_fio
#include "implicit.h"
#include "mxcent.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,DTOL=1.0D-16,D2=2.0D0) 
C
#include "dcborb.h"
#include "dcbexp.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dcbham.h"
#include "dcbgen.h"
      DIMENSION EXPVAL(*)
      CHARACTER MXFORM*6,FMT*6
      CHARACTER WF*4
      INTEGER :: state

      type(file_info_t), save :: exacc_file
      type(file_info_t), save ::  x2cmat_file, reex_file
      logical :: tobe_dirac, tobe_respect, prp_hdf5=.false.
      logical :: fndlab
      integer :: i,j,p,q, indxpr,pq
      integer :: nmo, nao, nel, nocc_ccsd, ninact, chargenuc
      complex(8) :: expvali
      integer , dimension(1) :: nocc_cc
      integer, allocatable :: mo_occ(:)
      real(8), allocatable :: rcoeffvec(:), icoeffvec(:)
      
      integer              :: prp_fer_rep, prp_mat_dim
      logical              :: prp_triang
      logical, dimension(4):: quat_array_filled = .false.
      real(8), allocatable :: scr_mat1(:,:)
      real(8), allocatable :: scr_mat2(:)

      real(8), allocatable :: rkudmvec(:), ikudmvec(:)
      real(8), allocatable :: rprpvec(:), iprpvec(:)
      real(8), allocatable :: qx2cmatvec(:), qcoeffvec(:)
      real(8), allocatable :: aopropervec(:)
      real(8), allocatable :: qqx2cmatvec(:,:) 
      real(8), allocatable :: qx2cmat(:,:,:), qcoeffmat(:,:,:)
      complex(8), allocatable :: ckudmvec(:), ccoeffvec(:), cprpvec(:)
      complex(8), allocatable ::kudm_mo(:,:),kudm_tmp(:,:),kudm_ao(:,:)
      complex(8), allocatable :: cprpmat(:,:), ccoeffmat(:,:)
      character*(:),allocatable :: dataid
      character(len=4)  :: numstr
      character(len=12) :: prplab
      character(len=6)  :: reex(3)
      complex(8), parameter :: DC1 = cmplx(1.0d0,0.0d0)
      complex(8), parameter :: DC0 = cmplx(0.0d0,0.0d0)
      complex(8), parameter :: DCI = cmplx(0.0d0,1.0d0)
      complex(8), parameter :: DCM1 = cmplx(-1.0d0,0.0d0)
C
C
C     Determining whether DIRAC or ReSpect spinors are used
C     =====================================================
C
      tobe_dirac = .false.
      tobe_respect = .false. 
      inquire (file='RSD_MOS',exist=tobe_respect)
      if (tobe_respect) then
          write(*,*) "   Found ReSpect file: using it"
      else
          tobe_dirac = .true.
          write(*,*) "   Using DIRAC MOs"
      endif
C
C
C     Retrieving dimensionality info
C     ==============================
C

      if (tobe_dirac) then
          dataid='/input/aobasis/1/n_ao'
          nmo = nesh(1)
          call checkpoint_read(dataid,idata=nao)
      elseif (tobe_respect) then

          ! checking existence of reex density file
          if (reex_file%status == -1) then
              reex_file%type = 2
              reex_file%name = "REEX.h5"
              reex_file%status = 0
          endif
          call lab_read(reex_file,'nmo' ,nmo)
          call lab_read(reex_file,'nsph',nao)
          nao = int(nao / 2)
          nmo = int(nmo / 2)
      end if

      !counting number of electrons
      nel = 0
      chargenuc = 0
      call rmolchr(chargenuc)
      nel = chargenuc - kcharg
C
C
C     Retrieving MO coefficient matrix
C     ================================
C
      if (tobe_dirac) then

          ! allocating MO coefficent arrays
          allocate(qcoeffvec(nao * nmo * 4))
          allocate(qcoeffmat(nao, nmo, 4  ))
          allocate(ccoeffmat(nao*2, nmo*2 ))

          ! reading in coefficient matrix from file from file
          qcoeffvec = D0  
          call reacmo_new(cmo=qcoeffvec)
          qcoeffmat = D0
          qcoeffmat = reshape(qcoeffvec, shape=(/nao,nmo,4/),
     &                order=(/1,2,3/))
          deallocate(qcoeffvec)
          ccoeffmat = DC0

          ! converting from quaternion to complex format
          call qtoc(qcoeffmat,ccoeffmat,1,nao,nmo,nao,nmo,4)
          deallocate(qcoeffmat) 

      elseif (tobe_respect) then

          ! allocating coefficient arrays
          allocate(rcoeffvec(4*nao*nmo  ))
          allocate(icoeffvec(4*nao*nmo  ))
          allocate(ccoeffvec(4*nao*nmo  ))
          allocate(ccoeffmat(2*nao,2*nmo))

          ! reading in MO coefficients from file
          rcoeffvec = D0
          icoeffvec = D0
          call lab_read (reex_file,'sphMO _RE',rdata=rcoeffvec)
          call lab_read (reex_file,'sphMO _IM',rdata=icoeffvec)
          ccoeffvec = DC0
          ccoeffvec = cmplx(rcoeffvec,icoeffvec,kind=kind(1.0d0))
          deallocate(rcoeffvec)
          deallocate(icoeffvec)
          ccoeffmat = DC0
          ccoeffmat = reshape(ccoeffvec,shape=(/2*nao,2*nmo/),
     &                order=(/1,2/))
          deallocate(ccoeffvec)

      endif ! (tobe...)
C
C
C     Retrieving DHF or CCSD density
C     ==============================
C
      nmo = nmo * 2 !nao stays the same; alpha basis = beta basis

      ! allocating density arrays
      allocate(kudm_mo(nmo,nmo))
      allocate(kudm_tmp(nao*2,nmo))
      allocate(kudm_ao(nao*2,nao*2))
      
      ! check if WF = DHF otherwise go on for CCSD 
      kudm_mo = DC0
      if (wf .eq. 'DHF') then
          ! add hf contribution

          if (exacc_file%status == -1) then
             exacc_file%type = 2
             exacc_file%name = "EXACC.h5"
             exacc_file%status = 0
          end if

          call lab_read(exacc_file,'nocc',idata=nocc_ccsd)
          allocate(mo_occ(nocc_ccsd))
          call lab_read(exacc_file,'mo_occ',idata=mo_occ)

          ninact = nel - nocc_ccsd
          do i=1,ninact
              kudm_mo(i,i) = kudm_mo(i,i) + DC1
          enddo
          do i=1,nocc_ccsd
              j = mo_occ(i)
              kudm_mo(j,j) = kudm_mo(j,j) + DC1
          enddo

          write(*,*) "   Expectation values at the KU DHF level"

      else ! (wf .eq. 'CCSD') then

          write(*,*) "   Expectation values at the KU X2C-CCSD level"

          ! allocating CC density intermediate arrays        
          allocate(rkudmvec(nmo**2))
          allocate(ikudmvec(nmo**2))
          allocate(ckudmvec(nmo**2))

          ! checking existence of cc density file
          if (exacc_file%status == -1) then
              exacc_file%type = 2
              exacc_file%name = "EXACC.h5"
              exacc_file%status = 0
          end if

     
          ! reading in cc density from file
          rkudmvec = D0
          ikudmvec = D0
          call lab_read (exacc_file,'CCDENS_RE',rdata=rkudmvec)
          call lab_read (exacc_file,'CCDENS_IM',rdata=ikudmvec)
          ckudmvec = DC0
          ckudmvec = cmplx(rkudmvec,ikudmvec,kind=kind(1.0d0))
          deallocate(rkudmvec)
          deallocate(ikudmvec)
          kudm_mo = DC0
          kudm_mo(1:nmo,1:nmo) = 
     &           reshape(ckudmvec,shape=(/nmo,nmo/),order=(/2,1/))
          deallocate(ckudmvec)
         
          ! add hf contribution 
          call lab_read (exacc_file,'nocc',idata=nocc_ccsd)
          allocate(mo_occ(nocc_ccsd))
          call lab_read (exacc_file,'mo_occ',idata=mo_occ)
          ninact = nel - nocc_ccsd
          do i=1,ninact
              kudm_mo(i,i) = kudm_mo(i,i) + DC1
          enddo
          do i=1,nocc_ccsd
              j = mo_occ(i)
              kudm_mo(j,j) = kudm_mo(j,j) + DC1
          enddo
          
      endif ! (WF == CCSD... see upper WF specifier for density)

      ! transforming density matrix from mo to ao (or x2c) basis
      kudm_ao  = DC0
      kudm_tmp = DC0
      kudm_tmp = matmul(ccoeffmat,kudm_mo)
      deallocate(kudm_mo)
      kudm_ao = matmul(kudm_tmp,transpose(dconjg(ccoeffmat)))
      deallocate(kudm_tmp)
      deallocate(ccoeffmat)
C
C
C     Obtaining Expectation Values
C     ============================
C
      if (tobe_dirac) then

          ! allocating property integral arrays
          if (x2c) then
              allocate(qx2cmatvec(nao**2 * 4))
              allocate(qx2cmat(nao, nao, 4))
              allocate(qqx2cmatvec(nao**2,4))
          else if (nonrel) then
              allocate(scr_mat1(nao**2,4))
              allocate(scr_mat2(nao* nao * 4))
          else
              allocate(aopropervec(nao*(nao+1)/2))
          endif

          ! checking existence of property matrices file 
          if (prp_hdf5) then
          write(*,*) 'Opening X2CMAT HDF5 file for property integrals'
              if (x2cmat_file%status == -1) then
                  x2cmat_file%type = 2
                  x2cmat_file%name = "X2CMAT.h5"
                  x2cmat_file%status = 0
              end if
          endif

      elseif (tobe_respect) then

          ! allocating property integral arrays
          allocate(rprpvec(nao**2 * 4))
          allocate(iprpvec(nao**2 * 4))
          allocate(cprpvec(nao**2 * 4))
          reex(1) = 'sphRx '
          reex(2) = 'sphRy '
          reex(3) = 'sphRz '

      endif ! (tobe...)
      allocate(cprpmat(nao*2, nao*2))
      
      CALL TITLER('Expectation values','*',127)

      do i = 1, nexpp

          indxpr = lexpp(i)
          !reading in property matrices from file
          if (tobe_dirac) then
              ! Either X2C or NONREL Hamiltonians!
              if (x2c) then
                  !Getting integrals from X2CMAT!
                  call num2str(indxpr,numstr)
                  prplab = 'prpint2c'//numstr
                  qqx2cmatvec = D0
                  call x2c_read(prplab,qqx2cmatvec,nao*nao*4,lux2c)
                  qx2cmat = reshape(qqx2cmatvec, shape=(/nao,nao,4/),
     &                              order=(/1,2,3/))
              cprpmat = DC0
              ! converting property matrix from
              ! quaternion to complex format
              call qtoc(qx2cmat,cprpmat,1,nao,nao,nao,nao,4)

              else if (nonrel) then
                  ! Getting integrals from AOPROPER calling PRPMAO!
                  prp_fer_rep = iprpsym(indxpr) - 1
                  prp_triang  = abs(iprltyp(indxpr)).eq.1
                  quat_array_filled(1) = .true.
                  scr_mat1 = 0
                  scr_mat2 = 0
                  if (prp_triang) then
                      prp_mat_dim = nao*(nao+1)/2
                      call prpmao(luprop,indxpr,.true.,scr_mat1,
     &                prp_triang,prp_mat_dim,scr_mat1,scr_mat2,
     &                quat_array_filled,4)
                      call tri2sq(scr_mat1,scr_mat2,quat_array_filled,
     &                            prp_fer_rep,iprptim(indxpr))
                  else
                      prp_mat_dim = nao**2
                      call prpmao(luprop,indxpr,.true.,scr_mat1,
     &                prp_triang,prp_mat_dim,scr_mat1,scr_mat2,
     &                quat_array_filled,4)
                      do j = 1,4
                          if(quat_array_filled(j)) then
                              call dzero(scr_mat1(1,j),prp_mat_dim)
                          endif
                      enddo
                  endif
                  
                  cprpmat = DC0
                  ! converting property matrix from
                  ! quaternion to complex format
                  call qtoc(scr_mat1,cprpmat,1,nao,nao,nao,nao,4)


              else ! non-relativistic reading from AOPROPER!
 
                  call quit("Error in KUPRPEX1, no suitable
     &                       Hamiltonian found!")
              
              endif ! x2c or nonrel

          elseif (tobe_respect) then

              ! reading in property matrices from file
              rprpvec = D0
              iprpvec = D0
              call lab_read (reex_file,reex(i)//'_RE',rdata=rprpvec)
              call lab_read (reex_file,reex(i)//'_IM',rdata=iprpvec)
              cprpvec = DC0
              cprpvec = cmplx(rprpvec,iprpvec,kind=kind(1.0d0))
              cprpmat = DC0
              cprpmat = reshape(cprpvec,shape=(/2*nao,2*nao/),
     &                  order=(/1,2/))
              cprpmat = cprpmat * DCM1

          endif ! (tobe...)

          ! carry out the contraction between prop. and dens. matrices
          expvali = DC0
          do p=1,nao*2
              do q=1,nao*2
                  expvali = expvali + cprpmat(p,q) * kudm_ao(q,p)
              enddo
          enddo

          !Specifying time-reversability and picking Re or Im part.
          if (iprptim(i) == 1) then
              expval(i) = dreal(expvali)
          elseif (iprptim(i) == -1) then
              expval(i) = dimag(expvali)
          endif
          
          ! printing out the expectation value
          FMT = MXFORM(EXPVAL(I),15)
          WRITE(LUPRI,'(4X,A16,A3,3X,'//FMT//',1X,A4)')
     &       PRPNAM(INDXPR),' : ',EXPVAL(I),'a.u.'

      enddo ! (for i = 1, nexpp)

      CALL PRSYMB(LUPRI,'-',75,4)
      WRITE(LUPRI,'(4X,A)')
     &   'Done calculating expectation values'

      !deallocating property (intermediate) arrays
      if (tobe_dirac) then
          if (x2c) then
              deallocate(qx2cmatvec,qqx2cmatvec,qx2cmat)
          else if (nonrel)  then
              deallocate(scr_mat1,scr_mat2)
          else
              deallocate(aopropervec)
          endif
      elseif (tobe_respect) then
          deallocate(rprpvec,iprpvec,cprpvec)
      endif
      deallocate(cprpmat,kudm_ao)

      RETURN
      END


C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpexp */
      SUBROUTINE PRPEXP(WF,EXPVAL,WORK,LWORK,state)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C
C***********************************************************************
!      use mxcent_mod, only: MXCENT
!      use dcbexp_mod, only: NEXPP
!      use dcbbas_mod, only: N2BBASXQ
!      use dcborb_mod, only: N2ORBXQ, NASHT, NORBT
!      use dgroup_mod, only: NZ
#include "implicit.h"
#include "mxcent.h"
C
      DIMENSION EXPVAL(*),WORK(LWORK)
      CHARACTER WF*4
      INTEGER, OPTIONAL :: state 
      INTEGER :: state_local 
C
#include "dcbexp.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      CALL QENTER('PRPEXP')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('REAL',KCMO  ,N2BBASXQ  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG  ,NORBT     ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KEXPST,2*NEXPP   ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
CHJ TODO : check if 2 DMAT or NOPEN+1 DMAT is the correct ...
      CALL MEMGET('REAL',KDMAT ,max(N2ORBXQ*2,nasht**2*nz*2,N2BBASXQ*2),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPMAT ,max(N2ORBXQ,nasht**2*nz,N2BBASXQ),
     &            WORK,KFREE,LFREE)
#ifdef ANALYZE_PROPERTY_GRADIENT
      LUNIT = 14
      CALL OPNFIL(LUNIT,'PMAT','UNKNOWN','PAMEXP')
#endif
      if (present(state)) then
          state_local = state
      else
          state_local = 1
      end if
      CALL PRPEX1(WF,EXPVAL,WORK(KEXPST),WORK(KCMO),WORK(KEIG),
     &            WORK(KDMAT),WORK(KPMAT),WORK(KFIRST),
     &            WORK(KFREE),LFREE,state_local)
      
#ifdef ANALYZE_PROPERTY_GRADIENT
      CLOSE(LUNIT,STATUS = 'KEEP')
#endif
C
C     Memory deallocation
C
      CALL MEMREL('PRPEXP',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('PRPEXP')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex1 */
      SUBROUTINE PRPEX1(WF,EXPVAL,LEXPST,CMO,EIG,DMAT,PRPAO,FIRST,
     &                  WORK,LWORK,state)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C     Last modifications: M.Ilias, March2006, Strasbourg
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use mxcent_mod, only: MXCENT
!      use dgroup_mod, only: FREP, IPQTOQ, NFSYM, NZ
!      use dcborb_mod, only: I2ORBX, ICMOQ, NASHT, NCMOTQ, NOCC, NORB, NORBT
!      use dcbbas_mod, only: I2BASX, N2BBASXQ, NFBAS, NTBAS
!      use dcbexp_mod, only: IPREXP, LEXPP, NEXPP, ORBANA, PRJANA, PRPCAN
!      use dcbgen_mod, only: LUCOEF
!      use dcbxpr_mod, only: IPRPSYM, IPRPTIM, PRPNAM
!      use dcbdhf_mod, only: AOC, DF
!      use dcbprj_mod, only: NREFS
#include "implicit.h"
      INTEGER :: I, J
#include "mxcent.h"
      PARAMETER(D0 = 0.0D0,D1=1.0D0,DTOL=1.0D-16,D2=2.0D0)      
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbexp.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dcbdhf.h"
#include "dcbprj.h"
      LOGICAL LEXPST(2,*),FILEX,FIRST(NZ)
      DIMENSION EXPVAL(*),CMO(*),DMAT(N2BBASXQ,*),PRPAO(N2BBASXQ),
     &          EIG(*),WORK(*),NP(2,0:2)
      integer  :: lwork
      CHARACTER MXFORM*6,FMT*6
      CHARACTER WF*4
      INTEGER   GETIOPEN
      LOGICAL FNDLAB
      INTEGER :: state 
C
C     Initialize
C     ==========
C
#include "memint.h"
      KFRSAV = KFREE
      CALL DZERO(EXPVAL,NEXPP)
      N2EXP = 2*NEXPP
      CALL LSET(N2EXP,.FALSE.,LEXPST)
C
C     Generate density matrix
C     =======================
C
      IF (WF .EQ. 'DHF' .OR. WF .EQ. 'DFT') THEN
        CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IDUM,TOTERG,6)
        IF (.NOT. ORBANA) THEN
C.........accumulate all density into one matrix
          FILEX = AOC
          AOC = .FALSE.
          CALL DENMAT(DMAT,CMO,IPREXP)
          AOC = FILEX
        END IF
      ELSE IF (WF .EQ. 'MP2 ' .OR. WF .EQ. 'CCSD') THEN

         IF (ORBANA) 
     &        CALL QUIT('*** ERROR in PRPEX1 *** '
     &        // 'ORBANA not implemented for MP2')
C
C        Get relaxed density matrix computed by RELCCSD or EXACORR.
C
         CALL GET_CC_DENSITY (WF,'AO',DMAT,IERR)

      ELSE IF (WF .EQ. 'EOM ') THEN

         IF (ORBANA) 
     &        CALL QUIT('*** ERROR in PRPEX1 *** '
     &        // 'ORBANA not implemented for EOMCC')

C        Get relaxed density matrix computed for excited states

         write(*,*)' Number of EOM states',state

         CALL GET_EOM_DENSITY ('AO',STATE,DMAT,IERR)

      ELSE IF (WF .EQ. 'KRMC') THEN
C
         IF (ORBANA) 
     &        CALL QUIT('*** ERROR in PRPEX1 *** '
     &        // 'ORBANA not implemented for MCSCF')
C
         INQUIRE(FILE='KRMCSCF', EXIST=FILEX)
C
         IF (.NOT. FILEX) 
     &      CALL QUIT('*** ERROR in PRPEX1 *** ' //
     &        'file "KRMCSCF" not found!')
C
         call dzero(dmat,N2BBASXQ*2)

         CALL OPNFIL(LUCOEF,'KRMCSCF','OLD','PRPEXP')
         CALL REAKRMC(LUCOEF,'NEWORB  ',CMO,NCMOTQ)
         CLOSE(LUCOEF,STATUS='KEEP')

         INQUIRE(FILE='MCDENS', EXIST=FILEX)
C
         IF (.NOT. FILEX) 
     &      CALL QUIT('*** ERROR in PRPEX1 *** ' //
     &        'file "MCDENS" not found!')
         CALL OPNFIL(LUCOEF,'MCDENS','OLD','PRPEXP')
         CALL REAKRMC(LUCOEF,'MC  DENS',DMAT,nasht**2*nz)
!        total DV = DV_act + DV_core
         IF (NASHT .GT. 0) CALL DVSCT(DMAT(1,1),DMAT(1,2)) ! DV_act
         call dzero(PRPAO,N2BBASXQ)
         CALL DCSCT(PRPAO) ! DV_core
         call daxpy(N2BBASXQ,1.0D0,PRPAO,1,DMAT(1,2),1)
         CLOSE(LUCOEF,STATUS='KEEP')
C
C        Get AO density matrix.
!        input: C == U+ A U
C
         DO I = 1, NFSYM
            CALL QTRANS('MOAO','S',D0,
     &           NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &           DMAT(1+I2BASX(I,I),1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &           DMAT(1+I2ORBX(I,I),2),
     &           NORBT,NORBT,NZ,IPQTOQ(1,0),
     &           CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           CMO(1+ICMOQ(I)),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK,LWORK,0)
         END DO
C
C        Result is scaled by 2 in PRPEX2...
C
         CALL DSCAL(N2BBASXQ,0.5d0,DMAT,1)
      ELSE
         WRITE(LUPRI,'(1X,2A)')
     &        '*** ERROR in PRPEX1 *** Unknown wave function: ',WF
         CALL QUIT('*** ERROR in PRPEX1 *** Unknown wave function')
      END IF
C
C     Projection analysis
C     ===================
C
      IF(PRJANA) THEN
        INQUIRE(FILE='DFPRJC',EXIST=FILEX)
        IF(.NOT.FILEX) THEN
          WRITE(LUPRI,'(A)')
     &      '** WARNING ! No projection coefficients found !',
     &      '    - Skipping expectation value projection analysis !'
          PRJANA = .FALSE.
          KRVEC = KFREE
          KBVEC = KFREE
          KPRMT = KFREE
        ELSE
          CALL OPNFIL(LUCOEF,'DFPRJC','OLD','PRPEXP')
          READ(LUCOEF) NREFS
          REWIND LUCOEF
          CALL MEMGET('INTE',KNR,NFSYM*NREFS,WORK,KFREE,LFREE)
          CALL REAPRJ(LUCOEF,NP,WORK(KNR),KRVEC,KBVEC,
     &                WORK,KFREE,LFREE)
          CLOSE(LUCOEF,STATUS='KEEP')
          NBMAX = 0
          NPFMT = 0
          DO IFRP = 1,NFSYM
            NBMAX = MAX(NBMAX,NP(IFRP,0))
            NPFMT = NPFMT + NP(IFRP,0)*NP(IFRP,0)*NZ
          ENDDO
          NN = NBMAX*(NBMAX+1)/2
          CALL MEMGET('REAL',KPFMT,NPFMT,WORK,KFREE,LFREE)       
          CALL MEMGET('REAL',KDFMT,NPFMT,WORK,KFREE,LFREE)       
          CALL MEMGET('INTE',KIBUF,NN,WORK,KFREE,LFREE)
          CALL MEMGET('INTE',KINDX,NBMAX*2,WORK,KFREE,LFREE)
          NTOT = NREFS + 1
          NFDIM = NTOT*NTOT
          CALL MEMGET('REAL',KFMAT,NFDIM,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KPM  ,NTOT,WORK,KFREE,LFREE)
C.........Calculate density matrix in fragment basis
          CALL PRJEX4(WORK(KDFMT),WORK(KBVEC),NP(1,0),NP(1,2),IPREXP,
     &                WORK,KFREE,LFREE)          
        ENDIF
      ELSE
        KRVEC = KFREE
        KBVEC = KFREE
        KPFMT = KFREE
      ENDIF
C
C     Expectation values normal procedure
C     ===================================
C
      IF(.NOT.ORBANA) THEN

        CALL TITLER('Expectation values','*',127)

        DO I = 1,NEXPP
          INDXPR = LEXPP(I)
          IREP = IPRPSYM(INDXPR)-1
          ITIM = IPRPTIM(INDXPR)
          CALL PRPEX2(INDXPR,EXPVAL(I),.true.,LEXPST(1,I),
     &                PRPAO,FIRST,WORK,KFREE,LFREE,IPREXP)
          IF(LEXPST(1,I).OR.LEXPST(2,I)) GOTO 30
C.............calculate expectation value
          CALL PRPEX3(EXPVAL(I),INDXPR,DMAT,PRPAO,FIRST,
     &        WORK,KFREE,LFREE,IPREXP)
 30       CONTINUE

          IF(PRJANA) THEN
C.........Transform property matrix to fragment orbital basis
            CALL PRJEX1(WORK(KPFMT),PRPAO,WORK(KRVEC),
     &                  INDXPR,IREP,ITIM,NP(1,0),FIRST,IPREXP,
     &                  WORK,KFREE,LFREE)
          ENDIF
          FMT = MXFORM(EXPVAL(I),15)
          WRITE(LUPRI,'(4X,A16,A3,3X,'//FMT//',1X,A4,3X,2(3X,A,L1))')
     &       PRPNAM(INDXPR),' : ',EXPVAL(I),'a.u.',
     &     's0 = ',LEXPST(1,I),'t0 = ',LEXPST(2,I)
          CALL FLSHFO(LUPRI)
          IF(PRJANA.AND.ABS(EXPVAL(I)).GT.DTOL) THEN
            CALL PRJEX2(EXPVAL(I),WORK(KPFMT),WORK(KDFMT),WORK(KBVEC),
     &                  WORK(KIBUF),INDXPR,NP(1,0),NP(1,2),
     &                  WORK(KNR),NTOT,WORK(KFMAT),WORK(KPM),
     &                  WORK(KINDX),IPREXP,PRPAO,WORK,KFREE,LFREE)
            CALL PRSYMB(LUPRI,'-',75,4)
          ENDIF
          IF(PRPCAN) THEN
            CALL EXPCAN(INDXPR,PRPAO,CMO,IPREXP,WORK,KFREE,LFREE)
          ENDIF
        ENDDO

        CALL PRSYMB(LUPRI,'-',75,4)
        WRITE(LUPRI,'(4X,A)')
     &   's0 = T : Expectation value zero by point group symmetry.'
        WRITE(LUPRI,'(4X,A)')
     &   't0 = T : Expectation value zero by time reversal symmetry.'
C
C     Expectation values; individual orbital contributions
C     ====================================================
C
      ELSE
        CALL TITLER('Expectation values','*',127)
        DO I = 1,NEXPP
          INDXPR = LEXPP(I)
          IREP = IPRPSYM(INDXPR)-1
          ITIM = IPRPTIM(INDXPR)
C          CALL PRSYMB(LUPRI,'-',85,4)
          CALL HEADER('Operator '//PRPNAM(INDXPR)//':',0)            
          CALL HEADER('Expectation value for individual orbitals',-1)
          WRITE(LUPRI,'(/16X,A)') 'Matrix element     Occ.  '
          EXPVAL(I) = D0
          DO 10 IFRP = 1,NFSYM
            IF(NORB(IFRP).EQ.0) GOTO 10
C...........generate property matrix in SO-basis
            CALL PRPEX2(INDXPR,EXPBUF,.true.,LEXPST(1,I),PRPAO,FIRST,
     &                  WORK,KFREE,LFREE,IPREXP)
            IF(LEXPST(1,I).OR.LEXPST(2,I)) GOTO 20
            DO J = 1,NOCC(IFRP)
              CALL DENORB(DMAT,J,1,IFRP,CMO,IPREXP)
              EXPBUF = D0
C.............calculate expectation value
              CALL PRPEX3(EXPBUF,INDXPR,DMAT,PRPAO,FIRST,
     &          WORK,KFREE,LFREE,IPREXP)
              IOPEN = GETIOPEN(J,IFRP)
              EXPVAL(I) = EXPVAL(I) + EXPBUF*DF(IOPEN)
              FMT = MXFORM(EXPBUF,15)
             WRITE(LUPRI,'(4X,A3,I5,3X,'//FMT//',5X,F6.4,5X,'//FMT//')')
     &          FREP(IFRP),J,EXPBUF/D2,D2*DF(IOPEN),EXPBUF*DF(IOPEN)
              CALL FLSHFO(LUPRI)
            ENDDO

            IF(PRJANA) THEN
C...........Transform property matrix to fragment orbital basis
              CALL PRJEX1(WORK(KPFMT),PRPAO,WORK(KRVEC),
     &                    INDXPR,IREP,ITIM,NP(1,0),FIRST,IPREXP,
     &                    WORK,KFREE,LFREE)
            ENDIF
 10       CONTINUE
 20       CONTINUE
          CALL PRSYMB(LUPRI,'-',85,4)
          FMT = MXFORM(EXPVAL(I),15)
          WRITE(LUPRI,'(4X,A8,34X,'//FMT//',1X,A4,3X,2(3X,A,L1))')
     &     'Total : ',EXPVAL(I),'a.u.',
     &     's0 = ',LEXPST(1,I),'t0 = ',LEXPST(2,I)
          CALL PRSYMB(LUPRI,'-',85,4)
          IF(PRPCAN) THEN
            CALL EXPCAN(INDXPR,PRPAO,CMO,IPREXP,WORK,KFREE,LFREE)
          ENDIF
          IF(PRJANA.AND.ABS(EXPVAL(I)).GT.DTOL) THEN
            CALL PRJEX2(EXPVAL(I),WORK(KPFMT),WORK(KDFMT),WORK(KBVEC),
     &                  WORK(KIBUF),INDXPR,NP(1,0),NP(1,2),
     &                  WORK(KNR),NTOT,WORK(KFMAT),WORK(KPM),
     &                  WORK(KINDX),IPREXP,PRPAO,WORK,KFREE,LFREE)
            CALL PRSYMB(LUPRI,'-',85,4)
          ENDIF
        ENDDO
        WRITE(LUPRI,'(4X,A)')
     &   's0 = T : Expectation value zero by point group symmetry.'
        WRITE(LUPRI,'(4X,A)')
     &   't0 = T : Expectation value zero by time reversal symmetry.'
      ENDIF
      IF(PRJANA) THEN
        CALL MEMREL('PRPEX1.prjana',WORK,1,KFRSAV,KFREE,LFREE)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex2 */
      SUBROUTINE PRPEX2(INDXPR,EXPVAL,ONLY_TS,LEXPST,PRPAO,FIRST,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     if ONLY_TS = .true.
C        for totally symmetric operators return property matrix
C        for other symmetries return .FALSE. in
C            LEXPST(1) - molecular point group symmetry
C        or  LEXPST(2) - time reversal symmetry
C
C     if ONLY_TS = .false.
C        return property matrix anyway (i need it in OpenRSP --radovan)
C
C     On output: PRPAO - the complete propery matrix in SA-AO basis !
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C                   Strasbourg/March 2006 - Miro ILIAS
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use mxcent_mod, only: MXCENT
!      use dcbxpr_mod, only: IPRPSYM, IPRPTIM, PRPNAM
!      use dcbexp_mod, only: PRJANA
!      use dcbbas_mod, only: N2BBASX
!      use dgroup_mod, only: NZ
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D0 = 0.0D0, D1=1.0D0
C
#include "mxcent.h"
#include "dcbxpr.h"
#include "dcbexp.h"
#include "dcbbas.h"
#include "dgroup.h"
C
!     radovan:                    ONLY_TS: allow only totally symmetric
      LOGICAL FIRST(NZ),LEXPST(2),ONLY_TS
      DIMENSION PRPAO(N2BBASX,NZ),WORK(*)

      CALL QENTER('PRPEX2')

      EXPVAL = D0
C
C     Check symmetry
C     ==============
C
C     1. Molecular point group symmetry
C
      IREP = IPRPSYM(INDXPR)-1
      IF(IREP.NE.0) LEXPST(1) = .TRUE.
C
C     2. Time reversal symmetry
C
      ITIM = IPRPTIM(INDXPR)
      IF(ITIM.EQ.-1) LEXPST(2) = .TRUE.

CMI ... control output...
      IF (IPRINT.GE.5) THEN
       WRITE(LUPRI,'(/2X,A,/1X,A,I2,1X,A,A,2I2)')
     & '*** Output from PRPEX2  ***: ',
     & ' Operator index and name:',
     &  INDXPR,PRPNAM(INDXPR),
     & ' IREP,ITIM:',IREP,ITIM
       IF (LEXPST(1)) WRITE(LUPRI,'(1X,A)')
     & ' Operator not symmetric! Do not calculate exp.value!'
      ENDIF
C
C     If not totally symmetric, return with zero value
C
      IF (ONLY_TS .AND. LEXPST(1).AND.(.NOT.PRJANA)) GOTO 90 
  109 CONTINUE
!.q 
C
C     Operator to be read from file
C     =============================
C
      CALL PRPMSAO(INDXPR,PRPAO,.TRUE.,WORK,FIRST,
     &             WORK(KFREE),LFREE,IPRINT)
C     
C
 90   CONTINUE     
      CALL QEXIT('PRPEX2')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpex3 */
      SUBROUTINE PRPEX3(EXPVAL,INDXPR,DMAT,PRPAO,FIRST,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Calculate expectation values defined in /CBIEXP/
C
C     Written by Trond Saue May 27 1996
C     Last revision May 27 1996 - tsaue
C     MI entry for LAO-project, May 2003
C
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dcbxpr_mod, only: IPRPSYM, IPRPTIM, PRPNAM
!      use dcbbas_mod, only: N2BBASX, NTBAS
!      use dgroup_mod, only: IPQTOQ, NZ
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D0 = 0.0D0
C
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      LOGICAL FIRST(NZ),LEXPST(2)
      DIMENSION DMAT(N2BBASX,NZ),WORK(*),
     &          PRPAO(N2BBASX,NZ)

      CALL QENTER('PRPEX3')

      IREP = IPRPSYM(INDXPR)-1
      ITIM = IPRPTIM(INDXPR)
      IF (IPRINT.GE.10) THEN
        WRITE(LUPRI,'(/A,A16)')
     & 'PRPEX3: Entering total SA-AO matrix of property ',PRPNAM(INDXPR)
        write(lupri,*) '   proerty operator IREP=',IREP,' ITIM=',ITIM
        CALL PRQMAT(PRPAO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &                 IPQTOQ(1,IREP),LUPRI)
      ENDIF        

C========================================================================
C     Get expectation value, using ONLY the matrices containing
C     integrals (ensured by logical FIRST(1..NZ) )
C========================================================================
      DO IZ = 1,NZ
      IF(.NOT.FIRST(IZ)) THEN
        EXPVAL = EXPVAL + DDOT(N2BBASX,DMAT(1,IZ),1,PRPAO(1,IZ),1)
        IF (IPRINT.GE.5) THEN
          WRITE(lupri,*)
          WRITE(LUPRI,*)
     &       ' == PRPEX3: accumulated EXPVAL (ORBANA=.false.)'//
     &       ' after IZ=',IZ,
     &       '/',NZ,' exp.val=',2*EXPVAL
          WRITE(LUPRI,'(/2X,A,I1)') 'PRPEX3: PRPAO of IZ=',IZ
          CALL OUTPUT(PRPAO(1,IZ),1,NTBAS(0),1,NTBAS(0),
     &              NTBAS(0),NTBAS(0),1,LUPRI)
          WRITE(LUPRI,'(/2X,A,I1)') 'PRPEX3: DMAT of IZ=',IZ
          CALL OUTPUT(DMAT(1,IZ),1,NTBAS(0),1,NTBAS(0),
     &              NTBAS(0),NTBAS(0),1,LUPRI)
        ENDIF
      ENDIF
      ENDDO
      EXPVAL = EXPVAL + EXPVAL
C
      CALL QEXIT('PRPEX3')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prjex1 */
      SUBROUTINE PRJEX1(PFO,PAO,RVEC,INDXPR,IREP,ITIM,NCNT,
     &                  FIRST,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Transform property matrix to fragment orbital basis
C
CMI   On input:  PAO - property matrix in SA-AO basis
C
C     Written by Trond Saue Apr 27 2003 
C     Last modifications: M.Ilias, March 2006, Strasbourg
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dcbxpr_mod, only: IPRPTYP, PRPNAM
!      use dcbbas_mod, only: N2BBASX, NFBAS, NTBAS
!      use dgroup_mod, only: IPQTOQ, NFSYM, NZ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: II
      REAL(dp), PARAMETER :: D0 = 0.0D0, D1=1.0D0
C
#include "dcbxpr.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      LOGICAL FIRST(NZ)
      DIMENSION PFO(*),PAO(N2BBASX,NZ),NCNT(*),RVEC(*),WORK(*)

      CALL QENTER('PRJEX1')
C
C     Transform on the fly
C
      IOPTYP = IPRPTYP(INDXPR)
      NDIM = 0
      DO IFRP = 1,NFSYM
        NDIM = NDIM + NCNT(IFRP)*NCNT(IFRP)*NZ
      ENDDO
      CALL DZERO(PFO,NDIM)
      FADD = D0

      DO IZ = 1,NZ

      IF(.NOT.FIRST(IZ)) THEN
        II = IZ
        IPA = 1
        IPF = 1
        IPR = 1
        DO IFRP = 1,NFSYM
          CALL QTRANS('AOMO','S',FADD,
     &         NFBAS(IFRP,0),NFBAS(IFRP,0),NCNT(IFRP),NCNT(IFRP),
     &         PAO(IPA,II),NTBAS(0),NTBAS(0),1,IPQTOQ(IZ,IREP),
     &         PFO(IPF),NCNT(IFRP),NCNT(IFRP),NZ,IPQTOQ(1,IREP),
     &         RVEC(IPR),NFBAS(IFRP,0),NCNT(IFRP),NZ,IPQTOQ(1,0),
     &         RVEC(IPR),NFBAS(IFRP,0),NCNT(IFRP),NZ,IPQTOQ(1,0),
     &         WORK(KFREE),LFREE,IPRINT)
          IPA = IPA + NFBAS(IFRP,0)*(NTBAS(0)+1)
          IPF = IPF + NCNT(IFRP)*NCNT(IFRP)*NZ
          IPR = IPR + NFBAS(IFRP,0)*NCNT(IFRP)*NZ
        ENDDO
        FADD = D1
      ENDIF
      ENDDO
C
#ifdef ANALYZE_PROPERTY_GRADIENT
      LUNIT = 14
      WRITE(LUNIT) PRPNAM(INDXPR)
      IPF = 1
      DO IFRP = 1,NFSYM
        NDIM = NCNT(IFRP)*NCNT(IFRP)*NZ
        CALL WRITT(LUNIT,NDIM,PFO(IPF))
        IPF = IPF + NDIM
      ENDDO
#endif
C
      IF(IPRINT.GE.2) THEN
        IPF = 1
        DO IFRP = 1,NFSYM
          WRITE(LUPRI,'(A,A16,A,I3)')
     &    'Property matrix of ',PRPNAM(INDXPR),
     &    'in fragment basis..',IFRP
          CALL PRQMAT(PFO(IPF),NCNT(IFRP),NCNT(IFRP),
     &         NCNT(IFRP),NCNT(IFRP),NZ,
     &         IPQTOQ(1,IREP),LUPRI)
          IPF = IPF + NCNT(IFRP)*NCNT(IFRP)*NZ
CTROND
          LUBUF = 35
          CALL OPNFIL(LUBUF,'PRPMAT','UNKNOWN','PAMEXP')
          NOST = NCNT(IFRP)*NCNT(IFRP)
          CALL WRITT(LUBUF,NOST,PFO)
          CLOSE(LUBUF,STATUS='KEEP')
        ENDDO       
      ENDIF
C
      CALL QEXIT('PRJEX1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prjex2 */
      SUBROUTINE PRJEX2(EXPVAL,PFO,DFO,BVEC,IBUF,
     &                  INDXPR,NCNT,NMO,
     &                  NR,NTOT,FMAT,PM,INDX,
     &                  IPRINT,PBUF,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Make generalized population matrix
C
C     Written by Trond Saue Apr 27 2003 
C
CMI    MI/March 2006 ... seems that PTRI is needed...
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dgroup_mod, only: NFSYM, NZ
!      use dcbprj_mod, only: NREFS, REFFIL
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I, II, J
      REAL(dp), PARAMETER :: D0 = 0.0D0, D1=1.0D0, D2=2.0D0
C
!#include "dcbxpr.h"
#include "dgroup.h"
#include "dcbprj.h"
C
      DIMENSION PFO(*),DFO(*),NCNT(*),NMO(*),BVEC(*),IBUF(*),PBUF(*),
     &          NR(NREFS,NFSYM),FMAT(NTOT,NTOT),PM(NTOT),INDX(2,*),
     &          WORK(*)
      CHARACTER MXFORM*6,FMT*6
C
      KFRSAV=KFREE
      N2TOT = NTOT*NTOT
      CALL DZERO(FMAT,N2TOT)
      CALL DZERO(PM,NTOT)
C
C     Make population matrix in fragment orbitals
C
      IPD = 1
      DO 10 IFRP = 1,NFSYM
        IF(NMO(IFRP).EQ.0) GOTO 10
        II = 0
        DO I = 1,NREFS
          DO J = 1,NR(I,IFRP)
            II = II + 1
            INDX(1,II) = I
            INDX(2,II) = J
          ENDDO
        ENDDO
        I = NREFS + 1
        DO J = 1,NMO(IFRP)
          II = II + 1
          INDX(1,II) = I
          INDX(2,II) = J
        ENDDO
        NTRI = NCNT(IFRP)*(NCNT(IFRP)+1)/2
        CALL MEMGET('REAL',KPTRI,NTRI,WORK,KFREE,LFREE)
        CALL PRJEX3(EXPVL2,PFO(IPD),DFO(IPD),WORK(KPTRI),NCNT(IFRP),
     &              NR(1,IFRP),NMO(IFRP),IBUF,PBUF,FMAT,PM,INDX,
     &              NTOT,INDXPR,IFRP,IPRINT)
        CALL MEMREL('PRJEX2.ex3',WORK,1,KPTRI,KFREE,LFREE)
        IPD = IPD + NCNT(IFRP)*NCNT(IFRP)*NZ
 10   CONTINUE
      CALL DSCAL(NTOT,D2,PM,1)
C
C     Give fragment contributions
C
      CALL HEADER('* Total fragment contributions',-1)
      WRITE(LUPRI,'(A)') ' * Intraatomic contributions:'
      TOT = D0
      DO J = 1,NREFS
        TMP = D2*FMAT(J,J)
        TOT = TOT + TMP
        FAC = TMP/EXPVAL
        FMT = MXFORM(TMP,25)
        WRITE(LUPRI,'(I5,2X,A6,16X,'//FMT//',E16.5)') 
     &    J,REFFIL(J),TMP,FAC
        FMT = MXFORM(PM(J),25)
        WRITE(LUPRI,'(9X,A,'//FMT//',E16.5)')
     &    '- principal moment: ',PM(J),PM(J)/EXPVAL
        TMP = TMP - PM(J)
        FMT = MXFORM(TMP,25)
        WRITE(LUPRI,'(9X,A,'//FMT//',E16.5)')
     &    '- hybridization   : ',TMP,TMP/EXPVAL
      ENDDO
      FMT=MXFORM(TOT,25)
      WRITE(LUPRI,'(3X,A,19X,'//FMT//',E16.5)') 
     &      'Total: ',TOT,TOT/EXPVAL 
      WRITE(LUPRI,'(A)') ' * Interatomic contributions:'
      TOT = D0
      DO J = 1,NREFS
        DO I = 1,(J-1)
          TMP = D2*(FMAT(I,J)+FMAT(J,I))
          TOT = TOT + TMP
          FAC = TMP/EXPVAL
          FMT=MXFORM(TMP,25)
          WRITE(LUPRI,'(2(I5,2X,A6),3X,'//FMT//',E16.5)') 
     &      J,REFFIL(J),I,REFFIL(I),TMP,FAC
        ENDDO
      ENDDO
      FMT=MXFORM(TOT,25)
      WRITE(LUPRI,'(3X,A,19X,'//FMT//',E16.5)') 
     &      'Total: ',TOT,TOT/EXPVAL
      WRITE(LUPRI,'(A)') ' * Polarization contributions:'
      J = NREFS+1
      TOT = D0
      DO I = 1,NREFS
        TMP = D2*(FMAT(I,J)+FMAT(J,I))
        TOT = TOT + TMP
        FAC = TMP/EXPVAL
        FMT=MXFORM(TMP,25)
        WRITE(LUPRI,'(2(I5,2X,A6),3X,'//FMT//',E16.5)') 
     &      J,REFFIL(J),I,REFFIL(I),TMP,FAC
      ENDDO
      TMP = D2*FMAT(J,J)
      TOT = TOT + TMP
      FAC = TMP/EXPVAL
      FMT=MXFORM(TMP,25)
      WRITE(LUPRI,'(2(I5,2X,A6),3X,'//FMT//',E16.5)') 
     &    J,REFFIL(J),J,REFFIL(J),TMP,FAC
      FMT=MXFORM(TOT,25)
      WRITE(LUPRI,'(3X,A,19X,'//FMT//',E16.5)') 
     &     'Total: ',TOT,TOT/EXPVAL
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjex3 */
      SUBROUTINE PRJEX3(EXPVAL,PFO,DFO,PTRI,NCNT,
     &                  NR,NMO,IBUF,PBUF,FMAT,PM,INDX,NTOT,
     &                  INDXPR,IFRP,IPRINT)
C***********************************************************************
C
C     Calculate generalized population matrix
C
C     Written by T. Saue Apr 27 2003
C     Last modifications: MI, March2006/Strasbourg
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dgroup_mod, only: FREP, NZ
!      use dcbxpr_mod, only: PRPNAM
!      use dcbprj_mod, only: REFFIL
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I, II, J, JJ
      REAL(dp), PARAMETER :: D2=2.0D0, D0=0.0D0, DTOL=1.0D-8
C
#include "dgroup.h"
#include "dcbxpr.h"
#include "dcbprj.h"
      DIMENSION PFO(NCNT,NCNT,NZ),DFO(NCNT,NCNT,NZ),PTRI(*),
     &          NR(*),PBUF(*),IBUF(*),FMAT(NTOT,NTOT),PM(NTOT),INDX(2,*)
C
C
C     Make population matrix
C
      NTRI = NCNT*(NCNT+1)/2
      CALL DZERO(PTRI,NTRI)
      DO IZ = 1,NZ
        IJ = 0
        DO J = 1,NCNT
          JR = INDX(1,J)
          DO I = 1,(J-1)
            IR = INDX(1,I)
            IJ = IJ + 1
            TMP = D2*PFO(I,J,IZ)*DFO(I,J,IZ)
            FMAT(IR,JR) = FMAT(IR,JR) + TMP
            PTRI(IJ) = PTRI(IJ) + TMP
          ENDDO
          IJ = IJ + 1
          TMP = PFO(J,J,IZ)*DFO(J,J,IZ)
          PTRI(IJ) = PTRI(IJ) + TMP
          PM(JR) = PM(JR) + TMP
          FMAT(JR,JR) = FMAT(JR,JR) + TMP
        ENDDO
      ENDDO      
      CALL DSCAL(NTRI,D2,PTRI,1)
C      CALL DSCAL(NTOT,D2,PM,1)
      EXPVAL = DSUM(NTRI,PTRI,1)
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,'(A,A3/A,A16,E16.10)') 
     &    '* Contribution from fermion ircop ',FREP(IFRP),
     &    '  to expectation value of property ',PRPNAM(INDXPR),EXPVAL
        WRITE(LUPRI,'(A)') '* Corresponding population matrix:'
        CALL PRMUTC(PTRI,NCNT,1,LUPRI)      
      ENDIF
C
C     Make ordered list of contributions
C
!     IF(IPRINT.GE.1) THEN
! We always want this list /hjaaj Sep 2012
        CALL HEADER('Ordered list of contributions',-1)            
        DO IJ = 1,NTRI
          PBUF(IJ) = -ABS(PTRI(IJ))
        ENDDO
        CALL INDEXX (NTRI,PBUF,IBUF)
        TMP = D0      
        WRITE(LUPRI,'(2(A4,10X),2A16,A10)')
     &    'Ind1','Ind2','Contribution','Accumulated','Fraction'
        DO M = 1,NTRI
          IJ = IBUF(M)
          IF(ABS(PTRI(IJ)/EXPVAL).LT.DTOL) GOTO 10
          J  = INT(SQRT(dble(2*IJ) + 0.25D0) + 0.4999D0)
          I  = IJ - J*(J-1)/2
          JR = INDX(1,J)
          JJ = INDX(2,J)
          IR = INDX(1,I)
          II = INDX(2,I)
          TMP = TMP + PTRI(IJ)
          FAC = TMP/EXPVAL
          WRITE(LUPRI,'(2(A6,I5,3X),1P,2E16.8,E16.5)')
     &          REFFIL(JR),JJ,REFFIL(IR),II,PTRI(IJ),TMP,FAC
        ENDDO
 10     CONTINUE
!     ENDIF
C      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck prjex4 */
      SUBROUTINE PRJEX4(DFO,BVEC,NFO,NMO,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Make density matrix in fragment basis
C     It may be that the user has selected orbitals to analyze that
C     does not coincide with all occupied ones; this is taken into
C     account below      
C
C     NFO(IFRP) - number of fragment orbitals in fermion ircop IFRP
C     NMO(IFRP) - number of molecular orbitals in fermion ircop IFRP
C
C     Written by Trond Saue Oct 21 2003 
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dgroup_mod, only: IPQTOQ, NFSYM, NZ
!      use dcborb_mod, only: NESH, NPSH
!      use dcbprj_mod, only: VECPRJ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: K
      REAL(dp), PARAMETER :: D0 = 0.0D0, D1=1.0D0, D2=2.0D0
C
#include "dgroup.h"
#include "dcborb.h"
#include "dcbprj.h"      
      DIMENSION DFO(*),NFO(*),NMO(*),BVEC(*),WORK(*)
      DIMENSION NSTR(2,0:2,0:0),IQ(4)
      integer, allocatable :: JMOL(:,:)
      real(8), allocatable :: OCC(:)
      DATA IQ/1,2,3,4/      
C
C     Make density matrix in fragment orbitals
C
      IPD = 1
      IPB = 1
      DO 10 IFRP = 1,NFSYM
      IF(NMO(IFRP).EQ.0) GOTO 10
        NSTR(IFRP,1,0) = NESH(IFRP)
        NSTR(IFRP,2,0) = NPSH(IFRP)
        CALL ORBNUM(VECPRJ(IFRP),IFRP,KVEC,NSTR(1,0,0),
     &       WORK,KFREE,LFREE)
        IF(NSTR(IFRP,2,0).NE.0) THEN
          CALL QUIT('PRJEX4: Do not know how to handle '//
     &            'negative-energy orbs')
        ENDIF
        allocate(OCC(NMO(IFRP)))
        allocate(JMOL(3,NMO(IFRP)))
        CALL PRJINF(JMOL,IFRP,0,0,NSTR(1,0,0),WORK(KVEC),NPSH(IFRP))
        CALL SETOCC(NMO(IFRP),OCC,JMOL)
        deallocate(JMOL)
        FAC  = D0
        IPB1 = IPB
        DO K = 1,NMO(IFRP)
          CALL QGEMM(NFO(IFRP),NFO(IFRP),1,OCC(K),
     &     'N','N',IQ,BVEC(IPB1),NFO(IFRP),NMO(IFRP),NZ,
     &     'H','N',IQ,BVEC(IPB1),NFO(IFRP),NMO(IFRP),NZ,
     &          FAC,IQ,DFO(IPD),NFO(IFRP),NFO(IFRP),NZ)
          FAC  = D1
          IPB1 = IPB1 + NFO(IFRP)
        ENDDO
        deallocate(OCC)
        IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,'(A,I3)')
     &    'Density matrix in fragment basis..',IFRP
          CALL PRQMAT(DFO(IPD),NFO(IFRP),NFO(IFRP),
     &         NFO(IFRP),NFO(IFRP),NZ,
     &         IPQTOQ(1,0),LUPRI)
        ENDIF
        IPB = IPB + NFO(IFRP)*NMO(IFRP)*NZ
        IPD = IPD + NFO(IFRP)*NFO(IFRP)*NZ
 10   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Exploc */
      SUBROUTINE EXPCAN(INDXPR,PAO,CMO,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate property canonical occupied orbitals,
C     that is orbitals that diagonalize the matrix of 
C     property INDXPR
C
C     Input:
C       PAO - property matrix in AO-basis
C       CMO - molecular orbital coefficients
C
C     Written by T. Saue Jan 20 2005 
C
C***********************************************************************
!      use dcbbas_mod, only: I2BASX, NFBAS
!      use dcborb_mod, only: ICMOQ, IORB, NCMOTQ, NISH, NORB, NORBT
!      use dcbgen_mod, only: LUCOEF
!      use dcbxpr_mod, only: PRPNAM
!      use dgroup_mod, only: NFSYM, NZ
       use dircmo
       use labeled_storage
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D0 = 0.0D0
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dcbxpr.h"
#include "dgroup.h"
      DIMENSION PAO(*),CMO(*),WORK(*)
      integer, allocatable :: IBP(:)
      real(8), allocatable :: EIG(:)
      type(file_info_t)    :: fexpcan
      logical              :: tobe
C
      KFRSAV = KFREE
C
C     Read eigenvalues and boson irrep information
C
      allocate (eig(norbt))
      allocate (ibp(norbt))
      CALL REACMO(LUCOEF,'DFCOEF',DUM,EIG,IBP,TOTERG,12)
      DO IFRP = 1,NFSYM
C       Memory allocation
        CALL MEMGET('REAL',KPMO,NISH(IFRP)*NISH(IFRP)*NZ,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KSCO,NISH(IFRP)*NFBAS(IFRP,0)*NZ,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KUCO,NISH(IFRP)*NISH(IFRP)*NZ,
     &              WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KPEI,NISH(IFRP),WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KSEL,NISH(IFRP),WORK,KFREE,LFREE)
C       Select occupied molecular orbitals
        CALL OCCUPY(WORK(KSEL),NISH(IFRP))
        CALL SELCFS(CMO(ICMOQ(IFRP)+1),IFRP,WORK(KSCO),NISH(IFRP),
     &              WORK(KSEL),D0,NISH(IFRP),NFBAS(IFRP,0),NORB(IFRP))
        CALL MEMREL('EXPCAN.sel',WORK,KSEL,KSEL,KFREE,LFREE)
C       Generate canonical property orbitals
        CALL EXPCA1(INDXPR,IFRP,PAO(I2BASX(IFRP,IFRP)+1),
     &              WORK(KPMO),CMO(ICMOQ(IFRP)+1),EIG(1+IORB(IFRP)),
     &              WORK(KSCO),WORK(KUCO),WORK(KPEI),
     &              IPRINT,WORK,KFREE,LFREE)
C      
        CALL MEMREL('EXPCAN.1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      ENDDO
C     Write new coefficients to a separate file called EXPCAN.h5
      fexpcan%name = 'EXPCAN.h5'
      fexpcan%type = 1
      inquire (file=fexpcan%name,exist=tobe)
      if (.not.tobe) then
         fexpcan%status = -1 ! open for the first time
      else
         fexpcan%status = 0  ! add to existing data file
      endif
      call lab_write(fexpcan,trim(PRPNAM(INDXPR))//'_cmo',
     & rdata=cmo(1:ncmotq))
      call lab_write(fexpcan,trim(PRPNAM(INDXPR))//'_eig',rdata=eig)
      call lab_write(fexpcan,trim(PRPNAM(INDXPR))//'_ibp',idata=ibp)
C     Read old coefficients
      CALL REACMO(LUCOEF,'DFCOEF',CMO,DUM,IDUM,TOTERG,2)
      CALL MEMREL('EXPCAN.2',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Explo1 */
      SUBROUTINE EXPCA1(INDXPR,IFRP,PAO,PMO,CMO,EIG,SCO,UCO,PEIG,
     &                  IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Generate property canonical occupied orbitals,
C     that is orbitals that diagonalize the matrix of 
C     property INDXPR
C
C     Input:
C       PAO - property matrix in AO-basis
C       CMO - molecular orbital coefficients
C
C     Written by T. Saue Jan 20 2005 
C
C***********************************************************************
      use priunit_mod, only: LUPRI
!      use dcbbas_mod, only: NFBAS, NTBAS
!      use dcborb_mod, only: NISH, NORB, NPSH
!      use dgroup_mod, only: FREP, IPQTOQ, NZ
      use kinds_mod, only: dp
#include "implicit.h"
      INTEGER :: I
      REAL(dp), PARAMETER :: D0 = 0.0D0, D1 = 1.0D0
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      CHARACTER MXFORM*6,FMT*6
      DIMENSION PAO(*),PMO(*),SCO(*),UCO(*),CMO(*),EIG(*),
     &          PEIG(*),WORK(*)
C
      KFRSAV = KFREE
C
C     Generate property matrix
C
      CALL QTRANS('AOMO','S',D0,NFBAS(IFRP,0),NFBAS(IFRP,0),
     &            NISH(IFRP),NISH(IFRP),
     &            PAO,NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            PMO,NISH(IFRP),NISH(IFRP),NZ,IPQTOQ(1,0),
     &            SCO,NFBAS(IFRP,0),NISH(IFRP),NZ,IPQTOQ(1,0),
     &            SCO,NFBAS(IFRP,0),NISH(IFRP),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)
      IF(IPRINT.GE.1) THEN
        WRITE(LUPRI,*) '* Property matrix:'
        CALL PRQMAT(PMO,NISH(IFRP),NISH(IFRP),NISH(IFRP),NISH(IFRP),
     &            NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Diagonalize
C
      CALL QDIAG(NZ,NISH(IFRP),PMO,NISH(IFRP),NISH(IFRP),
     &           PEIG,1,UCO,NISH(IFRP),NISH(IFRP),
     &           WORK(KFREE),LFREE,IERR)
C
C     Write eigenvalues
C
      WRITE(LUPRI,*) '* Eigenvalues:'
      TMP = D0
      DO I = 1,NISH(IFRP)
        FMT = MXFORM(PEIG(I),15)
        WRITE(LUPRI,'(4X,A3,I5,3X,'//FMT//',1X,A4)')
     &        FREP(1),I,PEIG(I),'a.u.'
      ENDDO
      TMP = DSUM(NISH(IFRP),PEIG,1)
      FMT = MXFORM(TMP,15)
      WRITE(LUPRI,'(4X,A8,3X,'//FMT//',1X,A4)')
     &     'Total : ',TMP,'a.u.'
      IOFF = NPSH(IFRP)+1
      CALL DCOPY(NISH(IFRP),PEIG,1,EIG(IOFF),1)
C
C     Write eigenvectors
C
      IF(IPRINT.GE.2) THEN
        WRITE(LUPRI,*) '* Eigenvectors:'
        CALL PRQMAT(UCO,NISH(IFRP),NISH(IFRP),NISH(IFRP),NISH(IFRP),NZ,
     &            IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Backtransform
C
      IOFF = NFBAS(IFRP,0)*NPSH(IFRP)+1
      CALL QGEMM(NFBAS(IFRP,0),NISH(IFRP),NISH(IFRP),D1,
     &           'N','N',IPQTOQ(1,0),SCO,NFBAS(IFRP,0),NISH(IFRP),NZ,
     &           'N','N',IPQTOQ(1,0),UCO,NISH(IFRP),NISH(IFRP),NZ,
     &           D0,IPQTOQ(1,0),CMO(IOFF),NFBAS(IFRP,0),NORB(IFRP),NZ)
      RETURN
      END
