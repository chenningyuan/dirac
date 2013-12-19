***********************************************************************

      subroutine Q_correction(C,NBLOCK,IBLOCK,IBLTP,LUC,       ! C,LUSCR, ICISTR raus TODO??? 
     &                        LUSCR,ICISTR,IRC,WORK,relQC)
*
*
*''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
*''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!     in:
!         C: scratch space for ci vector
!         NBLOCK: # of ttsblocks
!         IBLOCK: array containing information about each ttsblock
!         IBLTP: ?
!         LUC: int filepointer
!         LUSCR: filepointer for scrathc file
!         ICISTR: 2
!         IRC: 1 (real) 2 (complex)
*
      use qcorr_cfg
      use qcorr_interface

      implicit real*8(A-H,O-Z)

      dimension IBLOCK(8,NBLOCK),IBLTP(*)
      dimension C(*), WORK(*)
      logical, intent(in)   :: relQC
      integer, allocatable  :: mat1(:,:)
      integer, allocatable  :: astring(:,:)
      integer, allocatable  :: bstring(:,:)
      integer, allocatable  :: occupation_class(:,:)
      real(8), allocatable  :: mat2(:,:)
      character (len=2)     :: ferm_lab
      character (len=200), allocatable   :: ref_wavefunction(:,:)
      real(8), allocatable  :: ref_wavefunction_coeff(:)
      real(8), allocatable  :: ref_wavefunction_coeff_i(:)         ! coeff von der ersten ci_iteration!!!
      integer               :: nref = 1                            ! for HF nref = 1, TODO: input parameter!!!
      integer               :: nstates_ci
      real(8), allocatable  :: e_corr(:), e_tot(:)
      real(8)               :: q_corr, edci
*
      call QENTER('QCORR')

!     initialize
      call init_qcorr(relQC)

      allocate(ref_wavefunction(10,2))
      allocate(ref_wavefunction_coeff(10))
      allocate(ref_wavefunction_coeff_i(10))
      ref_wavefunction_coeff      = 0.0d0 ! initialize all to 0.0d0
      ref_wavefunction_coeff_i    = 0.0d0 ! initialize all to 0.0d0
      ref_wavefunction_coeff_i(1) = 1.0d0
      if(relQC)then
        write(ref_wavefunction(1,1),"(4i4)"),1,2,3,4               ! for HF, TODO: multi reference case!!!
        write(ref_wavefunction(1,2),"(4i4)"),1,2,3,4
      else
        write(ref_wavefunction(1,1),"(4i4)"),1,2,5,9
        write(ref_wavefunction(1,2),"(4i4)"),1,2,5,9
      end if

*
** Specifications of internal space
*
      NTEST = 000
      NTEST = max(NTEST,IPRWFC)
*
      LBLK = -1
      IBOFF = 0
      call REWINE(LUSCR,-1)
      call REWINE(LUC,-1)
      call itods(-1,1,-1,LUSCR)
*
* Loop over real and imaginary parts
      do IRIR = 1,IRC,1                                      !loop über call

        if (NTEST.ge.50) write(6,*) 'right r/i is  ',IRIR
        if (IRIR.EQ.2) then
*
*. Skip EOV mark between real and imaginary part
          call IFRMDS(IONEM,1,-1,LUC)
        end if
        do IMK2_R = 1,NMS2VAL_q,1

          if(relQC)then
            IACTP = IST_FOR_DT_q(1,IMK2_R)
            IBCTP = IST_FOR_DT_q(2,IMK2_R)
          else
            IACTP = 1
            IBCTP = 2
          end if
*
*  Offset for IBLOCK(sym,type) information
          if (IMK2_R.eq.1) then
            IBOFF = 1
          else
            IBOFF = IBOFF + NBLK_MS2_q(IMK2_R-1)
          end if
*
          if(relQC)then 
            NBL_C = NBLK_MS2_q(IMK2_R)
          else
            NBL_C = NBLOCK
          end if

          call REWINE(LUSCR,-1)
          call COPNBLKD(LUC,LUSCR,C,NBL_C,0,LBLK)
          call ITODS(-1,1,-1,LUSCR)

          NAEL = NELEC_q(IACTP)
          NBEL = NELEC_q(IBCTP)
          write(6,*) 'nael, nbel ',nael,nbel
          allocate(astring(NAEL,MXNSTR_q))
          allocate(bstring(NBEL,MXNSTR_q))

*
          NOCTPA   = NOCTYP_q(IACTP)
          NOCTPB   = NOCTYP_q(IBCTP)
*
          IOCTPA   = IBSPGPFTP_q(IACTP)
          IOCTPB   = IBSPGPFTP_q(IBCTP)
*
          IASPGPTP = IACTP
          IBSPGPTP = IBCTP
*
*. Number of occupation classes
          IWAY = 1
          NEL = NAEL + NBEL
          
          if(relQC)then 
          write(6,*) 'nael, nbel before 0',nael,nbel
          call OCCLS_REL(IWAY,NOCCLS,IOCCLS,NEL,NGAS_q,
     &                  IGSOCC_q(1,1),IGSOCC_q(1,2))
          write(6,*) 'nael, nbel before 1',nael,nbel

          else
          write(6,*) 'nael, nbel before 0',nael,nbel
          call OCCLS(IWAY,NOCCLS,IOCCLS,NEL,NGAS_q,
     &           IGSOCC_q(1,1),IGSOCC_q(1,2),
     &           0,0)
          write(6,*) 'nael, nbel before 1',nael,nbel
          end if 
*
*. and then the occupation classes
          allocate(occupation_class(NGAS_q,NOCCLS))
          allocate(mat1(10,NOCCLS))
          allocate(mat2(10,NOCCLS))
          IWAY = 2
          
          if(relQC)then 
          write(6,*) 'nael, nbel before 2',nael,nbel
          call OCCLS_REL(IWAY,NOCCLS,occupation_class,NEL,NGAS_q,
     &                  IGSOCC_q(1,1),IGSOCC_q(1,2))
          write(6,*) 'nael, nbel before 3',nael,nbel

          else
          write(6,*) 'nael, nbel before 2',nael,nbel
          call OCCLS(IWAY,NOCCLS,occupation_class,NEL,NGAS_q,
     &           IGSOCC_q(1,1),IGSOCC_q(1,2),
     &           0,0)
          write(6,*) 'nael, nbel before 3',nael,nbel
          end if 
*

          call get_Q_correction(C,LUSCR,NBL_C,
     &                          WORK(KNSTSO_q(IACTP)),
     &                          WORK(KNSTSO2_q(IBCTP)),
     &                     NOCTPA,NOCTPB,MXPNGAS_q,IOCTPA,IOCTPB,
     &                     NBLOCK,IBLOCK(1,IBOFF),
     &                     THRES,MAXTRM,NAEL,NBEL,
     &                     astring, bstring,
     &                     IBLTP,NSMST,
     &                     mat1, mat2,
     &                     NELFSPGP_q,
     &                     NOCCLS,NGAS_q,occupation_class,
     &                     IASPGPTP,IBSPGPTP,
     &                     ICISTR,
     &                     IPRWFC, ref_wavefunction, 
     &                     ref_wavefunction_coeff,nref,
     &                     ref_wavefunction_coeff_i,NACTEL_q,
     &                     NORB_q,relQC)
*
          deallocate(mat2)
          deallocate(mat1)
!         NAEL/NBEL may be zero and therefore no allocation of string arrays...
          if(allocated(bstring)) deallocate(bstring)
          if(allocated(astring)) deallocate(astring)
!          deallocate(bstring)
!          deallocate(astring)
          deallocate(occupation_class)
        end do
      end do

      !nref = 1       ! for HF nref = 1 TODO: input parameter!!!
      do i = 1, nref
      
      write(6,'(A)') '                                       '
      write(6,'(A)') ' Reference wavefunction & coefficient: ' 
      write(6,'(A)') ' ===================================== '

      print *, "alpha ref",ref_wavefunction(1,1)
      print *, "beta ref",ref_wavefunction(1,2)

      print *,"coefficient of reference wf #",i," == ",
     &         ref_wavefunction_coeff(i)
      print *,"coefficient of HF reference wf             ==",
     &         ref_wavefunction_coeff_i(1)
      end do

      ! for +Q correction
      open(file="energies.CI",unit=10,status="old",
     &     form="unformatted",access="sequential")

      read(10) nstates_ci,reference_energy

      allocate(e_corr(nstates_ci))
      allocate(e_tot(nstates_ci))

      read(10) (e_tot(i),i=1,nstates_ci)
      close(10,status="keep")

      do i = 1, nstates_ci
        e_corr(i) = e_tot(i) - reference_energy
        
        call E_Q(e_corr,ref_wavefunction_coeff,nref,
     &           ref_wavefunction_coeff_i,NACTEL_q,Edc_1,
     &           Edrc_1, Epc_1, Emc_1, Edc_2, Edrc_2,
     &           Epc_2, Emc_2)

      ! output of the q-correction       
        write(6,"(/Af10.6)")  " Energie corrections: "
        write(6,'(A)')        " ==================== "

        write(6,"(/A,f10.6)") "Edc_1:"
        write(6,"(A,1f10.6,/A,1f10.6)") " original ecorr = ",
     &        e_corr(i), " +Q correction= ",Edc_1
        write(6,"(/A,1f10.6,A,1f10.6,A,1f10.6,/)"), 
     &        " final energies: original ", e_tot(i), ", +Q corrected ="
     &        ,e_tot(i) + Edc_1
 

        write(6,"(/A,f10.6)") "Edrc_1:"
        write(6,"(A,1f10.6,/A,1f10.6)") " original ecorr = ",
     &        e_corr(i), " +Q correction = ",Edrc_1
        write(6,"(/A,1f10.6,A,1f10.6,A,1f10.6,/)"), 
     &        " final energies: original ", e_tot(i), ", +Q corrected ="
     &        ,e_tot(i) + Edrc_1


        write(6,"(/A,f10.6)") "Epc_1:"
        write(6,"(A,1f10.6,/A,1f10.6)") " original ecorr = ",
     &        e_corr(i), " +Q correction = ",Epc_1
        write(6,"(/A,1f10.6,A,1f10.6,A,1f10.6,/)"), 
     &        " final energies: original ", e_tot(i), ", +Q corrected ="
     &        ,e_tot(i) + Epc_1


        write(6,"(/A,f10.6)") "Emc_1:"
        write(6,"(A,1f10.6,/A,1f10.6)") " original ecorr = ",
     &        e_corr(i), " +Q correction = ",Emc_1
        write(6,"(/A,1f10.6,A,1f10.6,A,1f10.6,/)"), 
     &        " final energies: original ", e_tot(i), ", +Q corrected ="
     &        ,e_tot(i) + Emc_1


        write(6,"(/A,f10.6)") "Edc_2:"
        write(6,"(A,1f10.6,/A,1f10.6)") " original ecorr = ",
     &        e_corr(i), " +Q correction= ",Edc_2
        write(6,"(/A,1f10.6,A,1f10.6,A,1f10.6,/)"), 
     &        " final energies: original ", e_tot(i), ", +Q corrected ="
     &        ,e_tot(i) + Edc_2
 

        write(6,"(/A,f10.6)") "Edrc_2:"
        write(6,"(A,1f10.6,/A,1f10.6)") " original ecorr = ",
     &        e_corr(i), " +Q correction = ",Edrc_2
        write(6,"(/A,1f10.6,A,1f10.6,A,1f10.6,/)"), 
     &        " final energies: original ", e_tot(i), ", +Q corrected ="
     &        ,e_tot(i) + Edrc_2


        write(6,"(/A,f10.6)") "Epc_2:"
        write(6,"(A,1f10.6,/A,1f10.6)") " original ecorr = ",
     &        e_corr(i), " +Q correction = ",Epc_2
        write(6,"(/A,1f10.6,A,1f10.6,A,1f10.6,/)"), 
     &        " final energies: original ", e_tot(i), ", +Q corrected ="
     &        ,e_tot(i) + Epc_2


        write(6,"(/A,f10.6)") "Emc_2:"
        write(6,"(A,1f10.6,/A,1f10.6)") " original ecorr = ",
     &        e_corr(i), " +Q correction = ",Emc_2
        write(6,"(/A,1f10.6,A,1f10.6,A,1f10.6,/)"), 
     &        " final energies: original ", e_tot(i), ", +Q corrected ="
     &        ,e_tot(i) + Emc_2

      end do


      deallocate(e_corr)
      deallocate(e_tot)

      deallocate(ref_wavefunction)
      deallocate(ref_wavefunction_coeff)
      deallocate(ref_wavefunction_coeff_i)

!
 
      call QEXIT('QCORR')
*
      end subroutine Q_correction

***********************************************************************

      subroutine get_Q_correction(C,LUC,NBL_C,
     &                            NSSOA,
     &                            NSSOB,
     &                       NOCTPA,NOCTPB,MXPNGAS,IOCTPA,IOCTPB,
     &                       NBLOCK,IBLOCK,
     &                       THRES,MAXTRM,NAEL,NBEL,
     &                       IASTR,IBSTR,IBLTP,NSMST,
     &                       NCPMT,WCPMT,NELFSPGP_q,NOCCLS,
     &                       NGAS_q,IOCCLS,
     &                       IASPGPTP,IBSPGPTP,
     &                       ICISTR,IPRNT,
     &                       ref_wavefunction,
     &                       ref_wavefunction_coeff,nref,
     &                       ref_wavefunction_coeff_i,N,
     &                       NORB,relQC)
*
***********************************************************************
*
      implicit real*8(A-H,O-Z)

      dimension C(*)
*. General input
      dimension NSSOA(NSMST,*), NSSOB(NSMST,*)
      dimension IASTR(NAEL,*),IBSTR(NBEL,*)
      dimension IBLTP(*)
      dimension NELFSPGP_q(MXPNGAS,*)
      dimension IOCCLS(NGAS_q,NOCCLS)
      logical, intent(in) :: relQC
      integer, intent(in) :: norb
*. Specific input
      dimension IBLOCK(8,NBL_C)
      character(len=200), intent(in)  :: ref_wavefunction(10,2)

*. Output
      dimension NCPMT(10,NOCCLS)
      dimension WCPMT(10,NOCCLS)
      real(8) , intent(inout) :: ref_wavefunction_coeff(10)
      real(8) , intent(inout) :: ref_wavefunction_coeff_i(10)       ! coeff von der ersten ci_iteration!!!
     
      character (len=4),  allocatable   :: conf_string_a(:)
      character (len=4),  allocatable   :: conf_string_b(:)
      character (len=10), allocatable   :: stringa(:)
      character (len=10), allocatable   :: stringb(:)
      character (len=200), allocatable  :: ref_wavefunction_tmp(:,:)
      logical fndlab, ex, ex_open

      allocate(ref_wavefunction_tmp(10,2))

*
      NTESTL = 0000
      NTEST = max(IPRNT,NTESTL)
*
      IDET      = 0

      if( ICISTR .ge. 2 ) call REWINE(LUC,-1)
      
      IDET     = 0
*
      write(6,*) "nbl_c is",nbl_c
*
      do JBLOCK = 1, NBL_C
*
         IATP = IBLOCK(1,JBLOCK)
         IBTP = IBLOCK(2,JBLOCK)
         IASM = IBLOCK(3,JBLOCK)
         IBSM = IBLOCK(4,JBLOCK)


*. Obtain alpha strings of sym IASM and type IATP
         if(relQC)then 
         call GETSTR_TOTSM_SPGP_REL(IASPGPTP,IATP,IASM,NAEL,
     &                              NASTR1,IASTR)

         else

         IDUM = 0
         call GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
         end if

*. Obtain beta strings of sym IBSM and type IBTP
         if(relQC)then 
         call GETSTR_TOTSM_SPGP_REL(IBSPGPTP,IBTP,IBSM,NBEL,
     &                              NBSTR1,IBSTR)

         else

         IDUM = 0
         call GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                           NORB,0,IDUM,IDUM)
         end if
*
         if(IBLTP(IASM).eq.2) then
            IRESTR_ANA = 1
         else
            IRESTR_ANA = 0
         end if
*
         NIA = NSSOA(IASM,IATP)
         NIB = NSSOB(IBSM,IBTP)
         IMZERO = 0
         write(6,*) "block jblock nia,nib",jblock,nia,nib
         if (ICISTR.ge.2) then
*. Read in a Type-Type-symmetry block
*? Perhaps use GSTTBL here?
            call IFRMDS(IDET,1,-1,LUC)
            call FRMDSC(C,IDET,-1,LUC,IMZERO,IAMPACK)
            IDET = 0
         end if
        write(6,*) ' IMZERO in ciq..ciq..  ',IMZERO
         if(IMZERO.ne.1) then
*
            IBBAS = 1
            IABAS = 1
*
            do IB = IBBAS,IBBAS+NIB-1                                  ! loop über beta strings
               if(IRESTR_ANA.eq.1.and.IATP.eq.IBTP) then
                  MINIA = IB - IBBAS + IABAS
               else
                  MINIA = IABAS
               end if
               do  IA = MINIA,IABAS+NIA-1                              ! loop über alpha strings 
*
                    IDET = IDET + 1
                   
                    write(6,*) 'idet, c',idet,c(idet),nael,nbel
                    write(6,"(60i4)")          ! save strings into array
     &              (iastr(iel,ia),iel=1,nael)                         
                    write(6,"(60i4)") 
     &              (ibstr(iel,ib),iel=1,nbel)
                    write(ref_wavefunction_tmp(1,1),"(60i4)")          ! save strings into array
     &              (iastr(iel,ia),iel=1,nael)                         
                    write(ref_wavefunction_tmp(1,2),"(60i4)") 
     &              (ibstr(iel,ib),iel=1,nbel)

                    ilena   = len(trim(ref_wavefunction_tmp(1,1)))      
                    ilenb   = len(trim(ref_wavefunction_tmp(1,2)))      
                                                                        
                    iref_tmp = 0

                    do 

                      iref_tmp = iref_tmp + 1                           ! exit if iref_tmp > nref
                      if(iref_tmp > nref) exit

                      ilena2  = len(trim(ref_wavefunction(iref_tmp,1)))
                      ilenb2  = len(trim(ref_wavefunction(iref_tmp,2)))

                    if(ref_wavefunction_tmp(1,1)(1:ilena) ==            ! compare alpha strings 
     &                 ref_wavefunction(iref_tmp,1)(1:ilena2) .and.
     &                 ref_wavefunction_tmp(1,2)(1:ilenb) ==            ! compare beta strings 
     &                 ref_wavefunction(iref_tmp,2)(1:ilenb2)
     &                )then 
                       ref_wavefunction_coeff(iref_tmp) = C(IDET)       ! if true, save coeff
   
                    exit
                    end if
                    end do


               end do
*     ^ End of loop over alpha strings
            end do
*     ^ End of loop over beta strings
         end if
*     ^ End of if statement for nonvanishing blocks
      end do
*     ^ End of loop over blocks

      end subroutine get_Q_correction

***********************************************************************

      subroutine E_Q(e_corr,ref_wavefunction_coeff,nref,
     &               ref_wavefunction_coeff_i,N,Edc_1,
     &               Edrc_1, Epc_1, Emc_1, Edc_2, Edrc_2,
     &               Epc_2, Emc_2)
*
*
* Davidson Correction
*
* Equations from: Chemical Reviews 2012, 112 page 216-218
*
* Edc  = eq. 34 p. 116
* Edrc = eq. 36 p. 116
* Epc  = eq. 38 p. 116
* Emc  = eq. 41 p. 116
*
* c01  = eq. 43 p. 117
* c02  = eq. 44 p. 117
*
***********************************************************************
*
      implicit none

      real*8,  intent(in)  :: e_corr
      real*8,  intent(in)  :: ref_wavefunction_coeff(nref)
      real*8,  intent(in)  :: ref_wavefunction_coeff_i(nref)
      real*8,  intent(out) :: Edc_1, Edrc_1, Epc_1, Emc_1
      real*8,  intent(out) :: Edc_2, Edrc_2, Epc_2, Emc_2
      real*8,  intent(in)  :: N
      integer, intent(in)  :: nref

      real*8  :: c01    = 0      ! c01    = sum ( cp ) **2
      real*8  :: c02    = 0      ! c02    = (sum ( cref * cp )) **2
      real*8  :: c02_2  = 0      ! c02_2  = c02 * c02, entspricht c01
      real*8  :: teta_1 = 0      ! teta_1 = arccos( sqrt(c01) )
      real*8  :: teta_2 = 0      ! teta_2 = arccos( sqrt(c02_2))
      integer :: i      = 0
      

      do i = 1, nref
      c01 = c01 + ref_wavefunction_coeff(i) ** 2
      end do

      teta_1 = acos(sqrt(c01))

      do i = 1, nref
      c02 = c02 + ref_wavefunction_coeff(i)*ref_wavefunction_coeff_i(i)
      end do

      c02_2  = c02 * c02
      teta_2 = acos(sqrt(c02_2))

      write(6,*) "c01",         c01
      write(6,*) "sqrt (c01)",  sqrt(c01)
      write(6,*) "c02",         c02
      write(6,*) "c02_2",       c02_2
     
      Edc_1  = (1 - c01) * e_corr            
      Edrc_1 = ((1-c01) / c01) * e_corr
      Epc_1  = (sqrt(N **2 + 2*N * (tan(2*teta_1))**2) - N) / 
     &         (2*(1/(2*teta_1)) - 1) * e_corr

      if(N > 3.0d0) Emc_1 = ((N - 2) * (N - 3) / (N * (N - 1))) *
     &                      ((1 - c01 **2) / c01 **2) * e_corr  

      Edc_2  = (1 - c02_2) * e_corr
      Edrc_2 = ((1-c02_2) / c02_2) * e_corr
      Epc_2  = (sqrt(N **2 + 2*N * (tan(2*teta_2))**2) - N) / 
     &         (2*(1/(2*teta_2)) - 1) * e_corr

      if(N > 3.0d0) Emc_2 = ((N - 2) * (N - 3) / (N * (N - 1))) *
     &                      ((1 - c02_2 **2) / c02_2 **2) * e_corr  
      
      end subroutine E_Q



