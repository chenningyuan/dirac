MODULE adc_fano

CONTAINS

  SUBROUTINE fanoadc(ckks,bufhp,ladc,krep,itapadc,iw,  &
                     eajl,oooo,vovo,mxno)
!
!  Purpose: Calculate decay widths using the FanoADC method.
!           This is the controlling module. All routines are to be found
!           in adc_fano_routines.
!           This module is supposed to be tidied up. No other routines are
!           to be added.
!
!  Author:  Elke Fasshauer

    use adc_cfg
    use adc_mat
    use adc_fano_real_routines
    use adc_fano_diag
    use adc_fano_routines
    use memory_allocator
    use adc_fano_exchange
    use adc_fano_matmul

    IMPLICIT NONE

!
!  Common Blocks
!
#include "../relccsd/symm.inc"

!
!  Data dictionary: Calling Variables
!
    INTEGER, INTENT(IN)   :: krep
    REAL(8), INTENT(IN), DIMENSION(no(krep)*no(krep))       :: ckks
    REAL(8), INTENT(IN), DIMENSION(no(krep)*nvoot(krep)) :: bufhp
    INTEGER, INTENT(IN)   :: ladc
    INTEGER, INTENT(IN)   :: itapadc
    INTEGER, INTENT(IN)   :: iw
    REAL(8), INTENT(IN)   :: eajl(*)
    REAL(8), INTENT(IN)   :: oooo(*), vovo(*)
    INTEGER, INTENT(IN)   :: mxno

!
!  Data dictionary: Local Variables (and to be shared with subroutines)
!
    INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: bigicra
    INTEGER                               :: nr_in, nr_fin
    INTEGER                               :: nok, najl
    INTEGER, ALLOCATABLE, DIMENSION(:)    :: finalpos
    INTEGER                               :: recl_condat_s
    INTEGER                               :: recl_fin
    CHARACTER(6)                          :: fn_conf='FANCNF'
    CHARACTER(6)                          :: fn_finconf='FINCNF'
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: buflanc
    INTEGER                               :: kfree, lfree, kfreebase
!    INTEGER                               :: ierror = 0
    INTEGER                               :: idpt
    INTEGER                               :: ixx

    REAL(8), ALLOCATABLE, DIMENSION(:)    :: e2h1p, t_mom

!
!  Interfaces
!
    INTERFACE
      SUBROUTINE diag_lanc(cc,kfree,lfree,iw,doincore)
        REAL(8), DIMENSION(*)  :: cc
        INTEGER                :: kfree, lfree, iw
        LOGICAL                :: doincore
      END SUBROUTINE diag_lanc
    END INTERFACE

!
!  Code
!
    WRITE(iw,*)
    WRITE(iw,*) '*********************************************'
    WRITE(iw,*) '*                 FANOADC                   *'
    WRITE(iw,*) '*       written by Elke Fasshauer           *'
    WRITE(iw,*) '*********************************************'
    WRITE(iw,*) 

    nok  = NO(krep)
    najl = NVOOT(krep)

!    ! Check physical possibility of initial state choice
    IF (reladc_fano_inrelsp > nok) THEN
      CALL QUIT('Unphysical choice of initial state')
    END IF

    ! Create lookup table bigicra
    CALL alloc(bigicra,9,NVOOT(krep), id='FanoADC: bigicra')
    CALL make_bigicra(krep,bigicra,iw)

    nr_in   = det_nr_configs(nok+1,ladc,bigicra,krep,1,iw)
    nr_fin  = det_nr_configs(nok+1,ladc,bigicra,krep,2,iw)
    fin_max = nr_fin

    IF (fin_max == 0) THEN
      CALL QUIT('No valid Fano final states given')
    END IF

    WRITE(iw,*) 'Number of initial state configurations: ', nr_in
    WRITE(iw,*) 'Number of  final  state configurations: ', nr_fin

    ! Create sorted h/2h1p block
    CALL fano_write_hhp(ckks,bufhp,ladc,krep,bigicra,   &
                             itapadc,iw)

    ! Create an array of final state positions in the 2h1p part
    CALL alloc(finalpos,fin_max, id='FanoADC: Array of final positions')
    CALL make_finalpos(krep,bigicra,finalpos,fin_max,iw)

    ! Create the sorted 2h1p/2h1p block
    CALL fano_makehphp(ckks,eajl,oooo,vovo,bigicra,ladc,krep,  &
                       finalpos,iw,mxno)

    CALL dealloc(bigicra)

    WRITE(iw,*)
    WRITE(iw,*) 'Diagonal elements of the sorted ADC matrix written to        : ', &
                 fn_adcdiag
    WRITE(iw,*) 'Non-diagonal elements of the sorted ADC matrix written to    : ', &
                 fn_adcmat
    WRITE(iw,*) 'Diagonal elements of initial state matrix written to         : ', &
                 fn_indiag
    WRITE(iw,*) 'Non-diagonal elements of initial state matrix written to     : ', &
                 fn_inmat
    WRITE(iw,*) 'Diag. and non-diag. elements of final state matrix written to: ', &
                 fn_finmat
    WRITE(iw,*)

    

    WRITE(iw,*) 'The matrix of final states is diag. by full diagonalization'
    !reladc_md_isfano   = .true.
    CALL FULLDIAR2(itapadc,fn_finmat,fano_intbuf,finnbufs,fin_max,krep)

    ! Create the file holding the configurations of the 1h and 2h1p configurations
    !CALL fano_wcondat_fin_s(itapfano+6,fn_finconf,recl_fin,krep,fin_max, &
    !                    finalpos,fin_max,iw)
    CALL fano_wcondat_s(itapfano+6,fn_conf,recl_condat_s,krep,ladc-fin_max, &
                        finalpos,fin_max,iw)
    !CALL WCONDAT_S(itapfano+6,fn_conf,recl_condat_s,krep,ladc)

    CALL dealloc(finalpos)

!!
!! This is only for testing purposes
!!
!    ! Write all characteristics of the matrix to be diagonalised
!    ! to adc_mat. We are using the lanczos diagonaliser.
!    reladc_md_iobase   = itapfano+6
!    reladc_md_ionizl   = 1
!    reladc_md_ioldnew  = 2 !new, diagonal is not written to extra file
!    reladc_md_intbuf   = fano_intbuf
!    reladc_md_desrep   = krep
!    reladc_md_rcw      = 1 !is real, for complex 2
!    reladc_md_lnzitr   = reladc_sipiter
!    reladc_md_matdim   = fin_max
!    reladc_md_irecl    = recl_fin
!    reladc_md_nmain    = nok
!    reladc_md_nbufs    = finnbufs
!    reladc_md_neigenv  = reladc_fano_eigv
!    reladc_md_fileadc  = fn_finmat
!    reladc_md_filecnf  = fn_finconf
!    reladc_md_nmspec   = 'FSPEF'
!
!    
!    CALL alloc(buflanc,reladc_lancmem, id='FanoADC Lanczos memory')!, stat=ierror)
!    
!    kfree     = 1
!    lfree     = reladc_lancmem
!    kfreebase = 1
!
!    CALL MEMOPT('REAL',IDPT,reladc_lancmem-10,buflanc,  &
!     &                KFREE,LFREE,.true.)
!    CALL MEMREL('IDPT',buflanc,KFREEBASE,KFREEBASE,KFREE,LFREE)
!
!    CALL DIAG_LANC(buflanc,kfree,lfree,iw,.false.)
!    ! For the next runs Lanczos is again in normal mode
!    CALL dealloc(buflanc)

!
! This is the actual code
!
    reladc_md_isfano   = .true.
    ! Write all characteristics of the matrix to be diagonalised
    ! to adc_mat. We are using the lanczos diagonaliser.
    reladc_md_iobase   = itapfano+6
    reladc_md_ionizl   = 1
    reladc_md_ioldnew  = 1 !old, diagonal is written to extra file
    reladc_md_intbuf   = fano_intbuf
    reladc_md_desrep   = krep
    reladc_md_rcw      = 1 !is real, for complex 2
    reladc_md_lnzitr   = reladc_sipiter
    reladc_md_matdim   = ladc-fin_max
    !reladc_md_matdim   = ladc
    reladc_md_irecl    = recl_condat_s
    reladc_md_nmain    = nok
    reladc_md_nbufs    = innbufs
    !reladc_md_nbufs    = adcnbufs
    reladc_md_neigenv  = reladc_fano_eigv
    reladc_md_fileadc  = fn_inmat
    !reladc_md_fileadc  = fn_adcmat
    reladc_md_filediag = fn_indiag
    !reladc_md_filediag = fn_adcdiag
    reladc_md_filecnf  = fn_conf
    reladc_md_nmspec   = 'FSPEC'

    
!    CALL FULLDIAR(itapfano,fn_indiag,fn_inmat,fano_intbuf,ladc-fin_max,krep)
!    CALL FULLDIAR(itapfano,fn_adcdiag,fn_adcmat,fano_intbuf,ladc,krep)

    CALL alloc(buflanc,reladc_lancmem, id='FanoADC Lanczos memory')!, stat=ierror)
    
    kfree     = 1
    lfree     = reladc_lancmem
    kfreebase = 1

    CALL MEMOPT('REAL',IDPT,reladc_lancmem-10,buflanc,  &
     &                KFREE,LFREE,.true.)
    CALL MEMREL('IDPT',buflanc,KFREEBASE,KFREEBASE,KFREE,LFREE)

    !CALL DIAG_LANC(buflanc,kfree,lfree,iw,.false.)
    ! For the next runs Lanczos is again in normal mode
    reladc_md_isfano = .false.
    CALL dealloc(buflanc)

    WRITE(iw,*)
    WRITE(iw,*) '------------------------------------------'
    WRITE(iw,*) 'All vectors and matrices are now present'
    WRITE(iw,*) 'The next step of the Fano run is the    '
    WRITE(iw,*) 'matrix vector multiplication.'
    WRITE(iw,*) '------------------------------------------'
    WRITE(iw,*)

    CALL fano_create_t_mom(fn_adcdiag,fn_adcmat,ladc,nr_in,nr_fin, &
                           adcfh,adcdiagfh,fano_intbuf,e2h1p,t_mom,iw)

  END SUBROUTINE fanoadc

END MODULE adc_fano
