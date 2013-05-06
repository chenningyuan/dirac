MODULE adc_fano

CONTAINS

  SUBROUTINE fanoadc(ckks,bufhp,ladc,krep,itapadc,iw,  &
                     eajl,oooo,vovo,mxno,nr_lanc_leigv)
!
!  Purpose: Calculate decay widths using the FanoADC method.
!           This is the controlling module. All routines are to be found
!           in adc_fano_routines.
!           This module is supposed to be tidied up. No other routines are
!           to be added.
!
!  Author:  Elke Fasshauer

    use adc_cfg
    use adc_fano_real_routines
    use adc_fano_diag
    use adc_fano_routines
    use memory_allocator

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
    REAL(8), INTENT(IN), DIMENSION(nvoot(krep)*nvoot(krep)) :: bufhp
    INTEGER, INTENT(IN)   :: ladc
    INTEGER, INTENT(IN)   :: itapadc
    INTEGER, INTENT(IN)   :: iw
    REAL(8), INTENT(IN)   :: eajl(*)
    REAL(8), INTENT(IN)   :: oooo(*), vovo(*)
    INTEGER, INTENT(IN)   :: mxno
    INTEGER, INTENT(IN)   :: nr_lanc_leigv !number of requested long eigenvectors
    

!
!  Data dictionary: Local Variables (and to be shared with subroutines)
!
    INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: icra
    INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: bigicra
    INTEGER                               :: nr_in, nr_fin
    INTEGER                               :: nok, najl
    INTEGER, ALLOCATABLE, DIMENSION(:)    :: finalpos
    INTEGER                               :: recl_condat_s
    CHARACTER(8)                          :: fn_conf='FANCNF'

!
!  Code
!
    WRITE(iw,*)
    WRITE(iw,*) '*********************************************'
    WRITE(iw,*) '*                 FANOADC                   *'
    WRITE(iw,*) '*       written by Elke Fasshauer           *'
    WRITE(iw,*) '*********************************************'
    WRITE(iw,*) 

    IF (nr_lanc_leigv == 0) THEN
      CALL quit('FanoADC requires long eigenvectors')
    END IF

    ! Create lookup table bigicra
    CALL alloc(bigicra,9,NVOOT(krep), id='FanoADC: bigicra')
    CALL make_bigicra(krep,bigicra,iw)
    CALL alloc(icra,6,nvoot(krep), id='FanoADC: icra - integer')


    nok  = NO(krep)
!    najl = NOOVT(krep)
    
    nr_in   = det_nr_configs(nok+1,ladc,bigicra,krep,1,iw)
    nr_fin  = det_nr_configs(nok+1,ladc,bigicra,krep,2,iw)
    fin_max = nr_fin
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

    CALL dealloc(icra)
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
    CALL FULLDIAR2(itapadc,fn_finmat,fano_intbuf,finnbufs,fin_max,krep)

    ! Create the file holding the configurations of the 1h and 2h1p configurations
    CALL fano_wcondat_s(itapfano+6,fn_conf,recl_condat_s,krep,ladc-fin_max, &
                        finalpos,fin_max,iw)

    CALL dealloc(finalpos)

  END SUBROUTINE fanoadc

END MODULE adc_fano
