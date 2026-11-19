module module_cmo_handler
    use module_cmo_variables
    implicit none

    integer :: idims(12) ! CMO data dimensions
    real(8) :: toterg ! total energy, not used in this module but required by dircmo

    public :: ivo_cmo_read, ivo_cmo_write
contains
    subroutine ivo_cmo_read
        use dircmo, only: reacmo_new
        implicit none
        integer :: ncoef, neig

        idims = 0
        ! first, read CMO data dimensions
        call reacmo_new(idims=idims)
        A = idims(2) ! A is nfsym
        B = idims(1) ! B is nz
        ncoef = idims(3) ! ncoef is the number of coefficients, see reacmo_new in dircmo.F90
        neig = idims(5) ! neig is the number of orbitals, see reacmo_new in dircmo.F90
        positronic_mo = (/idims(7), idims(10)/)
        electronic_mo = (/idims(8), idims(11)/)
        mo = positronic_mo + electronic_mo
        basis_ao = (/idims(9), idims(12)/)
        basis_all = (positronic_mo + electronic_mo)*basis_ao
        allocate (BUF(ncoef), eval(neig), syminfo(neig))
        call reacmo_new(cmo=BUF, eig=eval, ibeig=syminfo, toterg=toterg)
    end subroutine ivo_cmo_read

    subroutine ivo_cmo_write
        use dircmo, only: wricmo
        use module_global_variables, only: rank
        implicit none
        integer :: iunit ! dummy io unit

        ! only master process writes the CMO data to checkpoint file
        if (rank == 0) then
            call wricmo(iunit=iunit, cmo=BUF, eig=eval, ibeig=syminfo, toterg=toterg, &
                        title='CMO data written by ivo_cmo_write', idims=idims)
        end if
    end subroutine ivo_cmo_write
end module module_cmo_handler
