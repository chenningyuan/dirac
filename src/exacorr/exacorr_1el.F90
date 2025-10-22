!Datatype for keeping effecive one-electron datatypes. Should likely not be a separate module after deprecating exatensor.
       module exacorr_1el

        implicit none
        private

        integer, parameter, private:: REALD=8 ! default real for this module

!       One electron integrals
        type, public  :: one_el_t
         integer     :: n_spinor=-1           ! number of spinors, used for consistency checks
         integer     :: n_prop=0              ! number of property matrices that were read
         real(8)     :: e_core                ! energy of the core electrons plus nuclear repulsion energy
         complex(8), pointer :: h_core(:,:)   ! one-body part of Hamiltonian (kinetic energy + nuclear attraction + core shielding)
         complex(8), pointer :: h_prop(:,:,:) ! property integrals to enable adding a finite field perturbation
         character*8,pointer :: property_labels(:)  ! strings that should correspond to the label as written to MDPROP (or other interface files)
        end type one_el_t

       end module exacorr_1el
