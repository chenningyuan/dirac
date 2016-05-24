 module symmetry_offset

 implicit none 

  public auto_symmetry_offset
  public auto_symmetry_offset_triangular

 private
    public Offset
    type Offset
      integer, dimension(32)          :: oneDirac,oneNonDirac
      integer, dimension(32,32)       :: twoDirac,twoNonDirac
    endtype Offset

  contains

  subroutine auto_symmetry_offset(a,row,column,shift_for_boson_row,shift_for_boson_column)
#include "symm.inc"

! In this subroutine my target is to generate the offset arrays for symmetry.  
   integer,intent(in)  :: row(nrep),column(nrep)
   integer             :: irep, jrep, ijrep
   logical,intent(in)  :: shift_for_boson_row, shift_for_boson_column! a logical variable indicating whether the shift for bosonic irrreps is necessary 
   type(offset),intent(out) :: a
   integer :: shift1,shift2

! 1> case 1: Offset for (2,2) type integrals.
!            here we have dealt with four different type of offsets 
!            oneDirac ->  bra(k)*ket(l) = ket(kl) !multb(*,*,2)
!                      or,ket(k)*bra(l) = bra(kl)
!            twoDirac ->  bra(k)*ket(l) = ket(k,l)!multb(*,*,1) 
!                      or,ket(k)*bra(l) = bra(k,l)

!            oneNonDirac ->  bra(k)*bra(l) = bra(kl) !multb(*,*,2)
!                      or,ket(k)*ket(l) = ket(kl)
!            twoNonDirac ->  bra(k)*bra(l) = bra(k,l)!multb(*,*,1) 
!                      or,ket(k)*ket(l) = ket(k,l)

!  If the irreps are bosonic we have used a shift of nrep (the total number of irreps) i.e, for case1.
!  But, for fermionic irreps the shift is not necessary, that is the case2   

   shift1 = 0
   shift2 = 0
   if(shift_for_boson_row)    shift1=nrep 
   if(shift_for_boson_column) shift2=nrep 

   a%oneDirac(1:nrep) = 0
   a%oneNonDirac(1:nrep) = 0
   a%twoDirac(1:nrep,1:nrep) = 0
   a%twoNonDirac(1:nrep,1:nrep) = 0

    do jrep = 1,nrep
      do irep = 1,nrep
        ijrep = multb(irep+shift1,jrep+shift2,2)
        a%twoDirac(irep,jrep) = a%oneDirac(ijrep) 
        a%oneDirac(ijrep) = a%oneDirac(ijrep)+row(irep)*column(jrep) 

        ijrep = multb(irep+shift1,jrep+shift2,1)
        a%twoNonDirac(irep,jrep) = a%oneNonDirac(ijrep) 
        a%oneNonDirac(ijrep) = a%oneNonDirac(ijrep)+row(irep)*column(jrep) 
      enddo 
    enddo 

   end subroutine

!=============================================================================
  subroutine auto_symmetry_offset_triangular(a,row,column)
#include "symm.inc"

! In this subroutine my target is to generate the offset arrays for symmetry.  
   integer,intent(in)  :: row(nrep),column(nrep)
   integer             :: irep, jrep, ijrep
   type(Offset),intent(out) :: a

   a%oneNonDirac(1:nrep) = 0
   a%twoNonDirac(1:nrep,1:nrep) = 0

    do jrep = 1,nrep
        ijrep = multb(jrep,jrep,1)
        a%twoNonDirac(jrep,jrep) = a%oneNonDirac(ijrep) 
        a%oneNonDirac(ijrep) = a%oneNonDirac(ijrep)+row(jrep)*(column(jrep)-1)/2 

      do irep = jrep+1, nrep
        ijrep = multb(irep,jrep,1)
        a%twoNonDirac(irep,jrep) = a%oneNonDirac(ijrep) 
        a%oneNonDirac(ijrep) = a%oneNonDirac(ijrep)+row(irep)*column(jrep) 
      enddo 
    enddo 

   end subroutine


  end module
