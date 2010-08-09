!
!      Copyright (c) 2008 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the
!      research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.sdu.dk
!
! Contains a module for common matrix and vector operations/procedures.
!
! S. Knecht, 30-03-2010
!

module common_matvec_op
  implicit none

  private
  public lowgen90
  public calc_diff_a_b

contains

!**********************************************************************
  subroutine lowgen90(smat,ndim_s,vmat,nv,ldm,sstol,eigvl,iord,wrtunit)
!**********************************************************************
!
!     Generate Lowdin-type matrix
!     (Lowdin's canonical orthonormalization)
!
!     V_ij = U_ij/SQRT(s_j) where U_ij is the matrix that diagonalizes
!     the overlap matrix S and s_j is eigenvalue j.
!     Linear dependence is removed by eliminating columns of V
!     corresponding to eigenvalues of the overlap matrix below
!     the given threshold.
!
!     lowgen written by T.Saue 1995
!
!     adaptation for F90 by S. Knecht April 2010
!
!**********************************************************************
     real(8), intent(inout) :: smat(ldm,*)
     real(8), intent(inout) :: vmat(ldm,*)
     real(8), intent(inout) :: eigvl(*)
     real(8), intent(in)    :: sstol
     integer, intent(inout) :: nv
     integer, intent(in)    :: ldm
     integer, intent(in)    :: ndim_s
     integer, intent(in)    :: iord
     integer, intent(in)    :: wrtunit
!----------------------------------------------------------------------
     real(8)                :: d1
     real(8)                :: scaling_factor
     integer                :: ierr_qdiag
     integer                :: i
!**********************************************************************
!
       d1         =  1.0D0
       ierr_qdiag =  0
 
#ifdef MOD_LOG
      call matexport_text('% entering lowgen90')
      call matexport_double2('S',smat,ldm,ndim_s,ndim_s)
#endif

!
!     diagonalize overlap matrix
      call qdiag90(1,ndim_s,smat,ldm,ndim_s,eigvl,1,vmat,ldm,ndim_s,ierr_qdiag)
      if(ierr_qdiag /=0 )then
         call quit('*** error in lowgen90: diagonalization of overlap matrix failed. ***')
      end if
!     qdiag gives a reverse value ordering, so fix that here.
!     it would be better though to rewrite the code below instead.. 
      if( iord > 0 ) call order3(vmat,eigvl,ldm,ndim_s,ndim_s,-1)
!
!    ... ensure V_ij = U_ij/SQRT(s_j)
      nv = 0
      do i = 1, ndim_s
        if(eigvl(i)<= sstol) goto 10
        nv  = nv + 1
        scaling_factor = d1/sqrt(eigvl(i))
        call dscal(ndim_s,scaling_factor,vmat(1,i),1)
      end do
 10   continue

#ifdef MOD_LOG
      call matexport_double2('V',vmat,ldm,ns,nv)
      call matexport_text('% leaving lowgen90')
#endif

  end subroutine lowgen90

!**********************************************************************
  subroutine calc_diff_a_b(mat1,      &
                           mat2,      &
                           ilen,      &
                           f1lab,     &
                           f2lab)
!**********************************************************************
!
! calculate the maximum deviation between two matrix elements of matrix a and b 
! which are read from files.
!
!----------------------------------------------------------------------
     real(8), intent(inout) :: mat1(*)
     real(8), intent(inout) :: mat2(*)
     integer, intent(in)    :: ilen
     character  (len=7)     :: f1lab
     character  (len=7)     :: f2lab
!----------------------------------------------------------------------
     real(8)                :: maxdiff
     integer                :: i, j
!**********************************************************************
 
       open(99,file=f1lab,status='unknown',form='unformatted')
       read(99) (mat1(j),j=1,ilen)
       close(99,status="keep")
       open(99,file=f2lab,status='unknown',form='unformatted')
       read(99) (mat2(j),j=1,ilen)
       close(99,status="keep")
       print '(2x,a,a,a,a)', '* calc_diff_a_b: file 1 = ',f1lab,' and file 2 = ',f2lab

       maxdiff = 0.0D0
       do i = 1, ilen
         if(i < 10 ) print '(2x,a,i0,a,e26.18,2x,e26.18)', '* calc_diff_a_b: test print - element i= ',i,'in matrix a and b are',mat1(i),mat2(i)
         maxdiff = max(maxdiff,abs(abs(mat1(i))-abs(mat2(i))))
       end do

       print '(2x,a,e26.18)', '* calc_diff_a_b: largest deviation between two matrix elements in mat a and b: ',maxdiff

  end subroutine calc_diff_a_b

end module
