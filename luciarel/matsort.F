*  Sort matrix elements according to new index pattern
*  Output matrix is resorted input matrix
*
      subroutine matsort(ZMAT,ZMATHELP,IINOUT,JINOUT,IDIM,JDIM)
*
      implicit real*8 (A-H,O-Z)
*
      dimension IINOUT(IDIM),JINOUT(JDIM)
      dimension ZMAT(IDIM,JDIM),ZMATHELP(IDIM,JDIM)
*
* IINOUT(INPUT INDEX) = OUTPUT INDEX
* JINOUT(INPUT INDEX) = OUTPUT INDEX
*
      do I=1,IDIM,1
        do J=1,JDIM,1
          if (IINOUT(I).gt.IDIM.or.JINOUT(J).gt.JDIM.
     &        or.IINOUT(I).lt.1.or.JINOUT(J).lt.1) then
            write(6,*) 'Problem in matsort'
            write(6,'(A7,I6,A3,I6)') 'IINOUT(',I,') =',IINOUT(I)
            write(6,'(A7,I6,A3,I6)') 'JINOUT(',J,') =',JINOUT(J)
            Call Abend2( 'Quitting.' )
          end if
          ZMATHELP(IINOUT(I),JINOUT(J)) = ZMAT(I,J)
        end do
      end do
*
      call copvec(ZMATHELP,ZMAT,IDIM*JDIM)
*
      return
      end
