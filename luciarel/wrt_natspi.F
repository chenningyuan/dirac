*  Write natural (optimized) spinors to files
*   LUMOUTR (real part)
*   LUMOUTI (imaginary part)
*
      subroutine wrt_natspi(ONATSPI,IDIM,LUMOUTR,
     &                  LUMOUTI,ISPSM)
*
      implicit real*8 (A-H,O-Z)
*
      dimension ONATSPI(2*IDIM,IDIM)
*
 100  format(4F18.14)
*
      write(6,*)
      write(6,*) 'Writing real part of natural spinors'
      write(6,*) 'to file NATSPIR.'
      write(LUMOUTR,'(A16,I2)') 'SPINOR SYMMETRY ',ISPSM
      write(LUMOUTI,'(A16,I2)') 'SPINOR SYMMETRY ',ISPSM
      do ICOL=1,IDIM,1
        write(LUMOUTR,'(A14,I4)') 'Spinor number ',ICOL
        write(LUMOUTR,100) (ONATSPI(I,ICOL),I=1,IDIM,1)
        write(LUMOUTI,'(A14,I4)') 'Spinor number ',ICOL
        write(LUMOUTI,100) (ONATSPI(I,ICOL),I=IDIM+1,2*IDIM,1)
      end do
*
      return
      end
