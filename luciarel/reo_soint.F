*  Reorder SO integrals according to matrices INDREO
*
      subroutine reo_soint(SOIN,HELPMAT,HELPMAT2,IDIM,
     &                     INDREO,NTRIAIN,NMAXBAS,IPRINT)
*
      implicit real*8 (A-H,O-Z)
*
      dimension SOIN(NTRIAIN),HELPMAT(IDIM,IDIM),
     &          HELPMAT2(IDIM,IDIM)
      dimension INDREO(NMAXBAS)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
      call unfold(HELPMAT,SOIN,IDIM,1,NTRIAIN,1)
      call matsort(HELPMAT,HELPMAT2,INDREO,INDREO,IDIM,IDIM)
*
      if (NTEST.ge.3) then
        write(6,*) 'Checking reordered SO integral matrix: '
        call wrtmat(HELPMAT,IDIM,IDIM,IDIM,IDIM)
      end if
*
      call fold(SOIN,HELPMAT,IDIM,1,NTRIAIN)
*
      return
      end
