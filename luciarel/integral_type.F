      subroutine integral_type(ISOB,NGAS,ITPSV,LUINGR)
      include 'implicit.inc'
      character*11 CELEC2
      character*13 CELEC
      dimension ISOB(NGAS,4)
      ISUM = 0
      do J=1,4,1
        do I=1,NGAS,1
          ISUM = ISUM + ISOB(I,J)
        end do
      end do
      if (ISUM.eq.2) CELEC = 'One-electron '
      if (ISUM.eq.4) CELEC = 'Two-electron '
      NAC = 0
      NBC = 0
      NAA = 0
      NBA = 0
      do I=1,NGAS,1
        NAC = NAC + ISOB(I,1)
        NBC = NBC + ISOB(I,2)
        NAA = NAA + ISOB(I,3)
        NBA = NBA + ISOB(I,4)
      end do
      if (NAC.eq.NAA.and.NBC.eq.NBA) CELEC2 = 'spin-free  '
      if (NAC.ne.NAA.or.NBC.ne.NBA)  CELEC2 = 'spin-orbit '
      if (ISUM.ne.ITPSV) 
     &     write(LUINGR,'(A,A,A)') CELEC,CELEC2,'integrals.'
      ITPSV = ISUM
      return
      end
