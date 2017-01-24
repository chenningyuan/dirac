!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PH_TESTER(desrep)
      IMPLICIT INTEGER (A-Z)

      Integer                           :: desrep

#include "../relccsd/symm.inc"   ! provide NREP, MULTB and S-arrays

      Do airep = 1,nrep
        locnum=0
        write(iw,*) 'destination airep:',airep
        DO IREP = 1, NREP
          DO AREP = 1, NREP
            xrep = MULTB(AREP,IREP,2)
            If(xrep.eq.airep) then
               write(*,*) airep,':   G_a:',arep,
     &           '   G_i: ',IREP
               do i=1,no(IREP)
               do a=1,nv(AREP)
                 write(*,*) '    a(absindx),i(absindx):',
     &             io(nrep)+iv(arep)+a,io(irep)+i
                 locnum=locnum+1
               enddo
               enddo
            Endif  !xrep.eq.airep
          ENDDO
        ENDDO
        write(*,*)
        write(*,*) 'Number of valid p/h combinations:',locnum
        write(*,*)
      Enddo    !airep


      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
