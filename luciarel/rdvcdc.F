      subroutine rdvcdc(VE,SCR,LU,LBLK,IPRT)
*
* Read vector from disc to array
*
* Complex version      Timo Fleig, Dec 04, 2001
*
      include 'implicit.inc'
*. Output
      DIMENSION VE(*),SCR(*)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRT)
*
      CALL REWINE(LU,LBLK)
*
*
*. Loop over blocks of file
*
      if (LBLK.gt.0) then
        LENGTH = LBLK
        call frmdsc(VE,LENGTH,-1,LU,IMZERO,IAMPACK)
        stop 'Take care: Check me first!!.'
      else
        NASS = 0
        IOFF = 1
        do IRI = 1,2,1
 1000     continue
          call ifrmds(LENGTH,1,LBLK,LU)
          if (LENGTH.ge.0) then
            call frmdsc(SCR,LENGTH,LBLK,LU,IMZERO,IAMPACK)
            if (NTEST.ge.10) then
              write(6,*) ' Record read in: '
              call wrtmat(SCR,1,LENGTH,1,LENGTH)
            end if
            call copvec(SCR,VE(IOFF),LENGTH)
            IOFF = IOFF + LENGTH
          end if
          if (LENGTH.ge.0) then
            NASS = NASS + LENGTH
            goto 1000
          end if
        end do
      end if
*
      if (NTEST.ge.5) then
        write(6,*) ' Read vector in rdvcdc (real):  '
        call wrtmat(VE(1),1,NASS/2,1,NASS/2)
        write(6,*) ' Read vector in rdvcdc (imag):  '
        call wrtmat(VE(1+NASS/2),1,NASS/2,1,NASS/2)
      end if
*
      return
      end
