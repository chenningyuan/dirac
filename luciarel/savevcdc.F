* Use vector partitioning information to save a vector
* to LUCIAREL file.
*
*    Timo Fleig, March 16, 2002
*
      subroutine savevcdc(VECT,LU,IPRNT)
*
#include "implicit.inc"
#include "mxpdim.inc"
*
#include <civecinf.h>
*
      dimension VECT(*)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
      if (NTEST.ge.5) then
        write(6,*)
        write(6,*) ' +++++++++++++++++++++++++++++++ '
        write(6,*) '  SAVEVCDC called.               '
        write(6,*) '   (Save complex vec. to file.)  '
        write(6,*) ' +++++++++++++++++++++++++++++++ '
        write(6,*) '  File structure information:'
        write(6,*) '  Number of MK values :',NMK2VAL
        JBATABS = 0
        do IMK2 = 1,NMK2VAL,1
          write(6,*) '  MK value no.: ',IMK2
          do JBATS = 1,NBAT_MS2_MC(IMK2)
            JBATABS = JBATABS + 1
            write(6,*) '   Absolute batch no.: ',JBATABS
            write(6,*) '   No. of elements in batch :',
     &                     LEBATS_MC(JBATABS)
          end do
        end do 
      end if
*
      JOFF = 1
      IRIOFF = 0
      LENGTH = 0
      IRILP = 2
      DO IRI = 1,IRILP,1
        JBATABS = 0
        if (NTEST.ge.10) write(6,*) '  REAL/IMAG = ',IRI
        DO IMK2 = 1, NMK2VAL
*
* Loop over batches of sigma for this MK2 projection value
          do JBATS = 1,NBAT_MS2_MC(IMK2)
            JBATABS = JBATABS + 1
*
            if (NTEST.ge.10) then
               WRITE(6,*)
               WRITE(6,*) ' Start of writing of batch ', JBATABS
               write(6,'(2X,A,I3)') ' for MK2 value number ',IMK2
               WRITE(6,*) 
            end if
*
            LS = LEBATS_MC(JBATABS)
            LENGTH = LENGTH + LS
            if (NTEST.ge.10.and.LS.gt.0) then
              WRITE(6,*) ' **********************************'
              WRITE(6,*) ' * Array containing current batch *'
              WRITE(6,*) ' **********************************'
              CALL WRTMAT(VECT(JOFF),1,LS,1,LS)
              JOFF = JOFF + LS
            end if
            KBATSEND = I1BATS_MC(JBATABS) + LBATS_MC(JBATABS)
            do ISBLK = I1BATS_MC(JBATABS),KBATSEND-1,1
              IOFF = IRIOFF + IBATS_MC(6,ISBLK)
              LEN_V = IBATS_MC(8,ISBLK)
              if (NTEST.ge.100) then
                write(6,*) '    Saving block no. ',ISBLK
                write(6,*) '     Block offset:   ',IOFF
                write(6,*) '     Block length:   ',LEN_V
              end if
              call itods(LEN_V,1,-1,LU)
              call todsc(VECT(IOFF),LEN_V,-1,LU)
            end do
*
          end do
*         ^ End loop over S batches of given MK2 value
        END DO
*       ^ End of loop over IMK2
*
* Write end of vector mark between real and imaginary part.
        call itods(-1,1,-1,LU)
*
        IRIOFF = LENGTH
      END DO
*     ^ End of loop over IRI
* Write end of vector mark.
      call itods(-1,1,-1,LU)
*
      RETURN
      END
