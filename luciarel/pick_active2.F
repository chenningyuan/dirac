* General routine for picking blocks of elements over
* (e.g. active) indices from a full lower triangle matrix over
* all AO symmetries. (Optional: square matrix)
*
*  Array structure is given by NFCT(NINNER,NOUTER) for rows
*                       and by NFCT2(NINNER2,NOUTER2) for columns
*  IINNER gives the index of the desired block(s) for row index
*  IINNER2 for column index
*   in the array RMATIN.
*  It is assumed that the outer index runs over all present blocks.
*
      subroutine pick_active2(NOUTER,NINNER,NOUTER2,NINNER2,
     &                        NFCT,NFCT2,
     &                        IINNER,IINNER2,
     &                        TRIAIN,TRIAOUT,
     &                        RMATIN,RMATOUT,NTRIAIN,NTRIAOUT,
     &                        IFACIN,IFACOUT,IDIMIN,
     &                        IDIMOUTR,IDIMOUTC,
     &                        IINFORM,IOUTFORM,IPRINT)
*
*  Modified from pick_active for substructure arrays
*  Timo Fleig, June 2000 (it's a hot summer!)
*
      implicit real*8 (A-H,O-Z)
*
      dimension NFCT(NINNER,NOUTER),NFCT2(NINNER2,NOUTER2)
      dimension TRIAIN(NTRIAIN),TRIAOUT(NTRIAOUT),
     &          RMATIN(IDIMIN,IDIMIN),
     &          RMATOUT(IDIMOUTR,IDIMOUTC)
*
*  IINFORM = 1: Input lower triangle
*            2: Input full square
*  IOUTFORM = 1: Output lower triangle
*             2: Output full square
*
      NTESTL = 05
      NTEST = max(NTESTL,IPRINT)
*
      if (NTEST.ge.100) then
        write(6,*) 'Checking input matrices in pick_active'
        if (IINFORM.eq.1) then
          write(6,*) 'Input array triangle'
          call wrtmat(TRIAIN,1,NTRIAIN,1,NTRIAIN)
        else if (IINFORM.eq.2) then
          write(6,*) 'Input array block'
          call wrtmat(RMATIN,IDIMIN,IDIMIN,IDIMIN,IDIMIN)
        else
          write(6,*) 'Problem in pick_active.'
          write(6,*) 'IINFORM has disallowed value ',IINFORM
          Call Abend2( 'Quitting.' )
        end if
      end if
*
      if (IINFORM.eq.1) then
        call unfold(RMATIN,TRIAIN,IDIMIN,1,NTRIAIN,IFACIN)
      end if
*
*  J refers to columns, I to rows
      do JRR=1,NOUTER2,1
        do IRR=1,NOUTER,1
*
*  Input start row
*  Output start row
          IINSTA = 1
          IOUTSTA = 1
          do KRR=1,IRR-1,1
            do INN=1,NINNER,1
              IINSTA = IINSTA + NFCT(INN,KRR)
            end do
            IOUTSTA = IOUTSTA + NFCT(IINNER,KRR)
          end do
          do INN=1,IINNER-1,1
            IINSTA = IINSTA + NFCT(IINNER,IRR)
          end do
          IINEND = IINSTA + NFCT(IINNER,IRR) - 1
*
*  Input start column
*  Output start column
          JINSTA = 1
          JOUTSTA = 1
          do LRR=1,JRR-1,1
            do JNN=1,NINNER2,1
              JINSTA = JINSTA + NFCT2(JNN,LRR)
            end do
            JOUTSTA = JOUTSTA + NFCT2(IINNER2,LRR)
          end do
          do JNN=1,IINNER2-1,1
            JINSTA = JINSTA + NFCT2(IINNER2,JRR)
          end do
          JINEND = JINSTA + NFCT2(IINNER2,JRR) - 1
*
*  Do the job ...
          if (NTEST.ge.10) then
            write(6,*)
            write(6,*) 'IRR  ......   ',IRR
            write(6,*) 'JRR  ......   ',JRR
            write(6,*) 'IINSTA ....    ',IINSTA
            write(6,*) 'IINEND ....    ',IINEND
            write(6,*) 'JINSTA ....    ',JINSTA
            write(6,*) 'JINEND ....    ',JINEND
            write(6,*) 'IOUTSTA ...    ',IOUTSTA
            write(6,*) 'JOUTSTA ...    ',JOUTSTA
          end if
          ICOUNT = 0
          JSTART = JOUTSTA
          ISTART = IOUTSTA
          do JCOL=JINSTA,JINEND,1
            do I=IINSTA,IINEND,1
              RMATOUT(ISTART,JSTART) = RMATIN(I,JCOL)
              ISTART = ISTART + 1
              ICOUNT = ICOUNT + 1
            end do
            ISTART = IOUTSTA
            JSTART = JSTART + 1
          end do
          if (NTEST.ge.10) then
            write(6,*) 'No. of integrals transferred: '
            write(6,*) ICOUNT
          end if
        end do
      end do
*
      if (IOUTFORM.eq.1) then
        call fold(TRIAOUT,RMATOUT,IDIMOUT,1,NTRIAOUT)
        if (NTEST.ge.3) then
          write(6,*) 'Checking output from pick_active '
          do INT=1,NTRIAOUT,1
            write(6,'(I6,1F18.14)') INT,TRIAOUT(INT)
          end do
        end if
      else if (IOUTFORM.eq.2) then
        if (NTEST.ge.3) then
          write(6,*) 'Checking output from pick_active '
          call wrtmat(RMATOUT,IDIMOUTR,IDIMOUTC,IDIMOUTR,IDIMOUTC)
        end if
      end if
*
      return
      end
