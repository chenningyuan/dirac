*  Reorder Two-electron integrals from MOLECULE symmetry blocks
*  ordering to LUCIAREL symmetry block ordering
*  (A_u and B_g have to be interchanged)
*
      subroutine reotwin_c2h(SWEINT,HELPVEC,MOLOFF,LUCOFF,
     &                       MOLNIT,LUCNIT,IGLSUM_M,IGLSUM_L,
     &                       NSMBLK,NSMBLK_L,ILOCLEN)
      implicit real*8 (A-H,O-Z)
*
      dimension SWEINT(IGLSUM_M),HELPVEC(IGLSUM_M)
      dimension MOLOFF(NSMBLK),LUCOFF(NSMBLK_L)
      dimension MOLNIT(NSMBLK),LUCNIT(NSMBLK_L)
*
      NTESTL = 500
      NTEST = max(NTESTL,IPRINT)
*
      if (IGLSUM_M.ne.IGLSUM_L) then
         write(6,*) 'Problem in rdmole2'
         write(6,*) 'IGLSUM_M =  ',IGLSUM_M
         write(6,*) 'IGLSUM_L =  ',IGLSUM_L
         Call Abend2( 'Quitting.' )
      end if
*
      if (2*NSMBLK.gt.ILOCLEN) then
        write(6,*) 'Problem in rdmole2'
        write(6,*) '2*NSMBLK =  ',2*NSMBLK,'larger than '
        write(6,*) 'ILOCLEN =  ',ILOCLEN
        Call Abend2( 'Quitting.' )
      end if
*
      do ISMBLK=1,NSMBLK,1
         if (NTEST.ge.100) write(6,*) 'Doing ISMBLK =  ',ISMBLK
         INOFF = MOLOFF(ISMBLK)
         if (ISMBLK.eq.1) IOUTOFF = LUCOFF(1)
         if (ISMBLK.eq.2) IOUTOFF = LUCOFF(5)
         if (ISMBLK.eq.3) IOUTOFF = LUCOFF(7)
         if (ISMBLK.eq.4) IOUTOFF = LUCOFF(9)
         if (ISMBLK.eq.5) IOUTOFF = LUCOFF(11)
         if (ISMBLK.eq.6) IOUTOFF = LUCOFF(15)
         if (ISMBLK.eq.7) IOUTOFF = LUCOFF(16)
         if (ISMBLK.eq.8) IOUTOFF = LUCOFF(18)
         if (ISMBLK.eq.9) IOUTOFF = LUCOFF(19)
         if (ISMBLK.eq.10) IOUTOFF = LUCOFF(14)
         if (ISMBLK.eq.11) IOUTOFF = LUCOFF(2)
         if (ISMBLK.eq.12) IOUTOFF = LUCOFF(10)
         if (ISMBLK.eq.13) IOUTOFF = LUCOFF(6)
         if (ISMBLK.eq.14) IOUTOFF = LUCOFF(12)
         if (ISMBLK.eq.15) IOUTOFF = LUCOFF(13)
         if (ISMBLK.eq.16) IOUTOFF = LUCOFF(3)
         if (ISMBLK.eq.17) IOUTOFF = LUCOFF(8)
         if (ISMBLK.eq.18) IOUTOFF = LUCOFF(17)
         if (ISMBLK.eq.19) IOUTOFF = LUCOFF(4)
         if (NTEST.ge.100) write(6,*) 'INOFF =  ',INOFF
         if (NTEST.ge.100) write(6,*) 'IOUTOFF =  ',IOUTOFF
         do I=1,MOLNIT(ISMBLK),1
            HELPVEC(IOUTOFF+I-1) = SWEINT(INOFF+I-1)
         end do
         if (NTEST.ge.100) then
            write(6,*) 'transferred integrals:'
            do IWR=1,MOLNIT(ISMBLK),1
               write(6,*) SWEINT(INOFF+IWR-1)
            end do
         end if
      end do
*
      do ICOP=1,IGLSUM_M,1
         SWEINT(ICOP) = HELPVEC(ICOP)
      end do
*
      if (NTEST.ge.100) then
         write(6,*) 'Reordered integrals on SWEINT:'
         do I=1,IGLSUM_M,1
            write(6,*) I,SWEINT(I)
         end do
      end if
*
      return
      end
