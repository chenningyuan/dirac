      subroutine geth0(H,LUONE,NSYM,NORB,NFRO_L,
     &                 FNONEMO,IPRINT)
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Load one electron integrals                                      *
*                                                                      *
*     Calling parameters:                                              *
*     H   : core Hamiltonian matrix                                    *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*                                                                      *
*   modified for LUCIAREL                                              *
*                                                                      *
*   Timo Fleig, September 1999                                         *
*                                                                      *
************************************************************************
*
      implicit Real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
      include 'SysDef.inc'
*
      dimension H(*)
      dimension NORB(NSYM),IOLIST(64),ITOC(64),NBAS_L(8),
     &          NORB_L(8),NFRO_L(8),NDEL_L(8)
*
      common/CECORE/ECORE
*
      character*6 FNONEMO
      character*4 NAME(2,MXPORB)
*
      NTESTL = 0
      NTEST = max(NTESTL,IPRINT)
*
*----------------------------------------------------------------------*
*     Start procedure:                                                 *
*     open the transformed one-electron integral file                  *
*----------------------------------------------------------------------*
      call daname(LUONE,FNONEMO)
*----------------------------------------------------------------------*
*     Set up the scatter/gather list and read table of contents        *
*----------------------------------------------------------------------*
      IDISK=0
      Call GSList(IOLIST,8,ITOC,64,POTNUC,RtoI,NSYM_L,1,
     &            NBAS_L,8,NORB_L,8,NFRO_L,8,NDEL_L,8,
     &            NAME,2*4*MXPORB/ItoB)
      call dafile(LUONE,4,IOLIST,IDUM,IDISK)
*
*----------------------------------------------------------------------*
*     Determine the number of integrals (symmetry blocked)             *
*----------------------------------------------------------------------*
      NorbTT=0
      Do iSym=1,nSym
        NorbTT=NorbTT+(Norb(iSym)*Norb(iSym)+Norb(iSym))/2
      End Do
*----------------------------------------------------------------------*
*     Load the core Hamiltonian matrix                                 *
*----------------------------------------------------------------------*
      IDISK = ITOC(2)
      Call DaFile(LuOne,2,H,RtoI*NorbTT,IDISK)
*
      ECORE = POTNUC
*
*  Check what has been loaded for consistency
      ISTOP = 0
      do ISM=1,NSYM,1
        if (NORB_L(ISM).ne.NORB(ISM)) then
          write(6,*) 'Problem in geth0.'
          write(6,*) 'NORB_L from TRAINT does not match NORB'
          write(6,*) 'NORB_L(',ISM,') = ',NORB_L(ISM)
          write(6,*) 'NORB(',ISM,') = ',NORB(ISM)
          ISTOP = 1
        end if
      end do
      if (ISTOP.eq.1) call Abend2( 'Quitting.' )
*
      if (NTEST.ge.1) then
        write(6,*)
        write(6,*) '*********************************************'
        write(6,*) 'Doing CI with FOCK matrix from MOTRA         '
        write(6,*) '*********************************************'
        write(6,*)
        write(6,*) 'Frozen orbitals: '
        do ISM=1,NSYM,1
          write(6,'(A7,I3,A4,I3)') 'NFRO_L(',ISM,') = ',NFRO_L(ISM)
        end do
      end if
*
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) 'Checking read integrals from file ',FNONEMO
        do INT=1,NORBTT,1
          write(6,'(I7,1F18.14)') INT,H(INT)
        end do
        write(6,*)
        write(6,*) 'Core energy incl. Fock part (h_ii,g_4in):  ',ECORE
      end if
*
      return
      end
