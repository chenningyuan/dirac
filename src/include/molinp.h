!
! File: molinp.h
!
! Contains a copy of the .mol input file in MLINE(1:NMLINE).
! This is done by RDLINE in "abacus/herrdn.F".
! NCLINE(i) points to the line defining atom center no. i
! NMLAU     points to the line in which it is defined if coordinates
!           are in AU (bohr) or in Angstroem ("line4")
!
      INTEGER    KMLINE, len_MLINE
      PARAMETER (KMLINE = 1000, len_MLINE = 80)

      CHARACTER*(len_MLINE) MLINE(KMLINE)
      COMMON /MOLINP_C/ MLINE

      INTEGER         NMLINE, NCLINE(MXCENT), NMLAU
      COMMON /MOLINP/ NMLINE, NCLINE, NMLAU
! -- End of molinp.h --
