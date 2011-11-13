*comdeck luci.inc $Revision: 1.3 $ 
C     
C     COMMON INFORMATION FOR LUCITA CI PROGRAM
C     
*                   Timo Fleig
*
      character*72 TITLE(MXTIT)
      character*6 TYPE,WAFFCT,PRILUC
      character*3 CALSIZ
      common/LUCIIF/TITLE,TYPE,WAFFCT,PRILUC,CALSIZ
      integer        IRREP,NROOT,NACTEL,MULTIP,MXHOL1,MXELR3,INACOB,
     &               IDEMOL,IRSTRT
      common/iLUCIIF/IRREP,NROOT,NACTEL,MULTIP,MXHOL1,MXELR3,INACOB,
     &               IDEMOL,IRSTRT
