C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : dcbkrmc.h
C SOURCE  : $Source$
C
C*** KR-MCSCF calculations ***
      PARAMETER (NMCFLAG = 20)
      LOGICAL MCFLAG(NMCFLAG),KRMC_SKIPEE,KRMC_SKIPEP,KRMC_NOPFQ
      CHARACTER KRMC_CIPROGRAM*8
      COMMON /DCLKRMC/ MCFLAG,KRMC_SKIPEE,KRMC_SKIPEP,KRMC_NOPFQ
      COMMON/DCOKRMC/
     +              NKRMCISH(2),NKRMCASH(2),NKRMCOCC(2),NKRMCSSH(2),
     +              NKRMCAS1(2),NKRMCAS2(2),NKRMCAS3(2),
     +              NKRMCESH(2),NKRMCPSH(2),NKRMCFRO(2),
     +              NKRMCELEC(2),NKRMCAELEC 
      COMMON /DCIKRMC/ IKRMC_ITRINT(2),
     &                 IKRMC_INTDEF,
     &                 IKRMC_INTFLG,
     &                 IKRMC_INTBUF,
     &                 MKRMC_MXMICRO,
     &                 MKRMC_MXMACRO,
     &                 IKRMC_ISTATE,
     &                 IKRMC_SYMMETRY,
     $     IPRKRMC
      COMMON /DCRKRMC/ DKRMC_CNVINT(2),DKRMC_THRGRD
      COMMON /DCCKRMC/ KRMC_CIPROGRAM
