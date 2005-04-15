C
C     Fine structure constant
C

CMI   PARAMETER (CVEL = 137.0359895D0, ALPHAC = 1.0D0/CVEL,
CMI  &                               ALPHA2 = ALPHAC*ALPHAC)

CMI ... use the CVEL from the DIRAC ! If user changes CVEL there 
CMI    (ie to the older value),
CMI    it has to be changed in AMFI as well !!!!

       REAL*8 CVEL,ALPHAC,ALPHA2
CMI parameter (speed2=CVEL*CVEL), parameter (speed4=speed2*speed2) 
       REAL*8 speed2,speed4
       COMMON /AMFICB/ CVEL,ALPHAC,ALPHA2,speed2,speed4
