#ifdef MOD_ECP
C    -- FILE: argoscom.h --
c    This is for the argos routine
      INTEGER ngenag, nsag, naordsag, nconsag, ngcsag, itolag,          &
     &        icutag, aointsag, only1eag, inruag, ncrsag, l1recag,      &
     &        l2recag, aoint2ag, fsplitag, nstag, ndag(0:7),            &
     &        ndptag, maxang(60), irowmax(60,60), ireparg(60,60,100)
!CSK   fixed alignment: REAL, INTERGER
      COMMON /argoscom/ ngenag, nsag, naordsag, nconsag, ngcsag, itolag,&
     &                  icutag, aointsag, only1eag, inruag, ncrsag,     &
     &                  l1recag, l2recag, aoint2ag, fsplitag, nstag,    &
     &                  ndag, ndptag, maxang, irowmax, ireparg
      INTEGER ndptags(0:7)
      data ndptags/ 0,0,0,0,1,1,1,7/
C    -- end of argoscom.h --
#endif
