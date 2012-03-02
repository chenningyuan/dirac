#ifdef MOD_ECP
!    -- FILE: argoscom.h --
!    This is for the argos routine
      INTEGER ngenag, nsag, naordsag, nconsag, ngcsag, itolag,          &
     &        icutag, aointsag, only1eag, inruag, ncrsag, l1recag,      &
     &        l2recag, aoint2ag, fsplitag, nstag, ndag(0:7),            &
     &        ndptag, MAXANG(60), irowmax(60,60), ireparg(60,60,100),   &
     &        IGRPAG, ECPNUC, ECPNONT(600)
!CSK   fixed alignment: REAL, INTERGER
      COMMON /argoscom/ ngenag, nsag, naordsag, nconsag, ngcsag, itolag,&
     &                  icutag, aointsag, only1eag, inruag, ncrsag,     &
     &                  l1recag, l2recag, aoint2ag, fsplitag, nstag,    &
     &                  ndag, ndptag, MAXANG, irowmax, ireparg,         &
     &                  IGRPAG, ECPNUC, ECPNONT
      INTEGER ndptags(0:7)
      data ndptags/ 0,0,0,0,1,1,1,7/
!    -- end of argoscom.h --
#endif
