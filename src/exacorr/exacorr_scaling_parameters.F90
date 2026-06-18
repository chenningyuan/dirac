module exacorr_scaling_parameters

   implicit none
 
   complex(8), public, parameter :: ZERO=(0.D0,0.D0)
 
   complex(8), public, parameter :: ONE=(1.0D0,0.D0)
   complex(8), public, parameter :: MINUS_ONE=(-1.D0,0.D0)

   complex(8), public, parameter :: IMAGINARY=(0.0D0,1.0D0)
 
   complex(8), public, parameter :: TWO=(2.0D0,0.D0)
   complex(8), public, parameter :: MINUS_TWO=(-2.D0,0.D0)
 
   complex(8), public, parameter :: THREE=(3.0D0,0.D0)
   complex(8), public, parameter :: MINUS_THREE=(-3.D0,0.D0)

   complex(8), public, parameter :: SIX=(6.0D0,0.D0)
   complex(8), public, parameter :: MINUS_SIX=(-6.D0,0.D0)
 
   complex(8), public, parameter :: ONE_HALF=(0.5D0,0.D0)
   complex(8), public, parameter :: MINUS_ONE_HALF=(-0.5D0,0.D0)
 
   complex(8), public, parameter :: ONE_QUARTER=(0.25D0,0.D0)
   complex(8), public, parameter :: MINUS_ONE_QUARTER=(-0.25D0,0.D0)
 
   complex(8), public, parameter :: ONE_EIGHT=(0.125D0,0.D0)
   complex(8), public, parameter :: MINUS_ONE_EIGHT=(-0.125D0,0.D0)

   complex(8), public, parameter :: ONE_SIXTH=ONE/(6.0D0)
   complex(8), public, parameter :: MINUS_ONE_SIXTH=MINUS_ONE/(6.0D0)

   complex(8), public, parameter :: ONE_TWELFTH=ONE/(12.0D0)
   complex(8), public, parameter :: MINUS_ONE_TWELFTH=MINUS_ONE/(12.0D0)

   complex(8), public, parameter :: ONE_16TH=ONE/(16.0D0)
   complex(8), public, parameter :: MINUS_ONE_16TH=MINUS_ONE/(16.0D0)

   complex(8), public, parameter :: ONE_36TH=ONE/(36.0D0)
   complex(8), public, parameter :: MINUS_ONE_36TH=MINUS_ONE/(36.0D0)

   real(8), public, parameter    :: ONE_R=1.0D0
   real(8), public, parameter    :: ONE_HALF_R=0.5D0
   real(8), public, parameter    :: ONE_QUARTER_R=0.25D0
   real(8), public, parameter    :: ONE_SIXTH_R=(1.D0/6.D0)
   real(8), public, parameter    :: ONE_36_R=(1.D0/36.D0)

end module exacorr_scaling_parameters
