module macro9
!
!***********************************************************************
!  ----- functions computing various forms of products of     -----
!  ----- complex numbers given their real and imaginary parts -----
!  ----- these functions should be automatically inlined      -----
!=======================================================================
!     ZPXY(a,b,x,y) computes  A*X  (A = a+Ib and X = x+Iy)
!     Z = R for real part
!         I pour imaginary part
!     XY defines how to use A=a+Ib and X=x+Iy
!   X,Y = N the number is in 'normal' form (x+Iy)
!         C the number is in 'conjugate' form (x-Iy)
!         So:
!         NN => (a+Ib)*(x+Iy)
!         CN => (a-Ib)*(x+Iy)
!         NC => (a+Ib)*(x-Iy)
!         CC => (a-Ib)*(x-Iy)
!***********************************************************************
!
   implicit none
   contains
   function RPNN(ajp,bjp,xjp,yjp) result(value)
      double precision, intent(in) :: ajp,bjp,xjp,yjp
      double precision :: value
      value = ajp*xjp-bjp*yjp
      return
   end function RPNN
   function IPNN(ajp,bjp,xjp,yjp) result(value)
      double precision, intent(in) :: ajp,bjp,xjp,yjp
      double precision :: value
      value = ajp*yjp+bjp*xjp
      return
   end function IPNN
   function RPCN(ajp,bjp,xjp,yjp) result(value)
      double precision, intent(in) :: ajp,bjp,xjp,yjp
      double precision :: value
      value = ajp*xjp+bjp*yjp
      return
   end function RPCN
   function IPCN(ajp,bjp,xjp,yjp) result(value)
      double precision, intent(in) :: ajp,bjp,xjp,yjp
      double precision :: value
      value = ajp*yjp-bjp*xjp
      return
   end function IPCN
   function RPNC(ajp,bjp,xjp,yjp) result(value)
      double precision, intent(in) :: ajp,bjp,xjp,yjp
      double precision :: value
      value = ajp*xjp+bjp*yjp
      return
   end function RPNC
   function IPNC(ajp,bjp,xjp,yjp) result(value)
      double precision, intent(in) :: ajp,bjp,xjp,yjp
      double precision :: value
      value = bjp*xjp-ajp*yjp
      return
   end function IPNC
   function RPCC(ajp,bjp,xjp,yjp) result(value)
      double precision, intent(in) :: ajp,bjp,xjp,yjp
      double precision :: value
      value = ajp*xjp-bjp*yjp
      return
   end function RPCC
   function IPCC(ajp,bjp,xjp,yjp) result(value)
      double precision, intent(in) :: ajp,bjp,xjp,yjp
      double precision :: value
      value = -(ajp*yjp+bjp*xjp)
      return
   end function IPCC
end module macro9
