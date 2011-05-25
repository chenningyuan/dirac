c
c     *** atomic parameter ***
      integer    mxatom,    mxiqm,   mxjco,    mxnucij,    mxnrcij,
     &           mxnbf,     mxngcs,    mxicsu,     mxictu,     mxnst
      parameter (mxatom=10, mxiqm=6, mxjco=10, mxnucij=20, mxnrcij=20,
     &           mxnbf=20,  mxngcs=50, mxicsu=100, mxictu=100, mxnst=20)
c
c     *** numeric parameter ***
      real       D0,       D1,       tolsq
      parameter (D0=0.0D0, D1=1.0D0, tolsq=1.0d-12) 
c
c     *** read lines ***
      integer   mxmline
      parameter (mxmline=1000)
      integer   nmline
      character mline(mxmline)*90
c
c     *** files ***
      integer lufile2,  ! molecule.inp
     &        lufile3,  ! argosin
     &        lufile4,  ! cnvrtin
     &        lufile6   ! library(ecp & bas)
c
c     *** general variables ***
      integer   i,j,k,l,m
