! *** .Mol file *** ( Card 4 & 5 )
      integer   RECP_NONTYP, kcharg, symtxt, RECP_NONT(mxatom)
      real*8    thrs, q(mxatom), qexp1 
      character crt*1, kasym(3,3)*1, id3*1,
     &          namatm(mxatom)*3, namnar(mxatom,mxatom)*4
 
! *** argosin file ***
      character mtype(mxatom)*3 
      integer   nf(mxatom), nc(mxatom,mxatom), 
     &          igen, igenp(8), igenq(3), jgen, kgen, kgen1, kgen2,
     &          isend(8,mxatom,8,mxatom),mj2(mxatom,mxatom),
     &          ngenl, ngenp, mcons, igcs
      real*8    chg(mxatom), 
     &          xyz1(3,8,mxatom),        xyz2(3,8,mxatom,mxatom),
     &          xyz3(3,8,mxatom,mxatom), xyz4(3)


