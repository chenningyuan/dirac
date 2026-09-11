!-----------------------------------------
! Module with variables previously defined
! at nucata.h include file
!-----------------------------------------
   module nucdata_mod
!
   implicit none
   save
!
   ! Element symbols (1–138)
   character(len=4), parameter :: NUCLABEL(1:138) = (/ &
        'H  ', 'He ', 'Li ', 'Be ', 'B  ',                      &
        'C  ', 'N  ', 'O  ', 'F  ', 'Ne ',                      &
        'Na ', 'Mg ', 'Al ', 'Si ', 'P  ',                      &
        'S  ', 'Cl ', 'Ar ', 'K  ', 'Ca ',                      &
        'Sc ', 'Ti ', 'V  ', 'Cr ', 'Mn ',                      &
        'Fe ', 'Co ', 'Ni ', 'Cu ', 'Zn ',                      &
        'Ga ', 'Ge ', 'As ', 'Se ', 'Br ',                      &
        'Kr ', 'Rb ', 'Sr ', 'Y  ', 'Zr ',                      &
        'Nb ', 'Mo ', 'Tc ', 'Ru ', 'Rh ',                      &
        'Pd ', 'Ag ', 'Cd ', 'In ', 'Sn ',                      &
        'Sb ', 'Te ', 'I  ', 'Xe ', 'Cs ',                      &
        'Ba ', 'La ', 'Ce ', 'Pr ', 'Nd ',                      &
        'Pm ', 'Sm ', 'Eu ', 'Gd ', 'Tb ',                      &
        'Dy ', 'Ho ', 'Er ', 'Tm ', 'Yb ',                      &
        'Lu ', 'Hf ', 'Ta ', 'W  ', 'Re ',                      &
        'Os ', 'Ir ', 'Pt ', 'Au ', 'Hg ',                      &
        'Tl ', 'Pb ', 'Bi ', 'Po ', 'At ',                      &
        'Rn ', 'Fr ', 'Ra ', 'Ac ', 'Th ',                      &
        'Pa ', 'U  ', 'Np ', 'Pu ', 'Am ',                      &
        'Cm ', 'Bk ', 'Cf ', 'Es ', 'Fm ',                      &
        'Md ', 'No ', 'Lr ', 'Rf ', 'Db ',                      &
        'Sg ', 'Bh ', 'Hs ', 'Mt ', 'Ds ',                      &
        'Rg ', 'Cn ', 'Nh ', 'Fl ', 'Mc ',                      &
        'Lv ', 'Ts ', 'Og ', 'Uee', 'Ubn',                      &
        'Ubu', 'Ubb', 'Ubt', 'Ubq', 'Ubp',                      &
        'Ubh', 'Ubs', 'Ubo', 'Ube', 'Utn',                      &
        'Utu', 'Utb', 'Utt', 'Utq', 'Utp',                      &
        'Uth', 'Uts', 'Uto' /)
!
   ! Nuclear mass refers to most abundant or most stable isotope
   real(kind=8), parameter :: NUCMASS(1:138) = (/ &
         1.d0,   4.d0,   7.d0,   9.d0,  11.d0,      &
        12.d0,  14.d0,  16.d0,  19.d0,  20.d0,      &
        23.d0,  24.d0,  27.d0,  28.d0,  31.d0,      &
        32.d0,  35.d0,  40.d0,  39.d0,  40.d0,      &
        45.d0,  48.d0,  51.d0,  52.d0,  55.d0,      &
        56.d0,  59.d0,  58.d0,  63.d0,  64.d0,      &
        69.d0,  74.d0,  75.d0,  80.d0,  79.d0,      &
        84.d0,  85.d0,  88.d0,  89.d0,  90.d0,      &
        93.d0,  98.d0,  98.d0, 102.d0, 103.d0,      &
       106.d0, 107.d0, 114.d0, 115.d0, 120.d0,      &
       121.d0, 130.d0, 127.d0, 132.d0, 133.d0,      &
       138.d0, 139.d0, 140.d0, 141.d0, 144.d0,      &
       145.d0, 152.d0, 153.d0, 158.d0, 159.d0,      &
       162.d0, 162.d0, 168.d0, 169.d0, 174.d0,      &
       175.d0, 180.d0, 181.d0, 184.d0, 187.d0,      &
       192.d0, 193.d0, 195.d0, 197.d0, 202.d0,      &
       205.d0, 208.d0, 209.d0, 209.d0, 210.d0,      &
       222.d0, 223.d0, 226.d0, 227.d0, 232.d0,      &
       231.d0, 238.d0, 237.d0, 244.d0, 243.d0,      &
       247.d0, 247.d0, 251.d0, 252.d0, 257.d0,      &
       258.d0, 259.d0, 262.d0, 261.d0, 262.d0,      &
       263.d0, 262.d0, 265.d0, 266.d0, 281.d0,      &
       280.d0, 285.d0, 284.d0, 289.d0, 288.d0,      &
       292.d0, 293.d0, 294.d0, 303.d0, 304.d0,      &
       330.d0, 334.d0, 337.d0, 340.d0, 344.d0,      &
       347.d0, 350.d0, 354.d0, 357.d0, 361.d0,      &
       364.d0, 367.d0, 371.d0, 374.d0, 378.d0,      &
       381.d0, 385.d0, 388.d0 /)
!
   end module nucdata_mod
