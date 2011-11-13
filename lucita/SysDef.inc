*----------------------------------------------------------------------*
* Define data type conversion constants                                *
* ItoB : integer --> byte                                              *
* RtoB : real*8  --> byte                                              *
* RtoI : real*8  --> integer                                           *
* (units=Bytes)                                                        *
*----------------------------------------------------------------------*
      Integer    ItoB,      RtoB,      RtoI
      Character Arch_Type*13

      Parameter( ItoB = 8 , RtoB = 8 , RtoI = 2  )

      Data Arch_Type/'IRIX64'/
