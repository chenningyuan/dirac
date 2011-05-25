! argos_files.h

      integer   nunits
      parameter(nunits=4)
      integer        iunits
      common /units/ iunits(nunits)
      integer                  nlist
      equivalence ( iunits(1), nlist )
      integer                  aoints
      equivalence ( iunits(2), aoints )
!     integer                  ninput
!     equivalence ( iunits(4), ninput )
 

