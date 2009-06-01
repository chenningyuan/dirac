
#if defined STANDALONE_HSFS
 Program HIGH_SECTOR_FSCC_DRIVER
#else
 Subroutine HIGH_SECTOR_FSCC_DRIVER
#endif

  use high_sect_execs

#if defined STANDALONE_HSFS
  call get_line_arguments
#endif

  call launch_hsfs_methods

 End


