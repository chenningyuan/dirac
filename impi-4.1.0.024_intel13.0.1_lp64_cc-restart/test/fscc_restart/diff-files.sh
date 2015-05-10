#!/bin/bash

for my_file in "fsccsd_restart_s00_Mg" "fsccsd_restart_s01_Mg" "fsccsd_restart_s02_Mg" "fsccsd_restart_s10_Mg" "fsccsd_restart_s11_Mg" "fsccsd_restart_s20_Mg" ; do
   for my_ext in "out" "out.0" "out.1" "out.2" "out.3" "out.4" "out.5" "out.6" "out.7"; do

      diff mpi-{1,4}p/$my_file"."$my_ext > diff-files/$my_file"."$my_ext".diff"

   done
done

diff mpi-{1,4}p/fsccsd_original_Mg.out  > diff-files/fsccsd_original_Mg.diff
diff mpi-{1,4}p/fsccsd_Mg.out           > diff-files/fsccsd_Mg.diff
