#!/bin/bash

for my_file in "fsccsd_restart_s00" "fsccsd_restart_s01" "fsccsd_restart_s02" "fsccsd_restart_s10" "fsccsd_restart_s11" "fsccsd_restart_s20" ; do
   for my_ext in "_a_Mg.out" "_b_Mg.out"; do

      diff {,result/}$my_file""$my_ext > diff-files/$my_file"."$my_ext".diff"

   done
done

diff mpi-{1,4}p/fsccsd_original_Mg.out  > diff-files/fsccsd_original_Mg.diff
diff mpi-{1,4}p/fsccsd_Mg.out           > diff-files/fsccsd_Mg.diff
