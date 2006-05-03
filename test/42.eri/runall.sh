#!/bin/csh 

foreach MOLSTEM ( \
                  He     \
#                 He.vdz \
#                 He.vtz \
#                 He.bas1d \
                )

   foreach INPSTEM ( \
                     He.disk           \
                     He.direct         \
                     He.disk.twocomp   \
                     He.direct.twocomp \
                   )

      foreach INT ( \
                    eri \
                    her \
                  )


         if ( $INT == "eri" ) then
            set MOL=$MOLSTEM.sph
            set INP=$INPSTEM.eri
         else
            set MOL=$MOLSTEM.crt
            set INP=$INPSTEM.her
         endif 

#        ../../pam $MOL $INP 
         ../../pam $MOL $INP -get gmon.out
         mv gmon.out gmon.$MOL.$INP.out
         gprof ../../dirac.x gmon.$MOL.$INP.out > $MOL.$INP.profile
         python gprof2tree.py --once --children _dirctl_ < $MOL.$INP.profile > pamdrv.$MOL.$INP.tree 

   end
end
    
