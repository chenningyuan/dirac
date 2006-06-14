#!/bin/csh 

foreach MOLSTEM ( \
#                 He     \
#                 He.vdz \
#                 He.vtz \
#                 He.bas1d \
#                 He7.sto \
#                 He7.vdz \
#                 He7.vtz \
#                 Ne.sto \
#                 Ne.vdz \
#                 Ne.vtz \
                  Ne7.vtz \
#                 He7.vtz \
                )

   foreach INPSTEM ( \
#                    disk           \
#                    disk.twocomp   \
#                    large.disk.twocomp   \
                     large.disk.twocomp.contracted \
#                    large.disk.iotc.contracted \
#                    disk.twocomp.contracted \
#                    direct         \
#                    direct.twocomp \
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

         ../../pam -mw 250 $MOL $INP 
#        ../../pam $MOL $INP -get gmon.out
#        mv gmon.out gmon.$MOL.$INP.out
#        gprof ../../../dirac.x gmon.$MOL.$INP.out > $MOL.$INP.profile
#        python gprof2tree.py --once --children _dirctl_ < $MOL.$INP.profile > pamdrv.$MOL.$INP.tree 

      end
   end
end
    
