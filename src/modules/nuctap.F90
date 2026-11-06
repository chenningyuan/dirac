!-----------------------------------------
! Module with variables previously defined
! at nuctap.h include file (COMMON /NUCTAP/, /NUCFIL/)
!-----------------------------------------
   module nuctap_mod
   implicit none
   integer :: LUMC, LUGDR, LUGDI, LUGDT, LURDR, LURDI, LURDT, LUY, LUSFDA, LUNR1, LUNR3, &
              LUNR5, LUWLK, LUKAP, LUXYDA, LUFJ1, LUINTJ, LUTGD, LUTRD, LUZGD, LUZRD, &
              LUTRJ, LUTLM, LUDFCK, LRGDR, LRGDI, LRGDT, LRRDR, LRRDI, LRRDT, LRSFDA, &
              LRKAP, LRXYDA, LRTGD, LRTRD, LRTLM, LRDFCK, NBGDR, NBGDI, NBGDT, NBRDR, &
              NBRDI, NBRDT, NBSFDA, NBKAP, NBXYDA, NBTGD, NBTRD, NBTLM, NBDFCK, LUFDC, &
              LFDCRC
   character(len=10) :: ABAMC, ABAGDR, ABAGDI, ABAGDT, ABARDR, ABARDI, ABARDT, ABAY, ABASF, &
                        ABAKAP, ABAXY, ABATGD, ABATRD, ABAZGD, ABAZRD, ABATRJ, ABAWLK, &
                        ABAIRC, ABATLM, ABADFK
   end module nuctap_mod
