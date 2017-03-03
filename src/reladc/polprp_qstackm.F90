module qstack

       interface

         FUNCTION qstack_init (i1,i2) bind(c, name="qstack_init")
           use, intrinsic ::  iso_c_binding
           integer (c_int)         :: qstack_init
           integer (c_int)         :: i1,i2
         END FUNCTION
       
         FUNCTION qstack_push (i1, i2, ra1) bind(c, name="qstack_push")
           use, intrinsic ::  iso_c_binding
           integer (c_int)         :: qstack_push
           integer (c_int)         :: i1, i2
           real (c_double)         :: ra1(*)
         END FUNCTION
       
         FUNCTION qstack_popf (i1, ra1) bind(c, name="qstack_popf")
           use, intrinsic ::  iso_c_binding
           integer (c_int)         :: qstack_popf
           integer (c_int)         :: i1
           real (c_double)         :: ra1(*)
         END FUNCTION
       
         FUNCTION qstack_popl (i1, ra1) bind(c, name="qstack_popl")
           use, intrinsic ::  iso_c_binding
           integer (c_int)         :: qstack_popl
           integer (c_int)         :: i1
           real (c_double)         :: ra1(*)
         END FUNCTION
       
         FUNCTION qstack_meminf () bind(c, name="qstack_meminf")
           use, intrinsic ::  iso_c_binding
           real (c_double)         :: qstack_meminf
         END FUNCTION
       
         SUBROUTINE qstack_shutdown() bind(c, name="qstack_shutdown")
           use, intrinsic :: iso_c_binding
         END SUBROUTINE

       end interface

end module
