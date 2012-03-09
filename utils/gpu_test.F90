!----------------------------------------------
!  simple test of Fortran GPU CUDA libraries
!----------------------------------------------

Program GPU_math_tests

external cublas_init, cublas_set_matrix, cublas_get_matrix,cublas_set_vector,cublas_get_vector
external cublas_shutdown, cublas_alloc
integer cublas_alloc, cublas_set_matrix, cublas_get_matrix, cublas_get_vector


  print *,'hello'
  call cublas_init
  call cublas_shutdown


End Program GPU_math_tests



