program test_kernel_matvec2
  use kmv2
  implicit none
  include 'mpif.h'
  real *8 :: t0, t1
  call initialize()
! call kernel_matvec2 (nu   , dt   , i0, in  , j0, jn  , ns)
  t0 = MPI_wtime()
  call kernel_matvec2 (1.0_8, 2.0_8, 4 , l_ni-3, 4 , l_nj-3, 3 )
  t1 = MPI_wtime()
  print *, 'time=',nint(t1-t0)*1000,' msecs'
end program test_kernel_matvec2
