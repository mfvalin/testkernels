program test_kernel_matvec2
  use kmv2
  implicit none
  call initialize()
! call kernel_matvec2 (nu   , dt   , i0, in  , j0, jn  , ns)
  call kernel_matvec2 (1.0_8, 2.0_8, 1 , l_ni, 1 , l_nj, 3 )
end program test_kernel_matvec2
