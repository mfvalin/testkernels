program test_kernel_matvec2
  use kmv2
  implicit none
  include 'mpif.h'
  integer :: ier, rank, irep
  real *8 :: t0, t1
  call mpi_init(ier)
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ier)
  call initialize()
  call mpi_barrier(MPI_COMM_WORLD,ier)
! call kernel_matvec2 (nu   , dt   , i0, in  , j0, jn  , ns)
  t0 = MPI_wtime()
  call mpi_barrier(MPI_COMM_WORLD,ier)
  do irep = 1, 1
    call kernel_matvec2 (1.0_8, 2.0_8, 1 , l_ni, 1 , l_nj, 3 )
  enddo
  call mpi_barrier(MPI_COMM_WORLD,ier)
  t1 = MPI_wtime()
  if(rank == 0) print *, 'time=',nint((t1-t0)*1000),' msecs'
  call mpi_finalize(ier)
end program test_kernel_matvec2
