subroutine test(AB,A,b,ni,nj,nk)
  integer :: ni,nj,nk
  real, dimension(ni,nj,nk) :: a, b
  real, dimension(2,ni,nj,nk) :: AB
  integer :: i,j,k
#if defined(BAD)
  ab(1,:,:,:) = a(:,:,:)
  ab(2,:,:,:) = b(:,:,:)
#else
  do k=1,nk
  do j=1,nj
  do i=1,ni
    ab(1,i,j,k) = a(i,j,k)
    ab(2,i,j,k) = b(i,j,k)
  enddo
  enddo
  enddo
#endif
  return
end
