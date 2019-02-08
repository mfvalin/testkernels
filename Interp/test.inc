  interface
    function nanocycles() result (nano) bind(C,name='Nanocycles')
    import :: C_LONG_LONG
    integer(C_LONG_LONG) :: nano
    end function nanocycles
  end interface
  integer, parameter :: NI = 300
  integer, parameter :: NJ = 200
  integer, parameter :: NK = 110
  integer, parameter :: NR = 10
  integer :: my_proc = 0
  integer*8 :: t0, t1
  real*8 :: dx, dy, dz
  real*8, dimension(8) :: px, py, pz
  real*8, dimension(8) :: pxy
  real*8, dimension(2) :: xy
  real, dimension(-2:NI+3,-1:NJ+2,NK), target :: f1, f2, f3
  real, dimension(3,-2:NI+3,-1:NJ+2,NK), target :: f123
  real*8, dimension(NK) :: zlev
  integer :: i, j, k, ii, jj, kk, rr
  real*8 fx, fy, fz, fxyz, x, y, z, ovni, ovnj, ovnk, xx, yy, zz
  real*4, dimension(3) :: r1, r2, r3, e
  type(C_PTR) :: lv
  real, dimension(3,NI,NJ,NK) :: pxpypz
  real, dimension(NI,NJ,NK) :: expected, d, dmin, dmax, dlin
  real, dimension(3,NI,NJ,NK) :: d3
  real*8, dimension(NK) :: levels
  real :: delta, error, delta1, delta2, delta3
  real*8 :: avg, avg1, avg2, avg3
  integer :: exact, minmaxerr, ierr

  fx(x) = (x+1.0)*(x+1.1)*(x+1.2)*(x+1.3)
  fy(y) = (y+1.05)*(y+1.15)*(y+1.25)*(y+1.35)
  fz(z) = (z*1.25+1.02)
!   fxyz(x,y,z) = x + y + z
!   fxyz(x,y,z) = x * y * z
!   fxyz(x,y,z) = fx(x*ovni) * fy(y*ovnj) * fz(z*ovnk)
  fxyz(x,y,z) = fx(x*ovni) + fy(y*ovnj) * fz(z*ovnk)
#if defined(USE_MPI)
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,my_proc,ierr)
#endif
  px = 0
  py = 0
  pz = 0
  ovni = 1.0_8 / (NI - 1)
  ovnj = 1.0_8 / (NJ - 1)
  ovnk = 1.0_8 / (NK - 1)
  
!   dx = .375
!   dy = .25
!   dz = .125
!   call tricubic_coeffs_d(px,py,pz,dx,dy,dz)
!   print 103,'px   = ',px(1:4),'py   = ',py(1:4),'pz   = ',pz(1:4)
!   print 103,'pzl  = ',pz(5:8)
!   print 103,'sums = ', sum(px(1:4)),sum(py(1:4)),sum(pz(1:4)),sum(pz(5:8))

  pxpypz = 0.0

  do k = 1, NK
    levels(k) = k
    do j = 1, NJ
      do i = 1, NI
        dx = .375 + i * .001
        dy = .251 + j * .001
        dz = .125 + k * .003
        xx = (i)
        yy = (j)
	zz = (k)
	f1(i,j,k) = real(fxyz(xx,yy,zz))
	f2(i,j,k) = f1(i,j,k) + 10.0
	f3(i,j,k) = f1(i,j,k) + 20.0
	f123(1,i,j,k) = f1(i,j,k)
	f123(2,i,j,k) = f2(i,j,k)
	f123(3,i,j,k) = f3(i,j,k)
! 	xx = min( max(i,2), ni - 2) + dx + 3         ! offset x de 3
! 	yy = min( max(j,2), nj - 2) + dy + 2         ! offset y de 2
	xx = min( max(i,2), ni - 2) + dx         ! offset x de 3
	yy = min( max(j,2), nj - 2) + dy         ! offset y de 2
	if(k > 1) then
	  zz = zz + dz
	else
	  zz = zz - dz
	endif
	pxpypz(1,i,j,k) = xx + 3        ! offset x de 3
	pxpypz(2,i,j,k) = yy + 2        ! offset y de 2
	pxpypz(3,i,j,k) = zz
! 	expected(i,j,k) = real(fxyz(xx-3,yy-2,zz)) ! compensation d'offset
	expected(i,j,k) = real(fxyz(xx,yy,zz))
      enddo
    enddo
  enddo
  dmin = 0
  dmax = 0
  dlin = 0

  lv = vsearch_setup_plus(levels, NK, NI+6, NJ+4, -3, -2)    ! "NI" = NI + 6, "NJ" = nj + 4 offseti = -3, offsetj = -2
