   module kmv2
     integer, parameter :: P=1, WWW=2, WW=3, W=4, E=5, EE=6, EEE=7, SSS=8, SS=9, S=10, N=11, NN=12, NNN=13, UPWIND_SHAPE=13
     integer, parameter :: www_stag=1, ww_stag=2, w_stag=3, e_stag=4, eee_stag=5, ee_stag=6, sss_stag=7, &
                           ss_stag=8, s_stag=9, n_stag=10, nn_stag=11, nnn_stag=12, CENTERED_SHAPE=12
     integer, parameter :: halox = 2
     integer, parameter :: haloy = 2
     integer, parameter :: l_ni = 200
     integer, parameter :: l_nj = 150
     integer, parameter :: l_nk = 100
     integer, parameter :: l_minx = 1 - halox
     integer, parameter :: l_maxx = l_ni - halox
     integer, parameter :: l_miny = 1 - haloy
     integer, parameter :: l_maxy = l_nj - haloy
     integer, parameter :: nvars = 9
     integer, parameter :: nphi = 2
     save
    real, dimension(l_ni) :: fx
    real, dimension(l_nj) :: fy
    real, dimension(l_minx:l_maxx, l_miny:l_maxy, l_nk, nvars) :: vec
    real, dimension(l_minx:l_maxx, l_miny:l_maxy, l_nk, nvars) :: vis
    real, dimension(l_minx:l_maxx, l_miny:l_maxy, l_nk, nvars) :: prod
    real, dimension(l_minx:l_maxx, l_miny:l_maxy, l_nk, nvars, nphi) :: u
    real, dimension(l_minx:l_maxx, l_miny:l_maxy, l_nk) :: v_vec_on_u, u_vec_on_v
    real*8, dimension(nphi) :: V_aug
    real*8, dimension(l_ni) :: sum1,sum2,sum3
    real*8, dimension(:,:,:,:), allocatable :: stencil_eq1_h, stencil_eq1_u, stencil_eq1_v
    real*8, dimension(:,:,:,:), allocatable :: stencil_eq2_h, stencil_eq2_u
    real*8, dimension(:,:,:,:), allocatable :: stencil_eq3_h, stencil_eq3_v
    real*8, dimension(:,:,:,:), allocatable :: stencil_adv_tr
    real*8, dimension(:,:,:), allocatable :: stencil_eq2_v_interp, stencil_eq3_u_interp
    contains
    subroutine initialize
      implicit none
      allocate( stencil_eq1_h(l_minx:l_maxx, UPWIND_SHAPE  , l_miny:l_maxy, l_nk), &
                stencil_eq1_u(l_minx:l_maxx, CENTERED_SHAPE, l_miny:l_maxy, l_nk),&
                stencil_eq1_v(l_minx:l_maxx, CENTERED_SHAPE, l_miny:l_maxy, l_nk) )

      allocate ( stencil_eq2_h(l_minx:l_maxx, CENTERED_SHAPE, l_miny:l_maxy, l_nk), &
                 stencil_eq2_u(l_minx:l_maxx, UPWIND_SHAPE  , l_miny:l_maxy, l_nk), &
                 stencil_eq2_v_interp(l_minx:l_maxx, l_miny:l_maxy, l_nk) )

      allocate ( stencil_eq3_h(l_minx:l_maxx, CENTERED_SHAPE, l_miny:l_maxy,l_nk), &
                 stencil_eq3_u_interp(l_minx:l_maxx, l_miny:l_maxy, l_nk), &
                 stencil_eq3_v(l_minx:l_maxx, UPWIND_SHAPE, l_miny:l_maxy, l_nk) )

      allocate ( stencil_adv_tr(l_minx:l_maxx, UPWIND_SHAPE, l_miny:l_maxy, l_nk) )

      stencil_eq1_h = 1.1
      stencil_eq1_u = 1.2
      stencil_eq1_v = 1.3

      stencil_eq2_h        = 1.4
      stencil_eq2_u        = 1.4
      stencil_eq2_v_interp = 1.6

      stencil_eq3_h        = 1.7
      stencil_eq3_u_interp = 1.8
      stencil_eq3_v        = 1.9

      stencil_adv_tr       = 1.0

      sum1 = 0
      sum2 = 0
      sum3 = 0
      V_aug = 1.35
      vec = 1.234
      vis = 2.345
      prod = 0
      u = 0
      fx = 1.34
      fy = 1.12
      v_vec_on_u = 1.45
      u_vec_on_v = 1.56
    end subroutine initialize
     
   end module kmv2
