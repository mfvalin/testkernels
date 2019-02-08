  interface
    function vsearch_setup(levels, nk, ni, nj) result (ztab) bind(C,name='Vsearch_setup')
      import :: C_PTR, C_DOUBLE, C_INT
      real(C_DOUBLE), dimension(nk), intent(IN) :: levels
      integer(C_INT), intent(IN), value :: nk, ni, nj
      type(C_PTR) :: ztab
    end function vsearch_setup
    function vsearch_setup_plus(levels, nk, ni, nj, offseti, offsetj) result (ztab) bind(C,name='Vsearch_setup_plus')
      import :: C_PTR, C_DOUBLE, C_INT
      real(C_DOUBLE), dimension(nk), intent(IN) :: levels
      integer(C_INT), intent(IN), value :: nk, ni, nj, offseti, offsetj
      type(C_PTR) :: ztab
    end function vsearch_setup_plus
    subroutine tricublin_zyx1_n(d,f,pxyz,lv,n) bind(C,name='Tricublin_zyx1_n')
      import :: C_PTR
      real, dimension(*), intent(OUT)   :: d
      real, dimension(*), intent(IN)    :: f
      real, dimension(*), intent(IN)  :: pxyz
      type(C_PTR), intent(IN), value    :: lv
      integer, intent(IN), value        :: n
    end subroutine tricublin_zyx1_n
    subroutine tricublin_mono_zyx_n(d,l,mi,ma,f,pxyz,lv,n) bind(C,name='Tricublin_mono_zyx_n')
      import :: C_PTR
      real, dimension(*), intent(OUT)   :: d, l, mi, ma
      real, dimension(*), intent(IN)    :: f
      real, dimension(*), intent(IN)  :: pxyz
      type(C_PTR), intent(IN), value    :: lv
      integer, intent(IN), value        :: n
    end subroutine tricublin_mono_zyx_n
    subroutine tricublin_zyx3_n(d,f,pxyz,lv,n) bind(C,name='Tricublin_zyx3_n')
      import :: C_PTR
      real, dimension(*), intent(OUT)   :: d
      real, dimension(*), intent(IN)    :: f
      real, dimension(*), intent(IN)  :: pxyz
      type(C_PTR), intent(IN), value    :: lv
      integer, intent(IN), value        :: n
    end subroutine tricublin_zyx3_n
  end interface
