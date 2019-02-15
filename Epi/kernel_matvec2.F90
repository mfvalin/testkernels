   subroutine kernel_matvec2 (nu, dt, i0, in, j0, jn, ns)

   use kmv2
   implicit none

   real*8, intent(in ) :: nu, dt

   integer, intent(IN) :: i0, in, j0, jn, ns
   integer :: i, j, k, v, pp, slc
   real fyy

   do k = 1,l_nk
!       print *,'k=',k
      do j = j0, jn
         do i = i0, in
            sum1(i) = sum1(i) + ( vec(i, j, k, 1)   * vec(i, j, k, 1) )
            prod(i,j,k,1) =  stencil_eq1_h(i,WWW,j,k)      * vec(i-3,j,k,1) &
                           + stencil_eq1_h(i,WW ,j,k)      * vec(i-2,j,k,1) &
                           + stencil_eq1_h(i,W  ,j,k)      * vec(i-1,j,k,1) &
                           + stencil_eq1_h(i,P  ,j,k)      * vec(i  ,j,k,1) &
                           + stencil_eq1_h(i,E  ,j,k)      * vec(i+1,j,k,1) &
                           + stencil_eq1_h(i,EE ,j,k)      * vec(i+2,j,k,1) &
                           + stencil_eq1_h(i,EEE,j,k)      * vec(i+3,j,k,1)
            prod(i,j,k,2) =  stencil_eq2_h(i,www_stag,j,k) * vec(i-2,j,k,1) &
                           + stencil_eq2_h(i,ww_stag ,j,k) * vec(i-1,j,k,1) &
                           + stencil_eq2_h(i,w_stag  ,j,k) * vec(i  ,j,k,1) &
                           + stencil_eq2_h(i,e_stag  ,j,k) * vec(i+1,j,k,1) &
                           + stencil_eq2_h(i,ee_stag ,j,k) * vec(i+2,j,k,1) & 
                           + stencil_eq2_h(i,eee_stag,j,k) * vec(i+3,j,k,1)
            prod(i,j,k,1) = prod(i,j,k,1) &
                           + stencil_eq1_h(i,SSS,j,k)      * vec(i,j-3,k,1) &
                           + stencil_eq1_h(i,SS ,j,k)      * vec(i,j-2,k,1) &
                           + stencil_eq1_h(i,S  ,j,k)      * vec(i,j-1,k,1) &
                           + stencil_eq1_h(i,N  ,j,k)      * vec(i,j+1,k,1) &
                           + stencil_eq1_h(i,NN ,j,k)      * vec(i,j+2,k,1) &
                           + stencil_eq1_h(i,NNN,j,k)      * vec(i,j+3,k,1)
            prod(i,j,k,3) =  stencil_eq3_h(i,sss_stag,j,k) * vec(i,j-2,k,1) &
                           + stencil_eq3_h(i,ss_stag ,j,k) * vec(i,j-1,k,1) &
                           + stencil_eq3_h(i,s_stag  ,j,k) * vec(i,j  ,k,1) &
                           + stencil_eq3_h(i,n_stag  ,j,k) * vec(i,j+1,k,1) &
                           + stencil_eq3_h(i,nn_stag ,j,k) * vec(i,j+2,k,1) &
                           + stencil_eq3_h(i,nnn_stag,j,k) * vec(i,j+3,k,1)
         end do
      end do

      do j = j0, jn
         do i = i0, in
            prod(i,j,k,1) = prod(i,j,k,1) &
               + stencil_eq1_u(i,www_stag,j,k) * vec(i-3,j,k,2) &
               + stencil_eq1_u(i,ww_stag ,j,k) * vec(i-2,j,k,2) &
               + stencil_eq1_u(i,w_stag  ,j,k) * vec(i-1,j,k,2) &
               + stencil_eq1_u(i,e_stag  ,j,k) * vec(i  ,j,k,2) &
               + stencil_eq1_u(i,ee_stag ,j,k) * vec(i+1,j,k,2) &
               + stencil_eq1_u(i,eee_stag,j,k) * vec(i+2,j,k,2)
            prod(i,j,k,2) = dt * ( prod(i,j,k,2) &
               + stencil_eq2_u(i,P  ,j,k)      * vec(i  ,j,k,2) &
               + stencil_eq2_u(i,WWW,j,k)      * vec(i-3,j,k,2) &
               + stencil_eq2_u(i,WW ,j,k)      * vec(i-2,j,k,2) &
               + stencil_eq2_u(i,W  ,j,k)      * vec(i-1,j,k,2) &
               + stencil_eq2_u(i,E  ,j,k)      * vec(i+1,j,k,2) &
               + stencil_eq2_u(i,EE ,j,k)      * vec(i+2,j,k,2) &
               + stencil_eq2_u(i,EEE,j,k)      * vec(i+3,j,k,2) &
               + stencil_eq2_u(i,SSS,j,k)      * vec(i,j-3,k,2) &
               + stencil_eq2_u(i,SS ,j,k)      * vec(i,j-2,k,2) &
               + stencil_eq2_u(i,S  ,j,k)      * vec(i,j-1,k,2) &
               + stencil_eq2_u(i,N  ,j,k)      * vec(i,j+1,k,2) &
               + stencil_eq2_u(i,NN ,j,k)      * vec(i,j+2,k,2) &
               + stencil_eq2_u(i,NNN,j,k)      * vec(i,j+3,k,2) &
               + stencil_eq2_v_interp(i,j,k)   * v_vec_on_u(i,j,k) ) 
            sum1(i) = sum1(i) + ( vec(i, j, k, 2)   * vec(i, j, k, 2) ) * fx(i)
         end do
         do pp=1,nphi-1
            do i = i0, in
               prod(i,j,k,2) = prod(i,j,k,2) + u(i,j,k,2, nphi-pp+1)*nu*V_aug(pp)
            end do
         end do
         do i = i0, in
            sum2(i) = sum2(i) + ( vis(i, j, k, 2) * prod(i, j, k, 2) ) * fx(i)
            sum3(i) = sum3(i) + ( vec(i, j, k, 2) * prod(i, j, k, 2) ) * fx(i)
         end do
      end do

      do j = j0, jn
      fyy = fy(j)
         do i = i0, in
            prod(i,j,k,3) = prod(i,j,k,3) &
              + stencil_eq3_v(i,P  ,j,k)      * vec(i  ,j,k,3) &
              + stencil_eq3_v(i,WWW,j,k)      * vec(i-3,j,k,3) &
              + stencil_eq3_v(i,WW ,j,k)      * vec(i-2,j,k,3) &
              + stencil_eq3_v(i,W  ,j,k)      * vec(i-1,j,k,3) &
              + stencil_eq3_v(i,E  ,j,k)      * vec(i+1,j,k,3) &
              + stencil_eq3_v(i,EE ,j,k)      * vec(i+2,j,k,3) &
              + stencil_eq3_v(i,EEE,j,k)      * vec(i+3,j,k,3)
            prod(i,j,k,1) = dt * ( prod(i,j,k,1) &
              + stencil_eq1_v(i,sss_stag,j,k) * vec(i,j-3,k,3) &
              + stencil_eq1_v(i,ss_stag ,j,k) * vec(i,j-2,k,3) &
              + stencil_eq1_v(i,s_stag  ,j,k) * vec(i,j-1,k,3) &
              + stencil_eq1_v(i,n_stag  ,j,k) * vec(i,j  ,k,3) &
              + stencil_eq1_v(i,nn_stag ,j,k) * vec(i,j+1,k,3) &
              + stencil_eq1_v(i,nnn_stag,j,k) * vec(i,j+2,k,3) ) 
            prod(i,j,k,3) = dt * ( prod(i,j,k,3) &
              + stencil_eq3_v(i,SSS,j,k)      * vec(i,j-3,k,3) &
              + stencil_eq3_v(i,SS ,j,k)      * vec(i,j-2,k,3) &
              + stencil_eq3_v(i,S  ,j,k)      * vec(i,j-1,k,3) &
              + stencil_eq3_v(i,N  ,j,k)      * vec(i,j+1,k,3) &
              + stencil_eq3_v(i,NN ,j,k)      * vec(i,j+2,k,3) &
              + stencil_eq3_v(i,NNN,j,k)      * vec(i,j+3,k,3) &
              + stencil_eq3_u_interp(i,j,k)   * u_vec_on_v(i,j,k) ) 
            sum1(i) = sum1(i) + ( vec(i, j, k, 3)   * vec(i, j, k, 3) ) * fyy
         end do
         do pp=1,nphi-1
            do i = i0, in
               prod(i,j,k,1) = prod(i,j,k,1) + u(i,j,k,1, nphi-pp+1)*nu*V_aug(pp)
               prod(i,j,k,3) = prod(i,j,k,3) + u(i,j,k,3, nphi-pp+1)*nu*V_aug(pp)
            end do
         end do
         do i = i0, in
            sum2(i) = sum2(i) + ( vis(i, j, k, 1) * prod(i, j, k, 1) ) &
                              + ( vis(i, j, k, 3) * prod(i, j, k, 3) ) * fyy
            sum3(i) = sum3(i) + ( vec(i, j, k, 1) * prod(i, j, k, 1) ) &
                              + ( vec(i, j, k, 3) * prod(i, j, k, 3) ) * fyy
         end do
      end do
   end do

100   continue

   do v = 4,nvars
      do k = 1,l_nk
	do slc= 1, 1 !ns
	  do j = j0, jn
	      do i = i0, in
		!****************************************
		! Compute Jacobian-vector for tracers   *
		!****************************************
		prod(i,j,k,v) = dt * (  stencil_adv_tr(i,P  ,j,k) * vec(i  ,j,k,v) &
				      + stencil_adv_tr(i,WWW,j,k) * vec(i-3,j,k,v) &
				      + stencil_adv_tr(i,WW ,j,k) * vec(i-2,j,k,v) &
				      + stencil_adv_tr(i,W  ,j,k) * vec(i-1,j,k,v) &
				      + stencil_adv_tr(i,E  ,j,k) * vec(i+1,j,k,v) &
				      + stencil_adv_tr(i,EE ,j,k) * vec(i+2,j,k,v) &
				      + stencil_adv_tr(i,EEE,j,k) * vec(i+3,j,k,v) &
				      + stencil_adv_tr(i,SSS,j,k) * vec(i,j-3,k,v) &
				      + stencil_adv_tr(i,SS ,j,k) * vec(i,j-2,k,v) &
				      + stencil_adv_tr(i,S  ,j,k) * vec(i,j-1,k,v) &
				      + stencil_adv_tr(i,N  ,j,k) * vec(i,j+1,k,v) &
				      + stencil_adv_tr(i,NN ,j,k) * vec(i,j+2,k,v) &
				      + stencil_adv_tr(i,NNN,j,k) * vec(i,j+3,k,v) ) 
		sum1(i) = sum1(i) + ( vec(i, j, k, v)   * vec(i, j, k, v) )
	      end do
	      do pp=1,nphi-1
		do i = i0, in
		    prod(i,j,k,v) = prod(i,j,k,v) + u(i,j,k,v, nphi-pp+1)*nu*V_aug(pp)
		end do
	      end do
	      do i = i0, in
		sum2(i) = sum2(i) + ( vis(i, j, k, v) * prod(i, j, k, v) )
		sum3(i) = sum3(i) + ( vec(i, j, k, v) * prod(i, j, k, v) )
	      end do
	  end do
	end do
      end do
   end do


   return
   end subroutine kernel_matvec2
