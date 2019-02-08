! /* 
!  * This is free software; you can redistribute it and/or
!  * modify it under the terms of the GNU Lesser General Public
!  * License as published by the Free Software Foundation,
!  * version 2.1 of the License.
!  *
!  * This is distributed in the hope that it will be useful,
!  * but WITHOUT ANY WARRANTY; without even the implied warranty of
!  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  * Lesser General Public License for more details.
!  *
!  * You should have received a copy of the GNU Lesser General Public
!  * License along with this library; if not, write to the
!  * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
!  * Boston, MA 02111-1307, USA.
!  */
program test
  use ISO_C_BINDING
  implicit none
  include "f90_interfaces.inc"
#if defined(USE_MPI)
  include "mpif.h"
#endif
#include "test.inc"

  if(my_proc == 0) print *,'======================== 3 variables ==============================='

  call tricublin_zyx3_n(d3,f123(1,1,1,1),pxpypz,lv,NI*NJ*NK)
#if defined(USE_MPI)
  call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
  t0 = nanocycles()
  call tricublin_zyx3_n(d3,f123(1,1,1,1),pxpypz,lv,NI*NJ*NK)
  t1 = nanocycles()
  if(my_proc == 0) print *,"nanoseconds per point =",(t1-t0)/(NI*NJ*NK*3)

  exact = 0
  delta = 0.0
  delta1 = 0.0
  delta2 = 0.0
  delta3 = 0.0
  avg = 0.0
  avg1 = 0.0
  avg2 = 0.0
  avg3 = 0.0
  do k = 1, NK
    do j = 1, NJ
      do i = 1, NI
        error = abs(expected(i,j,k)+00.0 - d3(1,i,j,k))
        if(error == 0.0) exact = exact + 1
        delta1 = max(delta1,error / expected(i,j,k))
        avg1 = avg1 + (error / expected(i,j,k))

        error = abs(expected(i,j,k)+10.0 - d3(2,i,j,k))
        if(error == 0.0) exact = exact + 1
        delta2 = max(delta2,error / (expected(i,j,k)+10.0))
        avg2 = avg2 + error / (expected(i,j,k)+10.0)

        error = abs(expected(i,j,k)+20.0 - d3(3,i,j,k))
        if(error == 0.0) exact = exact + 1
        delta3 = max(delta3,error / (expected(i,j,k)+20.0))
        avg3 = avg3 + error / (expected(i,j,k)+20.0)

        delta = max(delta1,delta2,delta3)
        if((delta > .000001 .or. i+j+k == 3) .and. (my_proc == 0)) then
          print *,'i,j,k',i,j,k
          print *,'pxpypz =',pxpypz(1,i,j,k),pxpypz(2,i,j,k),pxpypz(3,i,j,k)
          print *,'expected =',expected(i,j,k) + [0.0, 10.0, 20.0]
          print *,'result   =',d3(:,i,j,k)
!           if(delta > .000001) stop
          if(delta > .000001) goto 222
        endif
      enddo
    enddo
  enddo
  if(my_proc == 0) then
  print *,'exact =',exact,' out of',NI*NJ*NK*3
  print *,'%      ',real(exact)/real(NI*NJ*NK*3)*100
  print*,'maxerr =',delta1, delta2, delta3
  print*,'avgerr =',real([avg1,avg2,avg3] / (NI*NJ*NK))
  endif
222 continue
#if defined(USE_MPI)
  call mpi_finalize(ierr)
#endif
end program test
