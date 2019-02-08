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

  if(my_proc == 0) print *,'======================== 1 variable ================================'

  call tricublin_zyx1_n(d,f1(1,1,1),pxpypz,lv,NI*NJ*NK)
#if defined(USE_MPI)
  call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
  t0 = nanocycles()
  call tricublin_zyx1_n(d,f1(1,1,1),pxpypz,lv,NI*NJ*NK)
  t1 = nanocycles()
  if(my_proc == 0) print *,"nanoseconds per point =",(t1-t0)/(NI*NJ*NK)

  exact = 0
  delta = 0.0
  avg = 0.0
  do k = 1, NK
    do j = 1, NJ
      do i = 1, NI
        error = abs(expected(i,j,k) - d(i,j,k))
        if(error == 0.0) exact = exact + 1
        delta = max(delta,error / expected(i,j,k))
        avg = avg + (error / expected(i,j,k))
        if((delta > .000001 .or. i+j+k == 3) .and. (my_proc == 0)) then
          print *,'i,j,k',i,j,k
          print *,'pxpypz =',pxpypz(1,i,j,k),pxpypz(2,i,j,k),pxpypz(3,i,j,k)
          print *,'expected =',expected(i,j,k)
          print *,'result   =',d(i,j,k)
!           if(delta > .000001) stop
          if(delta > .000001) goto 111
        endif
      enddo
    enddo
  enddo
  if(my_proc == 0) then
  print *,'exact =',exact,' out of',NI*NJ*NK
  print *,'%      ',real(exact)/real(NI*NJ*NK)*100
  print*,'maxerr =',delta
  print*,'avgerr =',real(avg/(NI*NJ*NK))
  endif
111 continue
#if defined(USE_MPI)
  call mpi_finalize(ierr)
#endif
end program test
