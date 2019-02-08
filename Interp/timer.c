/* 
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#include <stdlib.h>
#include <stdint.h>
#include <sys/time.h>
uint64_t Nanocycles(void) {
  struct timeval t;
  uint64_t i;
  int status;
  status = gettimeofday(&t,NULL);
  i = t.tv_sec;
  i = i * 1000000;
  i = i + t.tv_usec;
  i = i * 1000;      // nanoseconds
  return i;
}
