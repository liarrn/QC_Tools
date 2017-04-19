           subroutine initial_str(elem,x,y,z)
           include "param.h"
           double precision x(100),y(100),z(100)
           character*2  elem(100)

           open (unit=30, file='GEOMETRY.xyz', status='unknown')

           read (30, *)
           read (30, *)
           do i=1, Natoms
              read(30,*) elem(i),x(i),y(i),z(i)
           end do


           close(30)
 
           return

           end

