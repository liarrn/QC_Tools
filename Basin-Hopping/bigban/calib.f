            subroutine calib(elem,x,y,z,rms,
     $                 xcen,ycen,zcen)
            include "param.h"
            
            double precision x(100),y(100),z(100),rms
            double precision xcen,ycen,zcen
            character*2      elem(100)

            xcen=0.0
            ycen=0.0
            zcen=0.0
            rms=0.0

            do i=1,Natoms
               xcen=xcen+x(i)
               ycen=ycen+y(i)
               zcen=zcen+z(i)
            enddo

            xcen=xcen/Natoms
            ycen=ycen/Natoms
            zcen=zcen/Natoms

            do i=1,Natoms
              rms=rms+sqrt((x(i)-xcen)**2+(y(i)-ycen)**2+(z(i)-zcen)**2)
            enddo

            rms=rms/Natoms


            do i=1,Natoms
               x(i)=-xcen+x(i)
               y(i)=-ycen+y(i)
               z(i)=-zcen+z(i)
            enddo

            xcen=0.0
            ycen=0.0
            zcen=0.0

            return

            end subroutine



