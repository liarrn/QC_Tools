              subroutine  bigban(elem,x,y,z,
     $                          xcen,ycen,zcen)
              include "param.h"
           
              double precision x(100),y(100),z(100)
              double precision xf(100),yf(100),zf(100)
              double precision xcen,ycen,zcen
              double precision multi,rij
              character*2      elem(100)

          

100           do i=1, Natoms
                 multi=(1+rand())*Fact
!                 xf(i)=1/(x(i)-xcen)*multi+x(i)
!                 yf(i)=1/(y(i)-ycen)*multi+y(i)
!                 zf(i)=1/(z(i)-zcen)*multi+z(i)
                 xf(i)=1/(x(i)-xcen)*multi*((-1)**int(2*rand()))+x(i)
                 multi=(1+rand())*Fact
                 yf(i)=1/(y(i)-ycen)*multi*((-1)**int(2*rand()))+y(i)
                 multi=(1+rand())*Fact
                 zf(i)=1/(z(i)-zcen)*multi*((-1)**int(2*rand()))+z(i)
!                 xf(i)=x(i)+multi*((-1)**int(2*rand()))
!                 multi=(1+rand())*Fact
!                 yf(i)=y(i)+multi*((-1)**int(2*rand()))
!                 multi=(1+rand())*Fact
!                 zf(i)=z(i)+multi*((-1)**int(2*rand()))
                 do while(abs(xf(i)).gt.rdmax) 
                   xf(i)=xf(i)*(rand()*0.3+0.5)
!                   xf(i)=mod(xf(i),int(rand()*rdmax))
                 enddo
                 do while(abs(yf(i)).gt.rdmax) 
                   yf(i)=yf(i)*(rand()*0.3+0.5)
!                   yf(i)=mod(yf(i),int(rand()*rdmax))
                 enddo
                 do while(abs(zf(i)).gt.rdmax) 
                   zf(i)=zf(i)*(rand()*0.3+0.5)
!                   zf(i)=mod(zf(i),int(rand()*rdmax))
                 enddo
              enddo


              do i=1,Natoms-1
               do j=i+1,Natoms
                rij=(xf(i)-xf(j))**2+(yf(i)-yf(j))**2+(zf(i)-zf(j))**2
                rij=sqrt(rij)
                if(rij.lt.rdmin) then
                  goto 100
                endif
               end do
              end do

              do i=1,Natoms
                x(i)=xf(i)
                y(i)=yf(i)
                z(i)=zf(i)
              enddo

 
              return

              end subroutine

