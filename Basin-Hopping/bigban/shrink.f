               subroutine shrink(elem,x,y,z,
     $                          xcen,ycen,zcen)
               include "param.h"

               double precision x(100),y(100),z(100)
               double precision xf(100),yf(100),zf(100)
               double precision xcen,ycen,zcen
               double precision multi,dis,ratio,rms
               character*2      elem(100)
               integer istep

!               multi=rand()*SFact
               ratio=Sfact

               istep=0

100            do i=1,Natoms
                  multi=(1+rand())*ratio
                  xf(i)=x(i)-(x(i)-xcen)*multi
                  yf(i)=y(i)-(y(i)-ycen)*multi
                  zf(i)=z(i)-(z(i)-zcen)*multi
               end do


               do i=1,Natoms-1
                do j=i+1,Natoms
                 dis=(xf(i)-xf(j))**2+(yf(i)-yf(j))**2
     $                +(zf(i)-zf(j))**2 
                 dis=sqrt(dis)
!                 k=0
                 if(dis.le.rdmin) then 
                    x(j)=x(j)+((-1.0)**(int(rand()+0.5)))*slen*rand()
                    y(j)=y(j)+((-1.0)**(int(rand()+0.5)))*slen*rand()
                    z(j)=z(j)+((-1.0)**(int(rand()+0.5)))*slen*rand()
!                   istep=istep+1
!                   if(istep.gt.10) then
!                    call bigban(elem,x,y,z,xcen,ycen,zcen)
!                    call calib(elem,x,y,z,rms,xcen,ycen,zcen)
!                     ratio=-Sfact
!                   endif
!                  multi=multi*0.5
!                  goto 100
!                 dis=(xf(i)-xf(j))**2+(yf(i)-yf(j))**2
!     $                +(zf(i)-zf(j))**2 
!                 dis=sqrt(dis)
!                 k=k+1
!                 if(k.ge.200) goto 999
                 endif
                enddo
               enddo

               do i=1,Natoms
                x(i)=xf(i)
                y(i)=yf(i)
                z(i)=zf(i)
               enddo
   

               return

               end subroutine
               
              
