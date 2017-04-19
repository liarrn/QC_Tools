             program  main
             include "param.h"
             
             double precision x(100),y(100),z(100)
             double precision xf(MXCY,100),yf(MXCY,100)
             double precision zf(MXCY,100)
             double precision energy,en(MXCY)
             double precision en1,eni,xfi,yfi,zfi
             double precision RMS,rms1,multi
             double precision xcen,ycen,zcen
             integer step
             character*2 elem(100)
             logical  conv

             do i=1,MXCY
               do j=1,100
                 xf(i,j)=0.0d0
                 yf(i,j)=0.0d0
                 zf(i,j)=0.0d0
               enddo
               en(i)=0.0d0
             enddo
             energy=0.0d0
        
             open(unit=01,file='END',status='unknown')
             write(01,*) 'ENDREAD'
             close(01)



             call initial_str(elem,x,y,z)
    

!             call quench_opt(elem,x,y,z,energy,conv)
                call IODMOL(x,y,z,energy,conv)
             
             call calib(elem,x,y,z,RMS,xcen,ycen,zcen)

             en(1)=energy
             do j=1,Natoms
               xf(1,j)=x(j)
               yf(1,j)=y(j)
               zf(1,j)=z(j)
             enddo
!            write(10,*) "1",en(1)
!            do i=1,Natoms
!               write(10,*) elem(1),xf(1,i),yf(1,i),zf(1,i)
!            enddo
!            write(10,*) "1",xcen,ycen,zcen,RMS

  
            do i=2,MXCY
 
             call bigban(elem,x,y,z,xcen,ycen,zcen)

             call calib(elem,x,y,z,rms1,xcen,ycen,zcen)

             step=0

 
             open(unit=02,file='test',status='unknown')
             write(02,*) i,"begin"
             do j=1,Natoms
                write(02,*) x(j),y(j),z(j)
             enddo
             close(02)

             do while((abs(rms1-RMS).gt.RAD)
     $                .OR.(conv.EQV..FALSE.)) 
                call shrink(elem,x,y,z,xcen,ycen,zcen)

             step=step+1
             open(unit=01,file='initial',status='unknown')
             write(01,*) "after shrink",i,step
             do j=1,Natoms
                write(01,*) x(j),y(j),z(j)
             enddo
             write(01,*) "center",xcen,ycen,zcen
             write(01,*)
             close(01)


!                call quench_opt(elem,x,y,z,energy,conv)
	        call IODMOL(x,y,z,energy,conv)
!                do while(conv.EQV..FALSE.) 
!                    call shrink(elem,x,y,z,xcen,ycen,zcen)
!                    call quench_opt(elem,x,y,z,energy,conv)
!                enddo

                call calib(elem,x,y,z,rms1,xcen,ycen,zcen)
 
             enddo
             
             
             en(i)=energy
             do j=1,Natoms
               xf(i,j)=x(j)
               yf(i,j)=y(j)
               zf(i,j)=z(j)
             enddo
!             RMS=rms1

            open(unit=10, file='out',status='unknown')
            do j=1,i
              write(10,*) j,en(j)
              do k=1,Natoms
                write(10,*) elem(k),xf(j,k),yf(j,k),zf(j,k)
              enddo
              write(10,*)
            enddo
            close(10)
!            stop

            enddo

            do i=1,MXCY-1
              do j=i+1,MXCY
                if(en(i).gt.en(j)) then
                  eni=en(j)
                  en(j)=en(i)
                  en(i)=eni
                  do k=1,Natoms
                    xfi=xf(j,k)
                    xf(j,k)=xf(i,k)
                    xf(i,k)=xfi
                    yfi=yf(j,k)
                    yf(j,k)=yf(i,k)
                    yf(i,k)=yfi
                    zfi=zf(j,k)
                    zf(j,k)=zf(i,k)
                    zf(i,k)=zfi
                  enddo
                endif
              enddo
            enddo

            open(20,file='Final',status='unknown')
            write(20,*) "Final results"
            do i=1,MXCY
              write(20,'(I4,f14.6)') i, en(i)
                do j=1,Natoms
                  write(20,*) elem(j),xf(i,j),yf(i,j),zf(i,j)
                enddo
            enddo

             close(20)

             end


