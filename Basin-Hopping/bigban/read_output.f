          subroutine read_output(elem,x,y,z,energy,conv)
          include "param.h"

          double precision a(100),b(100),x(100),y(100),z(100)
          double precision xi(100),yi(100),zi(100)
          double precision energy,eni
          character*2      elem(100)
          character line*100
          character SCF1*3,SCF2*5,SCF3*6,SCF4,test*6
          logical conv
             
          open(unit=30, file='exam1.log', status='unknown')

          eni=energy   
          do i=1,Natoms
             xi(i)=x(i)
             yi(i)=y(i)
             zi(i)=z(i)
          enddo

          do i=1,MXLN
             read(30,'(a)') line
!             if(opt1.EQ."Convergence ".and.opt2.EQ."failure".OR.
!     $          opt1.EQ."Error       ".and.opt2.EQ."termina") then
             k1=index(line,'SCF Done:')
             k2=index(line,'=')
             k3=index(line,'Optimization completed')
             k5=index(line,'Error ')

             if((k1.ne.0).and.(k2.ne.0)) then
               read(line((k2+1):(k2+18)),*) energy
               write(*,*) energy
             endif

             if(k3.ne.0) then
               conv=.TRUE.
               do j=1, MXLN
                read(30,'(a)') line
                k4=index(line,'Input orientation:')
                if(k4.ne.0) then
                 read(30,'(a)') line
                 read(30,'(a)') line
                 read(30,'(a)') line
                 read(30,'(a)') line
                 do k=1,Natoms
                  read(30,'(a)') a(i),elem(i),b(i),x(i),y(i),z(i)
                 enddo
                 goto 100
                endif
               enddo 
             endif

             if(k5.ne.0) then
               conv=.FALSE.
               goto 200
             endif

          end do

200       open(unit=20, file='error', status='unknown')
          write(20,*) "Converge not filled"
          energy=eni
          do i=1,Natoms
              x(i)=xi(i)
              y(i)=yi(i)
              z(i)=zi(i)
          end do
          goto 99
          

100      open(unit=90,file='Energy',status='unknown')
          write(90,*) "energy",energy

          close(30)
          close(20)
          close(90)
         
99          return
          
          end

           
