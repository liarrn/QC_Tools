! This program is written by Yi Gao on 2014/7/25.
! It is a program for searching for the best structures of Basin-Hopping
! Supported Cluster (BHSC). 
! The largest system is 500 atoms in total.
             program  main
             include "param.h"
             
             double precision x(500),y(500),z(500)
             double precision xf(MXCY,500),yf(MXCY,500)
             double precision zf(MXCY,500)
             double precision energy,en(MXCY)
             double precision en1,eni,xfi,yfi,zfi
             double precision RMS,rms1,multi
             double precision xcen,ycen,zcen
             integer step
             character*2 elem(500)
             logical  conv
             do i=1,MXCY
               do j=1,500
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
             call IODMOL(elem,x,y,z,energy,conv)
             
             en(1)=energy
             n=1
             
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
 
             call move(elem,x,y,z)
             
             step=0
             open(unit=02,file='test',status='unknown')
             write(02,*) i,"begin"
             do j=1,Natoms
                write(02,*) x(j),y(j),z(j)
             enddo
             close(02)
                call IODMOL(elem,x,y,z,energy,conv)
              m=i
              
              en(m)=energy
              do j=1,Natoms
               xf(m,j)=x(j)
               yf(m,j)=y(j)
               zf(m,j)=z(j)
              enddo
                           
              ratio=exp(-(en(m)-en(n))*2625499.62/(8.3144621*TEMP))
              
              if(en(m).gt.en(n).and.ratio.le.rand()) then
                do j=1,Natoms
                  x(j)=xf(n,j)
                  y(j)=yf(n,j)
                  z(j)=zf(n,j)
                enddo
              else
                n=m                
              endif
 
             
           do m=1,i-1
              do j=m+1,i
                if(en(m).gt.en(j)) then
                  eni=en(j)
                  en(j)=en(m)
                  en(m)=eni
                  do k=1,Natoms
                    xfi=xf(j,k)
                    xf(j,k)=xf(m,k)
                    xf(m,k)=xfi
                    yfi=yf(j,k)
                    yf(j,k)=yf(m,k)
                    yf(m,k)=yfi
                    zfi=zf(j,k)
                    zf(j,k)=zf(m,k)
                    zf(m,k)=zfi
                  enddo
                endif
              enddo
            enddo
             
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
! This subroutine is written by Yi Gao on 2014/7/25.
! This subroutine is to read the initial structure of the supported
! cluster.
! The largest size of the system is 500 atoms in total.
           subroutine initial_str(elem,x,y,z)
           include "param.h"
           double precision x(500),y(500),z(500)
           character*2  elem(500)
           open (unit=30, file='GEOMETRY.xyz', status='unknown')
           read (30, *)
           read (30, *) 
           do i=1, Natoms
              read(30,*) elem(i),x(i),y(i),z(i)
           end do
           close(30)
 
           return
           end
!  This subroutine is written by Yi Gao on 2014/7/25.
!  This subroutine is to generate the DMol3 input file.
!  The maximum is 500 atoms in total.
        subroutine IODMOL(elem,rxx,ryy,rzz,fp,conv)
        include "param.h"
!        implicit double precision (a-h, o-z)
        PARAMETER(eunit=4.3597482e-18,runit=0.5291772083e-10,
     *            amunit=1.66053873e-27,bzcon=1.380603e-23,
     *            Avcon = 6.02214199e23)
        logical task, conv
        double precision rxx(500), ryy(500), rzz(500)
        double precision fp
        real*8 leng
        character*2 elem(500)
          x=11.8360
          y=12.9938
          z=18.4771
          A=90.0000
          OPEN  (02, FILE = 'input.car')
          WRITE (02, '(a)')  '!BIOSYM archive 3'
          WRITE (02, '(a)')  'PBC=OFF'
          WRITE (02, '(a)')  'Materials Studio Generated CAR File'
          WRITE (02, '(a)')  '!DATE Wed Nov 29 13:31:15 2006'
!          WRITE (02, '(A,6f10.4,A)')  'PBC',x,y,z,A,A,A,' (P1)' 
          do ia = 1, Ncluster
            WRITE (02, '(A,3X,3f15.9,21X,A,f8.3)') elem(ia), 
     +             rxx(ia), ryy(ia), rzz(ia),elem(ia),0.000
          enddo
          do ia = Ncluster+1, Ncluster+Nrelax 
            WRITE (02, '(A,3X,3f15.9,21X,A,f8.3)') elem(ia), 
     +             rxx(ia), ryy(ia), rzz(ia),elem(ia),0.000
          enddo
          do ia = Ncluster+Nrelax+1, Natoms 
            WRITE (02, '(A,3X,3f15.9,21X,A,f8.3)') elem(ia), 
     +             rxx(ia), ryy(ia), rzz(ia),elem(ia),0.000
          enddo
          WRITE (02,'(a)')  'end'
          WRITE (02,'(a)')  'end'
          CLOSE (02)
          call system('rm tmpp')          
                  call system('./*.sh -np 24 input ')
          call system('cat input.outmol END > tmpp')
          call system('cp tmpp input.outmol')
          call EXTRACT(rxx,ryy,rzz,fp,conv)
        return
        end
! This subroutine is written by Yi Gao on 2014/7/25.
! This subroutine is to read the initial structure of the supported
! cluster.
! The largest size of the system is 500 atoms in total.
           subroutine move(elem,x,y,z)
           include "param.h"
           double precision x(500),y(500),z(500)
           character*2  elem(500)
           if(NSYS.eq.1) then
!  Pure cluster system
            
             call move_cluster(elem,x,y,z)
           elseif(NSYS.eq.2) then
!  Exchange atom between the cluster and the substrate
             call exchange_cs(elem,x,y,z)
           elseif(NSYS.eq.3) then
       
!  Cluster on Substrate
             call move_cluster_on_sub(elem,x,y,z)
           
           endif
           return
           end
!  This subroutine is written by Yi Gao on 2014/8/7.
!  This is a subroutine to move the atom of the cluster unbiased.    
!  The maximum number is 500 atoms in total.
              subroutine move_cluster(elem,x,y,z)
              include "param.h"
               
              double precision x(500),y(500),z(500)                
              double precision xx,yx,zx,r(500)             
              double precision multi,rij            
              double precision xcen_i,ycen_i,zcen_i
              double precision xcen_f,ycen_f,zcen_f
              double precision r_ave_i,r_ave_f
              character*2      elem(500)    
              integer  n_far,nf,nff
              
              xcen_i=0.0
              ycen_i=0.0
              zcen_i=0.0
              r_ave_i=0.0
              r_ave_f=0.0
              
              do i=1,Ncluster
                xcen_i=xcen_i+x(i)
                ycen_i=ycen_i+y(i)
                zcen_i=zcen_i+z(i)
              end do
              
              xcen_i=xcen_i/Ncluster
              ycen_i=ycen_i/Ncluster
              zcen_i=zcen_i/Ncluster
              
              
               do i=1, Ncluster
                r(i)=(x(i)-xcen_i)**2+(y(i)-ycen_i)**2+(z(i)-zcen_i)**2
                r(i)=sqrt(r(i))
                r_ave_i=r_ave_i+r(i)
               end do
              
               r_ave_i=r_ave_i/Ncluster
               r_ave_f=r_ave_i
              
              do while(abs(r_ave_f-r_ave_i).le.0.1) 
               r_ave_i=r_ave_f
               n_far=1
               do i=1, Ncluster-1
                do j=i+1, Ncluster
                  if(r(i).lt.r(j)) then
                   if(r(j).gt.r(n_far)) n_far=j
                  else
                   if(r(i).gt.r(n_far)) n_far=i
                  endif
                end do
               end do
 
              
               nf=0
               do while(nf.eq.0.or.nff.eq.0)
                nf=1
                nff=0
                fi=2*pi*rand()
                psi=2*pi*rand()
                              
                xx=xcen_i+r(n_far)*(1+(-1)**int(0.5+rand())*0.5*rand())
     $             *cos(fi)*cos(psi)
                yy=ycen_i+r(n_far)*(1+(-1)**int(0.5+rand())*0.5*rand())
     $             *cos(fi)*sin(psi)
                zz=zcen_i+r(n_far)*(1+(-1)**int(0.5+rand())*0.5*rand())
     $             *sin(fi)
              
                do i=1,Ncluster
                  if(i.ne.n_far) then
                    rij=(x(i)-xx)**2+(y(i)-yy)**2+(z(i)-zz)**2
                    rij=sqrt(rij)
                    if(rij.lt.rdmin) then
                      nf=0
                      goto 100
                    endif
                    if(rij.ge.rdmin.and.rij.le.rdmax) then
                      nff=1
                    endif
                  endif
                  
                end do
                
100            end do
               x(n_far)=xx
               y(n_far)=yy
               z(n_far)=zz
               xcen_f=0.0
               ycen_f=0.0
               zcen_f=0.0
               
               do i=1,Ncluster
                xcen_f=xcen_f+x(i)
                ycen_f=ycen_f+y(i)
                zcen_f=zcen_f+z(i)
               end do
              
               xcen_f=xcen_f/Ncluster
               ycen_f=ycen_f/Ncluster
               zcen_f=zcen_f/Ncluster              
               
               r_ave_f=0.0
               do i=1, Ncluster
                r(i)=(x(i)-xcen_f)**2+(y(i)-ycen_f)**2+(z(i)-zcen_f)**2
                r(i)=sqrt(r(i))
                r_ave_f=r_ave_f+r(i)
               end do  
               r_ave_f=r_ave_f/Ncluster
              
              end do             
                           
              end subroutine
        subroutine EXTRACT(rxx,ryy,rzz,fp,conv)
              include 'param.h'
!        implicit double precision (a-h, o-z)
        character line*100, e1*2, e2*1
        logical conv
        PARAMETER(eunit=4.3597482e-18,runit=0.5291772083e-10,
     *            amunit=1.66053873e-27,bzcon=1.380603e-23,
     *            Avcon = 6.02214199e23)
         double precision rxx(500), ryy(500), rzz(500)
         double precision rxxf(500), ryyf(500), rzzf(500)
         double precision fp
         integer k1,k2,k3,k4,k5
          conv=.FALSE.
          open(001, file = 'input.outmol' )
          do kkk = 1, MXLN
              read(001,'(a)') line
            k1 = index(line, 'Ef')
            k2 = index(line, 'Final Coordinates (Angstroms)')
            k3 = index(line, 'Message: DMol3 job finished successfully')
            k4 = index(line, 'ENDREAD')
            if ( k1 .ne. 0 )  then
               read(line(9:22),*) fp
            endif
            if ( k2 .ne. 0 ) then
               read(001,'(a)') line
               read(001,'(a)') line
               do i = 1, Natoms
                  read(001,'(a)') line
                  read(line(20:80),*) rxxf(i), ryyf(i), rzzf(i)
!                  read(line(10:45),*) rxx(i), ryy(i), rzz(i)
               enddo
            endif
            if ( k3 .ne. 0 ) then
               conv = .TRUE.
               do i=1, Natoms
                 rxx(i)=rxxf(i)
                 ryy(i)=ryyf(i)
                 rzz(i)=rzzf(i)
               enddo
               goto 100
            endif
            if ((index(line, 'Error') .ne. 0) .or.( k4 .ne. 0)) then
               conv = .FALSE.
               goto 100
            endif
          enddo
100       continue
        close( 001 )
        return
        end

