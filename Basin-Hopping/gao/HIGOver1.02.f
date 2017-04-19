! Hunting Isomers with Global Optimization (HIGO) version 1.0
! This program is written by Yi Gao on 2014/7/25 and updated by Beien ZHU for MS6.0 & MS7.0.
! It is a program for searching for the best structures of cluster by DFT-Bashing hoping method
! It can be used for globe optimization of isolated cluster and supported Cluster as well (BHSC). 
! The largest system is 500 atoms in total.
             program  main
             include "param.h"
             
             integer i, j, k
             double precision x(500),y(500),z(500)
             double precision xf(MXCY,500),yf(MXCY,500)
             double precision zf(MXCY,500)
             double precision energy,en(MXCY)
             double precision en1,eni,xfi,yfi,zfi
             double precision RMS,rms1,multi
             double precision xcen,ycen,zcen
             integer step
             integer CNACount(MXCY, 1000)
             integer CNATmp(1000)
             character*2 elem(500)
             logical  conv
             integer idenList(100)
             integer p, q
             logical cnaAlike
             
             
             do i = 1, MXCY
                  do j = 1, 1000
                        CNACount(i, j) = 0
                  end do
             end do
             
             
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
             
             
             if (CNA_TOGGLE .eq. 1) then
                  call CNA_CALC(elem, x, y, z, CNATmp)
                  call CNA_ADD(1, CNATmp, CNACount)
            open(unit=12, file='cna',status='unknown')
            do i = 1, 1000
                  if (CNACount(1, i) .ne. 0) then
                        write(12, *) i, CNACount(i, j)
                  end if
            end do
            close(12)
            end if
            
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
            cnaAlike = .FALSE.
             call move(elem,x,y,z, CNACount)
             
             step=0
             open(unit=02,file='test',status='unknown')
             write(02,*) i,"begin"
             do j=1,Natoms
                write(02,*) x(j),y(j),z(j)
             enddo
             close(02)
                call IODMOL(elem,x,y,z,energy,conv)
              
              m=i
            if (CNA_TOGGLE .eq. 1) then
                  cnaAlike = .FALSE.
                  p = 1
                  q = 1
                  call CNA_CALC(elem, x, y, z, CNATmp)
                  call CNA_CHECK(CNATmp, CNACount, idenList, idenNum)
                  call CNA_ADD(m, CNATmp, CNACount)
                  if (idenNum .eq. 0) then
                        cnaAlike = .FALSE.
                  else
                        do p = 1, idenNum
                    if (abs(en(idenList(p))-energy) .lt. En_Cut) then
                              if (rand() < CNAAcpR) then
                              cnaAlike = .TRUE.
                              else
                              cnaAlike = .FALSE.
                              end if
                              EXIT
                    end if
                    end do
                  end if
            end if
              
              
              
              en(m)=energy
              do j=1,Natoms
               xf(m,j)=x(j)
               yf(m,j)=y(j)
               zf(m,j)=z(j)
              enddo
                           
              ratio=exp(-(en(m)-en(n))*2625499.62/(8.3144621*TEMP))
              
C               if(en(m).gt.en(n).and.ratio.le.rand()) then
C                 do j=1,Natoms
C                   x(j)=xf(n,j)
C                   y(j)=yf(n,j)
C                   z(j)=zf(n,j)
C                 enddo
C               else
C                 n=m                
C               endif
 
            if(en(m).lt.en(n)) then
                  n = m
        else if (rand() .lt. ratio .and. cnaAlike .eqv. .FALSE.) then
                  n = m
            else
                  do j=1,Natoms
                        x(j)=xf(n,j)
                        y(j)=yf(n,j)
                        z(j)=zf(n,j)
                  end do
             end if
             
             
             
            open(unit=10, file='out',status='unknown')
            do j=1,i
                  write(10, *) Natoms
              write(10,*) j,en(j)
              do k=1,Natoms
                write(10,*) elem(k),xf(j,k),yf(j,k),zf(j,k)
              enddo
            enddo
            close(10)
            
            if (CNA_TOGGLE .eq. 1) then 
            open(unit=12, file='cna',status='unknown')
            do j=1,i
              write(12,*) j,en(j)
              do k=1,1000
                  if (CNACount(j, k) .ne. 0) then
                        write(12,*) k, CNACount(j, k)
                  end if
              enddo
              write(12,*)
            enddo
            close(12)
            end if
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
          if (PBC_TOGGLE .eq. 0) then 
		    WRITE (02, '(a)')  'PBC=OFF'
          else if (PBC_TOGGLE .eq. 1) then
		    WRITE (02, '(a)')  'PBC=ON'
		    WRITE (02, '(a)')  'PBC   ', A_AXIS, B_AXIS, C_AXIS, ALPHA_ANG, 
     + BETA_ANG, GAMMA_ANG, '(P1)'
          end if
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
                  call system('./*.sh -np 12 input ')
          call system('cat input.outmol END > tmpp')
          call system('cp tmpp input.outmol')
          call EXTRACT(rxx,ryy,rzz,fp,conv)
        return
        end
! This subroutine is written by Yi Gao on 2014/7/25.
! This subroutine is to read the initial structure of the supported
! cluster.
! The largest size of the system is 500 atoms in total.
           subroutine move(elem,x,y,z, CNACount)
           include "param.h"
           double precision x(500),y(500),z(500)
           character*2  elem(500)
           integer CNACount(MXCY, 1000)
           
           if(NSYS.eq.1) then
!  Pure cluster system
            
             call move_cluster(elem,x,y,z, CNACount)
C            elseif(NSYS.eq.2) then
!  Exchange atom between the cluster and the substrate
C              call exchange_cs(elem,x,y,z)
C            elseif(NSYS.eq.3) then
       
!  Cluster on Substrate
C              call move_cluster_on_sub(elem,x,y,z, CNACount)
           
           endif
           return
           end
!  This subroutine is written by Yi Gao on 2014/8/7.
!  This is a subroutine to move the atom of the cluster unbiased.    
!  The maximum number is 500 atoms in total.
              subroutine move_cluster(elem,x,y,z, CNACount)
              include "param.h"
               
              double precision x(500),y(500),z(500)                
              double precision xx,yx,zx,r(500)             
              double precision multi,rij            
              double precision xcen_i,ycen_i,zcen_i
              double precision xcen_f,ycen_f,zcen_f
              double precision r_ave_i,r_ave_f
              character*2      elem(500)    
              integer  n_far,nf,nff, cnaValidation
              integer CNATmp(1000)
              integer CNACount(MXCY, 1000)
              integer idenNum, idenList(100)
              double precision xo(500),yo(500),zo(500)
              integer CNAToggle
              
              CNAToggle = CNA_TOGGLE
              CNAToggle = 0
              
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
 
 

C				some changes made here 
               nf=0
               cnaValidation = 0
C                do while(nf.eq.0.or.nff.eq.0)
      do while((nf.eq.0).or.(nff.eq.0).or.(cnaValidation.eq.0))
                nf=1
                cnaValidation = 1
                nff=0
                fi=2*pi*rand()
                psi=2*pi*rand()
                              
                xx=xcen_i+r(n_far)*(1+(-1)**int(0.5+rand())*0.5*rand())
     $             *cos(fi)*cos(psi)
                yy=ycen_i+r(n_far)*(1+(-1)**int(0.5+rand())*0.5*rand())
     $             *cos(fi)*sin(psi)
                zz=zcen_i+r(n_far)*(1+(-1)**int(0.5+rand())*0.5*rand())
     $             *sin(fi)
              
              if (CNAToggle .eq. 1) then
              do i = 1, Ncluster
                  if (i .eq. n_far) then
                        xo(i) = xx
                        yo(i) = yy
                        zo(i) = zz
                  else
                        xo(i) = x(i)
                        yo(i) = y(i)
                        zo(i) = z(i)
                  end if
              end do
                  call CNA_CALC(elem, xo, yo, zo, CNATmp)
                  call CNA_CHECK(CNATmp, CNACount, idenList, idenNum)
                  if (idenNum .eq. 0) then
                  cnaValidation = 1
                  else
                  cnaValidation = 0
                  end if
              else if (CNAToggle .eq. 0) then
              cnaValidation = 1
              end if
              
              
                do i=1,Ncluster
                  if(i.ne.n_far) then
                    rij=(x(i)-xx)**2+(y(i)-yy)**2+(z(i)-zz)**2
                    rij=sqrt(rij)
C                     if(rij.lt.rdmin) then
C                       nf=0
C                       goto 100
C                     endif
C                     if(rij.ge.rdmin.and.rij.le.rdmax) then
                    if(rij.lt.rdmin) then
                      nf=0
                      EXIT
                    endif
                    if(rij.ge.rdmin.and.rij.le.rdmax) then
                      nff=1
                    endif
                  endif
                end do
C                 print *, nf, cnaValidation
                end do
                
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

        

      subroutine CNA_ADD(m, CNATmp, CNACount)
            include "param.h"
            integer CNACount(MXCY, 1000)
            integer CNATmp(1000)
            integer m
            integer i, j, k
            
            do i = 1, 1000
                  CNACount(m, i) = CNATmp(i)
            end do
      end 
        
        

      subroutine CNA_CHECK(CNATmp, CNACount, idenList, idenNum)
            include "param.h"
            integer CNACount(MXCY, 1000)
            integer CNATmp(1000)
            integer i, j, idenNum
            integer CNADiff
            integer idenList(100)
           
            idenNum = 0
            do i = 1, 100
                  idenList(i) = 0
            end do
            do i = 1, MXCY
                  CNADiff = 0
                  do j = 1, 1000
                  CNADiff = CNADiff + abs(CNACount(i, j) - CNATmp(j))
                  end do
                  if (CNADiff .lt. CNACut) then
                        idenNum = idenNum + 1
                        idenList(idenNum) = i
                  end if
            end do
           
           
        end

		
      subroutine CNA_CALC(elem, x, y, z, CNATmp)
            include "param.h"
            integer CNATmp(1000)
            double precision x(500),y(500),z(500)
            character*2  elem(500)
            integer i, j, k
            integer nij, ix1, ix2, ix3
            integer set(1000), Aset(1000, 2), setm, setn
            real*8 rxij, ryij, rzij, rij, rxmn, rymn, rzmn, rmn
            integer maxlength, m, n, p, q, a, b, flag, c
            integer SubBond(1000), nBond(1000)
            integer CNAType
            
            do i = 1, 1000
                  CNATmp(i) = 0
            end do
            nij = 0
            do i = 1, Natoms-1
                  do j = i + 1, Natoms
                        rxij=x(i)-x(j)
                        ryij=y(i)-y(j)
                        rzij=z(i)-z(j)
                        rij=rxij**2+ryij**2+rzij**2
                        rij=sqrt(rij)
                        if (rij .gt. CUTOFF_DIS) then
                              CYCLE
                        end if
                        nij=nij+1
                        ix1=0
                        ix2=0
                        ix3=0
                        
C                         the first index for CNA
                        do k=1,Natoms
                              if (k.ne.i .and. k.ne.j) then
                                    rxik=x(i)-x(k)
                                    ryik=y(i)-y(k)
                                    rzik=z(i)-z(k)
                                    rik=rxik**2+ryik**2+rzik**2
                                    rik=sqrt(rik)
                                    rxjk=x(j)-x(k)
                                    ryjk=y(j)-y(k)
                                    rzjk=z(j)-z(k)
                                    rjk=rxjk**2+ryjk**2+rzjk**2
                                    rjk=sqrt(rjk)
                   if (rik.lt.CUTOFF_DIS .and. rjk.lt.CUTOFF_DIS) then
                                          ix1=ix1+1
                                          set(ix1)=k
                                    end if
                              end if
                        end do
                        
C                         the second index for CNA
                        do m = 1, ix1-1
                              do n = m + 1, ix1
                              setm=set(m)
                              setn=set(n)
                              rxmn=x(setm)-x(setn)
                              rymn=y(setm)-y(setn)
                              rzmn=z(setm)-z(setn)
                              rmn=rxmn**2+rymn**2+rzmn**2
                              rmn=sqrt(rmn)
                              if (rmn.lt.CUTOFF_DIS) then
                                    ix2=ix2+1
                                    Aset(ix2,1)=setm
                                    Aset(ix2,2)=setn
                              end if
                              end do
                        end do
                        
                        
C                         the third index for CNA
            maxlength = 0
            do p = 1, Natoms
                  SubBond(p) = p
                        nBond(p) = 0
            end do
            do p = 1, ix2
                  m = Aset(p, 1)
                        n = Aset(p, 2)
                        a = nBond(m)
                        b= nBond(n)
                        c = SubBond(m)
                  if ( SubBond(m) .eq. SubBond(n) ) then
                              do q = 1, Natoms
                                    if (SubBond(q) .eq. c) then
                                                nBond(q) = nBond(q) + 1
                                          end if
                              end do
                  else
                              do q = 1, Natoms
                              if (SubBond(q) .eq. c) then
                                    SubBond(q) = SubBond(n)
                              end if
                              if (SubBond(q) .eq. SubBond(n)) then
                                    nBond(q) = a + b + 1
                              end if
                        end do
                  end if
            end do
            do p = 1, Natoms
                  if (maxlength .lt. nBond(p)) then
                        maxlength = nBond(p)
                  end if
            end do 
            ix3 = maxlength
                        
                        
                        CNAType=ix1*100+ix2*10+ix3
                        if (CNAType >= 1000) then
                              CNAType = 999
                        end if
                        if(CNAType .ne. 0) then
                              CNATmp(CNAType)=CNATmp(CNAType)+1
                        end if
                        
                        
                        
                  end do
            end do
            
            
        end