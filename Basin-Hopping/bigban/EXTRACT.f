	subroutine EXTRACT(rxx,ryy,rzz,fp,conv)
	include 'param.h'
!        implicit double precision (a-h, o-z)
        character line*100, e1*2, e2*1
        logical conv
        PARAMETER(eunit=4.3597482e-18,runit=0.5291772083e-10,
     *            amunit=1.66053873e-27,bzcon=1.380603e-23,
     *            Avcon = 6.02214199e23)

         double precision rxx(100), ryy(100), rzz(100)
         double precision rxxf(100), ryyf(100), rzzf(100)
         double precision fp
         integer k1,k2,k3,k4,k5

!        ru1 = runit*1.0e10

          conv=.FALSE.

          open( 001, file = 'input.outmol' )

          do kkk = 1, MXLN

              read(001,'(a)') line

            k1 = index(line, 'Total E (au)')
            k2 = index(line, 'Input Coordinates (Angstroms)')
!            k2 = index(line, 'ATOMIC  COORDINATES (au)')
            k3 = index(line, 'Message: DMol3 job finished successfully')
            k4 = index(line, 'ENDREAD')

!           if (index(line, 'opt==').ne.0) read(line(15:30),*) fp
            if ( k1 .ne. 0 )  then
               read(001,'(a)') line
               read(001,'(a)') line
               read(line(5:23),*) fp
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


