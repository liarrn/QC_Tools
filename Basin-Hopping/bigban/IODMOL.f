	subroutine IODMOL(rxx,ryy,rzz,fp,conv)
	include "param.h"
!        implicit double precision (a-h, o-z)
        PARAMETER(eunit=4.3597482e-18,runit=0.5291772083e-10,
     *            amunit=1.66053873e-27,bzcon=1.380603e-23,
     *            Avcon = 6.02214199e23)
        logical task, conv
        double precision rxx(100), ryy(100), rzz(100)
        double precision fp
        real*8 leng

!          leng = 1.0e10*runit

          OPEN  (02, FILE = 'input.car')
          WRITE (02, '(a)')  '!BIOSYM archive 3'
          WRITE (02, '(a)')  'PBC=OFF'
          WRITE (02, '(a)')  'Materials Studio Generated CAR File'
          WRITE (02, '(a)')  '!DATE Wed Nov 29 13:31:15 2006'

          do ia = 1, 20
            WRITE (02, '(A,2X,3f15.9,21X,A,f8.3)') 'Cu ',
     +             rxx(ia), ryy(ia), rzz(ia),'Cu',0.000
          enddo
          do ia = 21, 32 
            WRITE (02, '(A,2X,3f15.9,21X,A,f8.3)') 'C  ',
     +             rxx(ia), ryy(ia), rzz(ia),'C ',0.000
          enddo

          WRITE (02,'(a)')  'end'
          WRITE (02,'(a)')  'end'

          CLOSE (02)


          call system('./*.sh -np 8 input ')
          call system('cat input.outmol END > tmpp')
          call system('cp tmpp input.outmol')
          call system('rm tmpp')

          call EXTRACT(rxx,ryy,rzz,fp,conv)

        return
        end

