             subroutine quench_opt(elem,x,y,z,energy,conv)
             include "param.h"

             double precision x(100),y(100),z(100),energy
             character*2     elem(100)
             logical     conv

             call write_input(elem,x,y,z)

             call system('g03 exam1.gjf')

             call read_output(elem,x,y,z,energy,conv)

             return

             end subroutine
