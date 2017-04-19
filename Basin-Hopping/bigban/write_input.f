                subroutine write_input(elem,x,y,z)
                include "param.h"

                double precision x(100),y(100),z(100)
                character*2     elem(100)

                open(unit=90, file='exam1.gjf', status='unknown')

                write(90,*) '%mem=120mw'
                write (90, *) '%nproc=4'
                write (90, *) '#p opt hf/sto-3g'
                write (90, *) ' scf(conver=6,maxcyc=300) nosymm'
                write (90, *)
                write (90, *) 'test'
                write (90, *)
                write (90, *) '-2 1'
                do i=1, Natoms
                   write(90, *) elem(i), x(i), y(i), z(i)
                end do
                write (90, *)


                close (90)


                return
                end




