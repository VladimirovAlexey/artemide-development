!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Simple test-code for Snowlake. Checks compilation integrity and simple evaluation        !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program SnowflakeTEST
use HexGrid

implicit none

real*8::t1,t2
!$ real*8::omp_get_wtime
real*8,dimension(0:815)::grid,Dgrid
real*8::x1,x2,x3,f1,f2
integer::i,j,n
logical::intersect

call  Initialize_HexGrid("Snowflake.ini")

call cpu_time(t1)
!$ t1=omp_get_wtime()


grid=GETgrid(initialF)
Dgrid=Dgrid_dX2(grid)

OPEN(UNIT=44, FILE=trim("prprpr.dat"), ACTION="write", STATUS="replace")

do i=0,size(grid)
    call get_X123_from_1Dindex(i,x1,x2,x3)
    write(44,*) x1,x2,grid(i),Dgrid(i)

end do

!
! do i=-50,50
! do j=-50,50
!
! if(i==0 .and. j==0) cycle
! x1=i/51.d0
! x2=j/51.d0
! if(abs(x1+x2)>1.) cycle
! f1=GETinterpolation(x1,x2,grid)
! f2=GETinterpolation(x1,x2,Dgrid)
!
! write(44,*) x1,x2,f1,f2
! end do
! end do

call cpu_time(t2)
!$ t2=omp_get_wtime()
write(*,*) "Time for computation of evolution",t2-t1

write(*,*) "  "
CLOSE(44, STATUS='KEEP')
contains

!!!! some symmetric function
    function initialF(x,y)
        real*8,intent(in)::x,y
        real*8::initialF

        if(abs(x)<=1 .and. abs(y)<=1 .and. abs(x+y)<=1) then
            initialF=(1-x**2)*(1-y**2)*(1-(x+y)**2)*(1+3*x+4*y+0.2*x*(x+y))
        else
            initialF=0.d0
        end if
    end function initialF

end program SnowflakeTEST
