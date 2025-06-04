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
real*8::grid(0:3119),x1,x2,x3,f1,f2
integer::i,j,n
logical::intersect

call  Initialize_HexGrid("Snowflake.ini")

call cpu_time(t1)
!$ t1=omp_get_wtime()
!call ComputeEvolution(1.d0,25d0,alpha,U1=initialF,U2=initialA,G1=initialG,inputQ="T",inputG="T")

!grid=GETgrid(initialF)
! !
! ! x1=GETinterpolation(0.3d0,-0.2d0,grid)
! !
! ! write(*,*) x1,initialF(0.3d0,-0.2d0)
!
! call get_X123_from_1Dindex(536,x1,x2,x3)
! write(*,*) x1,x2,x3
!
! x3=GETinterpolatorB(536,x1,x2)
!
! write(*,*) x3
!
! stop
!
! x3=GETinterpolatorB(530,0.3d0,0.2d0)
!
! write(*,*) x3
!
! stop

n=431
call get_X123_from_1Dindex(n,x1,x2,x3)
write(*,*) n, "-->",x1,x2,x3

n=385
call get_X123_from_1Dindex(n,x1,x2,x3)
write(*,*) n, "-->",x1,x2,x3


call LimitsX2(431,x1,x2,intersect,f1,f2)
write(*,*) "--------->",intersect,f1,f2

call LimitsX2(385,x1,x2,intersect,f1,f2)
write(*,*) "--------->",intersect,f1,f2
stop



!!call X123_fromRP(0.045d0,-0.1d0,x1,x2,x3)
!n=340
!call get_X123_from_1Dindex(n,x1,x2,x3)
!write(*,*) "PHI---->",x1,x2,x3
!x3=GETinterpolatorB(n,x1,x2)
!write(*,*) "---->",x1,x2,x3


! !call X123_fromRP(0.2d0,phi,x1,x2,x3)
!
do n=0,100
 !call X123_fromRP(0.0107502d0,3+2*n/1000.d0,x1,x2,x3)
 !call get_X123_from_1Dindex(n,x1,x2,x3)

 x3=x1+(1-x2-x1)*n/100
 write(*,'("{",F12.8,",",F12.8,"},")') n/100.d0,GETinterpolatorB(385,x3,x2)

 !x3=x1+(1-x1)*n/100
 !write(*,'("{",F12.8,",",F12.8,"},")') n/100.d0,GETinterpolatorB(431,x3,x2)
end do

stop

!!!!! sum of interpolators over full grid at any point should be 1
x3=0.d0
do n=0,3119
x3=x3+GETinterpolatorB(n,0.1d0,0.1d0)
end do

write(*,*) "---->",x3
stop

OPEN(UNIT=44, FILE=trim("prprpr.dat"), ACTION="write", STATUS="replace")

do i=-50,50
do j=-50,50

if(i==0 .and. j==0) cycle
x1=i/51.d0
x2=j/51.d0
if(abs(x1+x2)>1.) cycle
f1=GETinterpolation(x1,x2,grid)
f2=0.d0
do n=0,size(grid)
    f2=f2+grid(n)*GETinterpolatorB(n,x1,x2)
end do
!x3=GETinterpolatorB(n,x1,x2)
write(44,*) x1,x2,f1,f2
!write(*,'("{",F12.6,",",F12.6,",",F12.8,"},")') x1,x2,x3
! write(*,'("{",F12.6,",",F12.6,",",F12.8,",",F12.6,",",F12.8,"},")') &
!     x1,x2!(x3-initialF(x1,x2))/(abs(initialF(x1,x2))+0.0000000001),x3,initialF(x1,x2)
end do
end do

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

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            initialF=(1-x**2)*(1-y**2)*(1-(x+y)**2)!*(1+10*log(abs(x)+2.)-4*y)/max(abs(x),abs(y),abs(x+y))*Exp(-x+y)
        else
            initialF=0.d0
        end if
    end function initialF

end program SnowflakeTEST
