!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Simple test-code for Snowlake. Checks compilation integrity and simple evaluation        !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program SnowflakeTEST
use HexGrid
use EvolutionKernels
use SnowFlake

implicit none

integer::i,j,n,k
real*8::x1,x2,r1,r2,xx1,xx2,t1,t2

call  SnowFlake_Initialize("TEST.ini","prog/")
n=14
k=0

call cpu_time(t1)

r1=G2_projM(n,k,0.2d0)
write(*,*) r1
r1=G2_projM(n,k,0.22d0)
write(*,*) r1
r1=G2_projM(n,k,0.24d0)
write(*,*) r1


!call NKtoX12(n,k,x1,x2)

!do i=-25,25
!do j=-25,25
!    xx1=x1+i*0.0003d0
!    xx2=x2+j*0.0003d0
!    r1=Windex(n,k,xx1,xx2)
!    r2=Windex_d3_f1(n,k,xx1,xx2)
!    write(*,'("{",F12.8,", ",F12.8,", ",F12.8,", ",F16.8,"},")') xx1,xx2,r1,r2
!end do
!end do


call cpu_time(t2)

write(*,*) ">>>>>>>> Timing: ", t2-t1

end program SnowflakeTEST
